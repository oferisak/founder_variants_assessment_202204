vep_server <- "https://rest.ensembl.org"
vep_hgvs_ext <- "/vep/human/hgvs"

# coordinates should be in the format: chrX:start..end
vep_convert_coordinates<-function(assemb_from,assemb_to,coordinates){
  ext <- glue("/map/human/{assemb_from}/{coordinates}/{assemb_to}?")
  response <- GET(paste(vep_server, ext, sep = ""))
  char <- rawToChar(response$content)
  res<-jsonlite::fromJSON(char)
  col_names <- paste0(c("chr_", "start_", "end_"), assemb_to)
  ret_error <- setNames(data.frame(matrix(NA, nrow = 1, ncol = 3)), col_names)
  
  if ('error' %in% names(res)){return(ret_error)}
  return(res$mappings$mapped%>%select(!!sym(glue('chr_{assemb_to}')):=seq_region_name,
                               !!sym(glue('start_{assemb_to}')):=start,
                                     !!sym(glue('end_{assemb_to}')):=end))
}

find_bad_request<-function(hgvs_notations){
  for (hn in hgvs_notations){
    message(hn)
    annotations<-list('hgvs_notations' = hn,
                      'refseq'=1,
                      'per_gene'=1,
                      'canonical'=1)
    json_body<-toJSON(annotations)
    response <- POST(paste(vep_server, vep_hgvs_ext, sep = ""), 
                     content_type_json(),
                     accept_json(),
                     body = json_body)
    
    stop_for_status(response)
    
  }
}

vep_cds_to_genomic_coordinates<-function(variant,build='grch37',distinct_genomic_effect=TRUE){
  gene_symbol_or_refseq<-ifelse(!is.na(variant$transcript) & variant$transcript!='','refseq','gene_symbol')
  variant_position<-vep_hgvs_notation(variant%>%pull(hgvs_notation),build,gene_symbol_or_refseq = gene_symbol_or_refseq)
  # if the transcript is not found, and there is a gene symbol, try using it instead
  if (('error' %in% colnames(variant_position)) && grepl('Could not get a Transcript object',variant_position$error) && variant$gene!=''){
    variant$hgvs_notation<-glue('{variant$gene}:{variant$variant}')
    gene_symbol_or_refseq<-'gene_symbol'
    variant_position<-vep_hgvs_notation(variant%>%pull(hgvs_notation),build='grch37',gene_symbol_or_refseq=gene_symbol_or_refseq)
  }
  # if the variant is an insertion and the error suggests ambiguity, try changing it to dup
  if (('error' %in% colnames(variant_position)) && grepl('(I|i)ns',variant$variant)){
      original_hgvs_notation<-variant$hgvs_notation
      variant$hgvs_notation<-stringr::str_replace(variant$hgvs_notation,'(I|i)ns','dup')
      message(glue('VEP: parsing failed, will retry with {variant$hgvs_notation} instead of {original_hgvs_notation}'))
      variant_position<-vep_hgvs_notation(variant%>%pull(hgvs_notation),build,gene_symbol_or_refseq = gene_symbol_or_refseq)
  }
  return(variant_position)
}

vep_hgvs_notation<-function(hgvs_notation,build='grch37',gene_symbol_or_refseq,distinct_genomic_effect=TRUE){
  if (build=='grch37'){
    server <- "https://grch37.rest.ensembl.org"
  }else{
    server<-"https://rest.ensembl.org"
  }
  message(glue('VEP: Parsing {hgvs_notation} in build {build} using {gene_symbol_or_refseq}'))
  vep_hgvs_ext <- "/vep/human/hgvs/:hgvs_notation"
  if (gene_symbol_or_refseq=='refseq'){
    ext <- glue("/vep/human/hgvs/{hgvs_notation}?refseq=true&canonical=true&ambiguous_hgvs=1&hgvs=1&pick_order=mane_plus_clinical")
  }else{
    ext <- glue("/vep/human/hgvs/{hgvs_notation}?refseq=true&canonical=true&ambiguous_hgvs=1&hgvs=1")
  }
  
  r <- tryCatch(GET(paste(server, ext, sep = ""), content_type("application/json")),
                error=function(e){as.character(glue('parsing error: {e}'))})
  # if the response is an error handle it
  if (class(r)!='response'){
    message(glue('VEP parsing error: {r}'))
    res<-data.frame(hgvs_notation=hgvs_notation,error=r)
    return(res)
  }
  char <- rawToChar(r$content)
  res<-jsonlite::fromJSON(char)

  if ('error' %in% names(res)){
    res<-data.frame(hgvs_notation=hgvs_notation,error=res$error)
    message(glue('VEP: {res$error}'))
    return(res)
  }
  transcript_cons<-bind_rows(res$transcript_consequences)#%>%distinct()
  transcript_cons$original_res_num<-rep(1:length(res$transcript_consequences),sapply(res$transcript_consequences,nrow))
  transcript_cons<-transcript_cons%>%mutate(transcript_id_no_isoform=stringr::str_replace(transcript_id,'\\..+',''))
  res$original_res_num<-1:nrow(res)
  if (gene_symbol_or_refseq=='gene_symbol'){
    
    input_gene_symbol<-stringr::str_replace(hgvs_notation,':.+','')
    # make sure you keep only the transcript cons related to the gene symbol used in the query
    transcript_cons<-transcript_cons%>%filter(gene_symbol==input_gene_symbol)
  }
  if (gene_symbol_or_refseq=='refseq'){
    input_transcript_id_no_isoform<-unglue::unglue_vec(hgvs_notation,'{transcript}.{isoform}:{cdna_change}')
    transcript_cons<-transcript_cons%>%filter(transcript_id_no_isoform==input_transcript_id_no_isoform)
    if (nrow(transcript_cons)==0){
      return(data.frame(hgvs_notation=hgvs_notation,error='Could not find the given transcript in the transcripts output by VEP'))
    }
  }
  # join the transcript cons with the original res
  res<-res%>%
    select(original_res_num,hgvs_notation=input,chr=seq_region_name,start,end,allele_string,strand,assembly_name)%>%
    separate(allele_string,into=c('ref','alt'),sep='/')%>%
    left_join(transcript_cons%>%
                select(-c(gene_symbol_source,biotype,gene_id,impact,strand,variant_allele)))%>%
    select(-original_res_num)%>%
    distinct()
  # fix cannonical
  if ('canonical'%in%colnames(res)){
    res$is_canonical<-ifelse(is.na(res$canonical),0,1)
  }else{
    res$is_canonical<-0
  }
  # match cds start
  if ('hgvsc' %in% colnames(res)){
    res<-res%>%
      mutate(original_cds_start=stringr::str_replace(hgvs_notation,'.+c\\.+','')%>%stringr::str_extract('\\d+'),
             annotation_cds_start=stringr::str_replace(hgvsc,'.+c\\.+','')%>%stringr::str_extract('\\d+'))%>%
      mutate(cds_start_match=ifelse(original_cds_start==annotation_cds_start,1,0),
             cds_start_match=ifelse(is.na(cds_start_match),0,cds_start_match))%>%
      slice_max(cds_start_match)
    # collect most common cds and protein changes
    res<-res%>%mutate(cds_change=stringr::str_replace(hgvsc,'[^:]+:',''))
    res<-res%>%left_join(
      res%>%count(cds_change)%>%rename(cds_count=n)
    )%>%mutate(is_most_common_cds_change=ifelse(cds_count==max(cds_count),1,0))
  }else{
    res<-res%>%
      mutate(is_most_common_cds_change=0)%>%
      distinct(chr,start,end,ref,alt,.keep_all = TRUE)
  }
  if ('hgvsp' %in% colnames(res)){
    res<-res%>%mutate(protein_change=stringr::str_replace(hgvsp,'[^:]+:',''))
    res<-res%>%left_join(
      res%>%count(protein_change)%>%rename(protein_count=n)
    )%>%
      mutate(is_most_common_protein_change=ifelse(protein_count==max(protein_count),1,0))
  }else{
    res$is_most_common_protein_change=0
  }
  # add genomic ref alt
  res<-res%>%mutate(
    genomic_ref=ifelse(strand==1,ref,chartr('ACGT','TGCA',ref)),
    genomic_alt=ifelse(strand==1,alt,chartr('ACGT','TGCA',alt))
  )
  # for results with the same chromosomal location select the ones that are most common cds change > most common protein change > is canonical
  if (distinct_genomic_effect){
    res <- res %>%
      group_by(chr, start, end, ref, alt) %>%
      arrange(desc(is_most_common_cds_change), desc(is_most_common_protein_change),desc(is_canonical)) %>%
      slice(1) %>%
      ungroup()
  }
  # keep only valid chromosomes
  res<-res%>%
    mutate(valid_chr=ifelse(grepl('HG|HS',chr),0,1))%>%
    slice_max(valid_chr)%>%
    select(-valid_chr)
  # compare the hgvs notation for ins/del/dup (which result in ambiguity) if they are given in the hgvs_notation make sure they match the ref alt
  if (max(nchar(res$ref))<500 & max(nchar(res$alt))<500 & grepl('ins|dup|del',unique(res$hgvs_notation))){
    res<-res%>%
      rowwise()%>%
      mutate(is_matching_del=ifelse(grepl(glue('del{ref}'),hgvs_notation),1,0),
             is_matching_ins=ifelse(grepl(glue('ins{alt}'),hgvs_notation),1,0),
             is_matching_dup=ifelse(grepl(glue('dup{alt}'),hgvs_notation),1,0))%>%
      ungroup()%>%
      slice_max(is_matching_del+is_matching_ins+is_matching_dup)
  } 
  # cleanup
  unnecessary_cols<-c('distance','polyphen_score','sift_score','polyphen_prediction','sift_prediction','cds_end','cds_start','flags','hgvs_offset','bam_edit','refseq_offset')
  for (uc in unnecessary_cols){
    if (uc %in% colnames(res)){
      res<-res%>%select(-all_of(uc))
    }
  }
  message(glue('VEP: found {nrow(res)} matching variants'))
  return(res)
}

vep_basic_annotation<-function(hgvs_notations){
  message(glue('VEP: submitting {length(hgvs_notations)} variants to VEP annotation..'))
  all_res<-NULL
  annotations<-list('hgvs_notations' = hgvs_notations,
                    'refseq'=1,
                    'per_gene'=1,
                    'canonical'=1)
  json_body<-toJSON(annotations)
  response <- POST(paste(vep_server, vep_hgvs_ext, sep = ""), 
                   content_type_json(),
                   accept_json(),
                   body = json_body)
  
  stop_for_status(response)
  message('VEP: Retrieved without errors.')
  char <- rawToChar(response$content)
  all_res<-jsonlite::fromJSON(char)
  message(glue('Got response for {nrow(res)}/{length(hgvs_notations)}'))
  return(all_res)
}

vep_annotate_variants<-function(hgvs_notations,dbnsfp='REVEL_score,CADD_phred,REVEL_rankscore,VEST4_score,VEST4_rankscore,transcript_match=1'){
  notations_transcripts<-stringr::str_extract(hgvs_notations,'N[MRX]_[0-9\\.]+')
  message(glue('VEP: submitting {length(hgvs_notations)} variants to VEP annotation..'))
  annotations<-list('hgvs_notations' = hgvs_notations,
                    'SpliceAI'=1,
                    'refseq'=1,
                    'per_gene'=1,
                    'canonical'=1)
  if (!is.na(dbnsfp)){
    annotations<-c(annotations,list('dbNSFP'=unbox(dbnsfp)))
  }
  json_body<-toJSON(annotations)
  response <- POST(paste(vep_server, vep_hgvs_ext, sep = ""), 
                   content_type_json(),
                   accept_json(),
                   body = json_body)
  
  stop_for_status(response)
  message('VEP: Retrieved without errors.')
  char <- rawToChar(response$content)
  res<-jsonlite::fromJSON(char)
  message(glue('Got response for {nrow(res)}/{length(hgvs_notations)}'))
  # a function that given the vep results retrieves the annotation only corresponding to the given HGVS transcript
  pull_transcript_annotations<-function(res){
    output_transcript_lines<-NULL
    for (i in 1:nrow(res)){
      transcript_consequences<-res[i,]$transcript_consequences[[1]]
      notation_transcript<-res[i,'id']%>%stringr::str_extract('N[MRX]_[0-9\\.]+')
      notation_transcript_consequences<-transcript_consequences%>%filter(transcript_id==notation_transcript)%>%
        mutate(id=res[i,'id'])
      output_transcript_lines<-output_transcript_lines%>%bind_rows(notation_transcript_consequences)
    }
    return(output_transcript_lines)
  }
  notation_transcript_consequences<-pull_transcript_annotations(res)
  annotated_vars<-data.frame(hgvs_notation=res$id,
                             transcript_id=stringr::str_extract(res$id,'N[MRX]_[0-9\\.]+'),
                             assembly_name=res$assembly_name,
                             seq_region_name=res$seq_region_name,
                             start=res$start,end=res$end,
                             allele_string=res$allele_string,
                             strand=res$strand,
                             most_severe_consequence=res$most_severe_consequence)
  annotated_vars<-annotated_vars%>%left_join(notation_transcript_consequences,by=c('hgvs_notation'='id'))
  # fix the start end to be start>end
  annotated_vars<-annotated_vars%>%mutate(original_start=start,
                                          start=ifelse(start>end,end,start),
                                          end=ifelse(original_start>end,original_start,end))%>%select(-original_start)
  return(annotated_vars)
}

# a function that compares the alleles between the hgvs notation (input variant) and the allele string output by VEP
validate_vep_variants<-function(annotated_vars,keep=T){
  verified_annotated_vars<-
    annotated_vars%>%
    mutate(allele_to_verify=stringr::str_extract(hgvs_notation,'(del|dup|ins).+')%>%str_replace('del|dup|ins',''))%>%
    rowwise()%>%
    mutate(verified=ifelse((is.na(allele_to_verify)|allele_to_verify==''),TRUE,grepl(allele_to_verify,allele_string)))%>%
    ungroup()%>%
    mutate(verified=ifelse(is.na(gene_id),FALSE,verified))
  
  message(glue('Out of {nrow(verified_annotated_vars)} variants, {nrow(verified_annotated_vars%>%filter(verified))} were validated'))
  if (!keep){verified_annotated_vars<-verified_annotated_vars%>%filter(verified)}
  return(verified_annotated_vars)
}
