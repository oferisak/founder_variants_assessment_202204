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

vep_cds_to_genomic_coordinates<-function(hgvs_notation,build='grch37'){
  if (build=='grch37'){
    server <- "https://grch37.rest.ensembl.org"
  }else{
    server<-"https://rest.ensembl.org"
  }
  vep_hgvs_ext <- "/vep/human/hgvs/:hgvs_notation"
  message(glue('Parsing {hgvs_notation} with build {build}'))
  ext <- glue("/vep/human/hgvs/{hgvs_notation}?refseq=true&canonical=true&ambiguous_hgvs=1&pick=1")
  
  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  char <- rawToChar(r$content)
  
  res<-jsonlite::fromJSON(char)

  if ('error' %in% names(res)){
    res<-data.frame(hgvs_notation=hgvs_notation,error=res$error)
    message(glue('VEP: {res$error}'))
    return(res)
  }
  transcript_cons<-bind_rows(res$transcript_consequences)
  res<-res%>%
    select(hgvs_notation=input,chr=seq_region_name,start,end,allele_string,strand,assembly_name)%>%
    separate(allele_string,into=c('ref','alt'),sep='/')%>%
    bind_cols(transcript_cons%>%
                select(-c(gene_symbol_source,biotype,gene_id,impact,strand,variant_allele)))
    
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
