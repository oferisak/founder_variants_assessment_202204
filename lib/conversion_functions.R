chr_refseq2ucsc<-read.table('/media/SSD/Bioinformatics/Databases/refseq/chromosomes_refseq_to_ucsc_20211101.txt.gz',sep='\t',header=T,quote = '',comment.char = '')
# refseq to gene
refseq2gene<-read.table('/media/SSD/Bioinformatics/Databases/refseq/refseq_curated_to_gene.20231004.csv.gz',sep='\t',header=T,quote = '',comment.char = '')
colnames(refseq2gene)<-c('X.name','chrom','name2')
refseq2gene<-refseq2gene%>%
  filter(!grepl('_',chrom))%>%
  mutate(original_transcript_without_isoform=stringr::str_extract(X.name,'N[MRX]_[0-9]+'))

mutalyzer_cds_converter <- function(transcript_name,cds_change) {
  transcript_found<-F
  isoform_num<-0
  transcript_name_without_isoform<-transcript_name%>%str_replace('\\.\\d+','')
  if(nrow(refseq2gene%>%filter(grepl(transcript_name_without_isoform,X.name)))==0){stop(sprintf('The provided transcript name %s in not found in the refseq list',transcript_name))}
  
  gene_name<-unique(refseq2gene%>%filter(grepl(transcript_name_without_isoform,X.name))%>%pull(name2))
  other_refseq_names<-refseq2gene%>%filter(name2==gene_name)%>%pull(X.name)
  stop<-F
  transcript_options<-c(transcript_name,setdiff(other_refseq_names,transcript_name))
  for (accession in transcript_options){
    for(isoform_num in 1:10){
      message(glue('Trying to retrieve coordinates for {accession}:{cds_change}'))
      mutalyzer_query<-sprintf("https://mutalyzer.nl/position-converter?assembly_name_or_alias=GRCh37&description=%s:%s",
                               accession,
                               cds_change)
      mutalyzer_query<-read_html(URLencode(mutalyzer_query))
      error<-mutalyzer_query%>%html_nodes('p.alert.alert-danger')%>%html_text()
      #print(error)
      if (length(error)==0){
        transcript_found<-T
        stop<-T
        break
      }else if (grepl('version',error)){# the isoform is not found in mutalayzer, try a different isoform
        accession<-accession%>%str_replace('\\.\\d+',sprintf('.%s',isoform_num))
        message(sprintf('isoform not found, will try: %s',accession))
      }else{
        message('accession not found in mutalyzer, will try a different one')
        break
      }
    }
    if (stop){break}
  }
  if (!transcript_found){stop(sprintf('could not find transcript: %s\n',transcript_name))}
  genomic_coordinates<-mutalyzer_query%>%html_nodes("code:not(.example-input)")%>%html_text()
  refseq_coords<-mutalyzer_query%>%html_nodes("pre")%>%html_text2()%>%str_replace_all('\\s+','|')%>%str_match_all('NM_[^\\|]+')%>%unlist()
  output<-data.frame(type='genomic',coordinates=genomic_coordinates)%>%
    bind_rows(data.frame(type='cds',coordinates=refseq_coords))%>%
    mutate(coords2split=coordinates)%>%
    separate(coords2split,c('transcript_id','change'),sep = ':')%>%
    left_join(chr_refseq2ucsc%>%dplyr::select(X.chrom,name),by=c('transcript_id'='name'))
  
  return(output)
}

mutalyzer_to_bed<-function(mutalyzer_df,bed_prefix='mutalyzer_to_bed_file'){
  original_size<-nrow(mutalyzer_df)
  mutalyzer_df<-mutalyzer_df%>%
    mutate(grch38_position=as.numeric(stringr::str_replace_all(g_GRCH38,'^[^:]+:g.','')%>%str_extract('\\d+')))
  bed_file<-mutalyzer_df%>%
    filter(!is.na(g_GRCH38))%>%
    select(chrom,start=grch38_position,end=grch38_position,base_var,gene_name)
  bed_size<-nrow(bed_file)
  output_file<-glue('./output/{bed_prefix}.bed')
  message(glue('Converted {bed_size}/{original_size} variants.. writing into {output_file}'))
  write.table(bed_file,file = output_file,quote = F,sep = '\t',row.names = F,col.names = F)
  return(output_file)
}

vep_to_bed<-function(vep_df,bed_prefix='mutalyzer_to_bed_file',add_chr=F){
  original_size<-nrow(vep_df)
  bed_file<-vep_df%>%
    filter(!is.na(start))%>%
    select(chr,start,end,gene,variant,hgvs_notation)%>%
    distinct()
  # fix it so that the start is always before the end
  fixed_bed_file<-NULL
  for (i in 1:nrow(bed_file)){
    bed_line<-bed_file%>%slice(i)
    new_bed_line<-bed_line
    if (bed_line$start>bed_line$end){
      new_bed_line<-new_bed_line%>%mutate(start=bed_line$end,end=bed_line$start)
    }
    fixed_bed_file<-fixed_bed_file%>%bind_rows(new_bed_line)
  }
  bed_file<-fixed_bed_file
  if (add_chr){
    bed_file<-bed_file%>%mutate(chr=glue('chr{chr}'))
  }
  bed_size<-nrow(bed_file)
  output_file<-glue('./output/{bed_prefix}.bed')
  message(glue('Converted {bed_size}/{original_size} variants.. writing into {output_file}'))
  write.table(bed_file,file = output_file,quote = F,sep = '\t',row.names = F,col.names = F)
  return(output_file)
}
