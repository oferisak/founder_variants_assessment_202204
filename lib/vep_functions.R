vep_server <- "https://rest.ensembl.org"
vep_ext <- "/vep/human/hgvs"


vep_annotate_variants<-function(hgvs_notations){
  message(glue('VEP: submitting {length(hgvs_notations)} variants to VEP annotation..'))
  response <- POST(paste(vep_server, vep_ext, sep = ""), 
                   #query = list('refseq'=1,'hgvs'=1,'canonical'=1),
                   #content_type("application/json"),
                   content_type_json(),
                   accept_json(),
                   #accept("application/json"), 
                   body = toJSON(list('hgvs_notations' = hgvs_notations)))
  #body = glue('{ "hgvs_notations" : {{hgvs_notations}} }',.open = '{{',.close = '}}'))
  
  stop_for_status(response)
  message('VEP: Retrieved without errors.')
  char <- rawToChar(response$content)
  res<-jsonlite::fromJSON(char)
  annotated_vars<-data.frame(hgvs_notation=res$id,
                             assembly_name=res$assembly_name,
                             seq_region_name=res$seq_region_name,
                             start=res$start,end=res$end,
                             allele_string=res$allele_string,
                             strand=res$strand)
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
    ungroup()
  message(glue('Out of {nrow(verified_annotated_vars)} variants, {nrow(verified_annotated_vars%>%filter(verified))} were validated'))
  if (!keep){verified_annotated_vars<-verified_annotated_vars%>%filter(verified)}
  return(verified_annotated_vars)
}
