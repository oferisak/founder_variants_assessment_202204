chr_accessions<-readr::read_delim('./data/accessory_data/chr_accession.csv')

send_variant_to_mutalyzer_normalize<-function(variant){
  message(glue('Sending {variant} to mutalyzer..'))
  response_status=500
  if (grepl(' ',variant_for_mutalyzer)){
    response=list()
    response$status_code<-404
    response$error_details<-'URL using bad/illegal format or missing URL'
    return(response)
  }
  while(response_status==500){
    response <- GET(glue('https://v3.mutalyzer.nl/api/normalize/{variant}?only_variants=false'),timeout(120))
    response_status<-response$status_code
    if (response_status==500){message('Response status is 500 (internal server error) will retry..')}
  }
  if (response_status!=200){
    if (response_status==404){
      message('404 not found')
      response$error_details<-'404 not found'
    }else{
      char <- rawToChar(response$content)
      error_details<-jsonlite::fromJSON(char)
      error_details<-error_details$custom$errors$details
      message(error_details)
      response$error_details<-error_details
    }
  }
  return(response)
}

# a function that takes in the response from the mutalyzer variant normalization api and returns whether mutalyzer reports a difference
# between the reference and the variant sequence
is_difference_between_reference_and_variant<-function(response){
  char <- rawToChar(response$content)
  #Convert to df 
  df <- jsonlite::fromJSON(char)
  # grab the details
  details<-df$infos$details
  if (sum(grepl('There are differences',details))>0){return(TRUE)}
  #if (sum(grepl('No.+selector found',details))>0){return(TRUE)}
  return(FALSE)
}

parse_mutalyzer_normalize_response<-function(response,is_intronic=F){
  # Convert response to dataframe
  char <- rawToChar(response$content)
  #Convert to df 
  df <- jsonlite::fromJSON(char)
  norm_hgvs<-df$normalized_description
  
  if (is_intronic){
    genomic_locations=data.frame(g_GRCH38=df$equivalent_descriptions$g)
    protein<-NA
  }else{
    if (is.null(df$chromosomal_descriptions)){
      message('FAIL: Could not find chromosomal locations for variant..')
      genomic_locations=data.frame(c_GRCH37=NA,g_GRCH38=NA)
      protein<-NA
    }else{
      genomic_locations<-df$chromosomal_descriptions%>%pivot_wider(names_from = assembly,values_from = c(c,g))
      protein<-df$protein$description
    }
  }
  message(glue('Parsed input variant. genomic location: {paste0(genomic_locations,collapse=" ")}'))
  return(data.frame(norm_hgvs,genomic_locations,protein))
}

parse_mutalyzer_normalize_response_with_chr_location<-function(response){
  # Convert response to dataframe
  char <- rawToChar(response$content)
  #Convert to df 
  df <- jsonlite::fromJSON(char)
  norm_hgvs<-df$normalized_description
  genomic_locations=data.frame(g_GRCH38=df$equivalent_descriptions$g)
  message(glue('Parsed input variant. genomic location: {paste0(genomic_locations,collapse=" ")}'))
  return(data.frame(norm_hgvs,genomic_locations))
}
