chr_refseq2ucsc<-read.table('/media/SSD/Bioinformatics/Databases/refseq/chromosomes_refseq_to_ucsc_20211101.txt.gz',sep='\t',header=T,quote = '',comment.char = '')
# refseq to gene
refseq2gene<-read.table('/media/SSD/Bioinformatics/Databases/refseq/refseq_to_gene_20211101.txt.gz',sep='\t',header=T,quote = '',comment.char = '')

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