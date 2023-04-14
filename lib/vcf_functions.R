vcf_to_tidy<-function(input_vcfR,happy_individual=NA){
  message(glue('Converting VCF to tidy form..'))
  vcf_tidy<-vcfR2tidy(input_vcfR)
  gt_df<-vcf_tidy$gt
  if (!is.na(happy_individual)){
    gt_df<-gt_df%>%filter(Indiv==happy_individual)
  }
  info_df<-vcf_tidy$fix
  if (nrow(info_df)!=nrow(gt_df)){stop(glue('The number of rows in the fix part of the vcf is different from the number in the gt part'))}
  vcf_df<-info_df%>%
    bind_cols(gt_df%>%select(-c(ChromKey,POS,Indiv)))
  return(vcf_df)
}
