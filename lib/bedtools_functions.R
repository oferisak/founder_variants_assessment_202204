# intersects two files, and returns all the rows from a_file with the number of overlapping features from b_file (the -c argument)
bedtools_intersect<-function(a_file,b_file){
  intersect_command<-glue('bedtools intersect -c -a {a_file} -b {b_file}')
  message(glue('Running {intersect_command}'))
  intersect_output_raw<-system(intersect_command,intern = TRUE)
  intersect_output<-read.table(text=intersect_output_raw,sep='\t',header = F)
  return(intersect_output)
}
