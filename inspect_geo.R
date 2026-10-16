base<-'E:/NHANES/analysis/bioinfo/data/geo'
for(f in c('GSE98895_series_matrix.txt.gz','GSE9006-GPL96_series_matrix.txt.gz')) {
 con<-gzfile(file.path(base,f));s<-readLines(con);close(con)
 h<-s[grepl('!Sample_',s,fixed=TRUE)]
 writeLines(h,file.path('work',paste0(f,'_headers.txt')))
 cat('\n',f,'\n');print(substr(s[grepl('!Sample_data_processing',s,fixed=TRUE)],1,500))
 print(substr(s[grepl('!Sample_characteristics_ch1',s,fixed=TRUE)],1,400))
 idx<-grep('!series_matrix_table_begin',s,fixed=TRUE);print(substr(s[(idx+1):(idx+3)],1,400))
}
