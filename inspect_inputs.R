.libPaths(c('work/Rlib',.libPaths())); Sys.setlocale('LC_CTYPE','Chinese_China.utf8'); library(foreign)
paths<-c(list.files('E:/NHANES/analysis/bioinfo/data/gwas',pattern='gz$',full.names=TRUE),list.files('E:/NHANES/analysis/bioinfo/gwas',pattern='gz$',full.names=TRUE))
for(p in paths){cat('\n',p,'\n');con<-gzfile(p);print(readLines(con,n=2));close(con)}
for(s in c('1999-2000','2001-2002','2003-2004')){dr<-paste0('E:/NHANES/',s,'年');for(f in list.files(dr,pattern='xpt$',full.names=TRUE)){if(any(startsWith(basename(f),c('LAB','L10','L11','L13','L40','DEMO','BPQ')))){x<-read.xport(f);cat(s,basename(f),':',paste(names(x),collapse=','),'\n')}}}
