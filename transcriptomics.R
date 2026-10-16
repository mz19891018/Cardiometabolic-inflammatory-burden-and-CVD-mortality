.libPaths(c('work/Rlib',.libPaths()));library(limma);library(data.table);library(jsonlite)
base<-'E:/NHANES/analysis/bioinfo/data/geo';out<-'outputs/final_analysis';dir.create(out,recursive=TRUE,showWarnings=FALSE)
read_geo<-function(f){s<-readLines(gzfile(file.path(base,f)));begin<-grep('!series_matrix_table_begin',s,fixed=TRUE);end<-grep('!series_matrix_table_end',s,fixed=TRUE)
 y<-as.matrix(read.delim(text=paste(s[(begin+1):(end-1)],collapse='\n'),row.names=1,check.names=FALSE));storage.mode(y)<-'double'
 vals<-function(line)as.character(read.delim(text=line,header=FALSE,check.names=FALSE)[1,-1])
 titles<-vals(s[startsWith(s,'!Sample_title')]);gsm<-vals(s[startsWith(s,'!Sample_geo_accession')]);meta<-data.frame(sample=gsm,title=titles)
 for(line in s[startsWith(s,'!Sample_characteristics_ch1')]){v<-vals(line);label<-trimws(strsplit(v[which(nzchar(v))[1]],':',fixed=TRUE)[[1]][1]);if(nzchar(label))meta[[tolower(label)]]<-sub('^[^:]*: *','',v)}
 list(y=y,meta=meta,processing=vals(s[startsWith(s,'!Sample_data_processing')])[1])}
annot<-function(f){s<-readLines(gzfile(file.path(base,f)));beg<-which(startsWith(s,'ID\t'))[1];a<-read.delim(text=paste(s[beg:length(s)],collapse='\n'),check.names=FALSE,quote='',fill=TRUE);a[,c('ID','Gene symbol')]}
allsummary<-list();alleffect<-list()
for(gse in c('GSE98895','GSE9006')){
 z<-read_geo(if(gse=='GSE98895')'GSE98895_series_matrix.txt.gz'else'GSE9006-GPL96_series_matrix.txt.gz');y<-z$y;md<-z$meta
 if(gse=='GSE98895'){md$group<-ifelse(md[['disease state']]=='Metabolic Syndrome','MetS','Control');md$age<-as.numeric(md$age);keep<-rep(TRUE,nrow(md));anfile<-'GPL6947.annot.gz'}else{md$group<-ifelse(md$illness=='Healthy','Control',ifelse(md$illness=='Type 2 Diabetes','T2D',ifelse(endsWith(md$title,'T1D_New'),'T1D','Followup')));keep<-md$group!='Followup';md$age<-as.numeric(sub(' .*','',md$age));md$sex<-factor(md$gender);anfile<-'GPL96.annot.gz'}
 md<-md[keep,];y<-y[,md$sample,drop=FALSE];q<-quantile(y,c(0,.01,.25,.5,.75,.99,1),na.rm=TRUE)
 stopifnot(q[6]>100); y<-log2(pmax(y,1));# Both supplied matrices contain linear intensities, metadata checked.
 norm<-if(gse=='GSE98895')'Provided GenomeStudio quantile normalization; log2 only'else'MAS5 scaled intensities; log2 then quantile normalization of baseline arrays'
 if(gse=='GSE9006')y<-normalizeBetweenArrays(y,method='quantile')
 a<-annot(anfile);idx<-match(rownames(y),a$ID);gene<-a[['Gene symbol']][idx];valid<-!is.na(gene)&nzchar(gene)&!grepl('///',gene,fixed=TRUE);y<-y[valid,,drop=FALSE];gene<-gene[valid]
 # Label-independent representative probe: highest IQR. Filter low signal consistently across all baseline samples.
 iq<-apply(y,1,IQR,na.rm=TRUE);ord<-order(iq,decreasing=TRUE);take<-ord[!duplicated(gene[ord])];y<-y[take,,drop=FALSE];rownames(y)<-gene[take]
 if(gse=='GSE9006')y<-y[rowMeans(y>log2(50))>=.2,,drop=FALSE]
 y<-y[apply(y,1,sd)>0,,drop=FALSE]
 md$group<-factor(md$group,levels=if(gse=='GSE98895')c('Control','MetS')else c('Control','T1D','T2D'))
 design<-model.matrix(if(gse=='GSE98895')~0+group+age else ~0+group+age+sex,md)
 contr<-makeContrasts(contrasts=if(gse=='GSE98895')'groupMetS-groupControl'else c('groupT1D-groupControl','groupT2D-groupControl'),levels=design)
 fit0<-lmFit(y,design);fit<-eBayes(contrasts.fit(fit0,contr),trend=TRUE,robust=TRUE)
 write.csv(md,file.path(out,paste0(gse,'_metadata.csv')),row.names=FALSE)
 png(file.path(out,paste0(gse,'_QC.png')),width=1800,height=1200,res=170);par(mfrow=c(2,2),mar=c(4,4,3,1));boxplot(y,outline=FALSE,axes=FALSE,main=paste(gse,'log2 expression'));axis(2);pca<-prcomp(t(y),scale.=FALSE);plot(pca$x[,1:2],col=as.integer(md$group),pch=19,xlab='PC1',ylab='PC2',main='PCA by clinical group');legend('topright',levels(md$group),col=seq_along(levels(md$group)),pch=19,cex=.7);plot(md$age,pca$x[,1],pch=19,col=as.integer(md$group),xlab='Age',ylab='PC1',main='Age and expression');plotMDS(y,col=as.integer(md$group),pch=19,main='MDS');dev.off()
 for(k in seq_len(ncol(contr))){nm<-paste0(gse,'_',if(gse=='GSE98895')'MetS'else c('T1D','T2D')[k]);tt<-topTable(fit,coef=k,number=Inf,sort.by='P');tt$gene<-rownames(tt);write.csv(tt,file.path(out,paste0(nm,'_limma.csv')),row.names=FALSE)
  png(file.path(out,paste0(nm,'_volcano.png')),width=1500,height=1050,res=170);sig<-tt$adj.P.Val<.05&abs(tt$logFC)>=.5;plot(tt$logFC,-log10(pmax(tt$adj.P.Val,1e-300)),pch=16,cex=.5,col=ifelse(sig,'#bc3a40','#8c9caa'),xlab='Adjusted log2 fold change',ylab='-log10 FDR',main=paste(nm,'limma'));abline(v=c(-.5,.5),h=-log10(.05),lty=2,col='gray50');sel<-head(which(sig),8);text(tt$logFC[sel],-log10(tt$adj.P.Val[sel]),tt$gene[sel],pos=3,cex=.7);dev.off()
  if(file.exists('work/hallmark.gmt')){gmt<-strsplit(readLines('work/hallmark.gmt'),'\t');sets<-lapply(gmt,function(x)x[-c(1,2)]);names(sets)<-sapply(gmt,'[',1);ind<-ids2indices(sets,rownames(y));cam<-camera(y,ind,design,contrast=contr[,k],inter.gene.cor=.01);cam$pathway<-rownames(cam);write.csv(cam,file.path(out,paste0(nm,'_hallmark_camera.csv')),row.names=FALSE)}
  allsummary[[nm]]<-data.frame(comparison=nm,n=ncol(y),n_control=sum(md$group=='Control'),n_case=sum(md$group==if(gse=='GSE98895')'MetS'else c('T1D','T2D')[k]),genes=nrow(y),FDR05=sum(tt$adj.P.Val<.05),FDR05_FC05=sum(sig),up=sum(sig&tt$logFC>0),down=sum(sig&tt$logFC<0),normalization=norm,adjustment=if(gse=='GSE98895')'age'else'age and sex')
  alleffect[[nm]]<-tt[tt$gene%in%c('IL1B','PTGS2','TNF','IL6','CPT1A','SCD','SOD2','EGR2','EGR3','FOSB'),c('gene','logFC','P.Value','adj.P.Val')];alleffect[[nm]]$comparison<-nm
 }
 # Unadjusted MetS sensitivity quantifies age confounding; it is not selected by significance.
 if(gse=='GSE98895'){de0<-model.matrix(~md$group);f0<-eBayes(lmFit(y,de0),trend=TRUE,robust=TRUE);tt0<-topTable(f0,coef=2,number=Inf);tt0$gene<-rownames(tt0);write.csv(tt0,file.path(out,'GSE98895_unadjusted_sensitivity.csv'),row.names=FALSE)}
 saveRDS(list(y=y,metadata=md,design=design,fit=fit),file.path('work',paste0(gse,'_processed.rds')))
 writeLines(c(z$processing,norm,paste(names(q),round(q,3),collapse='; ')),file.path(out,paste0(gse,'_preprocessing.txt')))
}
write.csv(do.call(rbind,allsummary),file.path(out,'transcriptome_summary.csv'),row.names=FALSE);write.csv(do.call(rbind,alleffect),file.path(out,'candidate_gene_results.csv'),row.names=FALSE)
