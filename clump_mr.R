.libPaths(c('work/Rlib',.libPaths()));library(data.table)
base<-'work/ld/EUR';plink<-normalizePath('work/plink/plink.exe',winslash='/');files<-list.files('work/mr',pattern='_candidates.csv$',full.names=TRUE);all<-list();flow<-list()
for(f in files){x<-fread(f);x<-x[order(p)][!duplicated(SNP)];nm<-sub('_candidates.csv','',basename(f),fixed=TRUE);p<-file.path('work/mr',paste0(nm,'_clump_input.tsv'));fwrite(x[,.(SNP,P=p)],p,sep='\t');prefix<-file.path('work/mr',paste0(nm,'_LD'))
 system2(plink,c('--bfile',base,'--clump',p,'--clump-p1','5e-8','--clump-p2','1','--clump-r2','0.001','--clump-kb','10000','--memory','2048','--threads','2','--out',prefix),stdout=FALSE,stderr=FALSE)
 z<-fread(paste0(prefix,'.clumped'));kept<-x[SNP%in%z$SNP];fwrite(kept,file.path('work/mr',paste0(nm,'_instruments.csv')));all[[nm]]<-kept$SNP;flow[[nm]]<-data.frame(exposure=nm,candidates=nrow(x),LD_independent=nrow(kept));cat(nm,nrow(x),'->',nrow(kept),'\n')
}
writeLines(unique(unlist(all)),'work/mr/query_snps.txt');write.csv(do.call(rbind,flow),'outputs/final_analysis/mr_instrument_flow.csv',row.names=FALSE)
# Joint lipid instrument selection uses the smallest P across all three exposures, followed by a single LD clump.
lip<-rbindlist(lapply(c('HDL','TG','LDL'),function(n)fread(file.path('work/mr',paste0(n,'_candidates.csv')))));lip<-lip[order(p)][!duplicated(SNP)];fwrite(lip[,.(SNP,P=p)],'work/mr/lipid_union_input.tsv',sep='\t')
system2(plink,c('--bfile',base,'--clump','work/mr/lipid_union_input.tsv','--clump-p1','5e-8','--clump-p2','1','--clump-r2','0.001','--clump-kb','10000','--memory','2048','--out','work/mr/lipid_union'),stdout=FALSE,stderr=FALSE)
z<-fread('work/mr/lipid_union.clumped');writeLines(z$SNP,'work/mr/lipid_union_snps.txt');writeLines(unique(c(unlist(all),z$SNP)),'work/mr/query_snps.txt')
