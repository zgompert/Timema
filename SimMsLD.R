## mel-stripe allele frequency
p_ms<-0.32 ## alt allele, q = 0.68

## read overall population allele frequency estiamtes
p<-scan("../data/est_pem_tcr_ecol_all_1X.txt")
p<-sample(p,7000000,replace=FALSE)## just use 7 million of the SNPs
L<-length(p) 


LDmat<-matrix(NA,nrow=5,ncol=10)## rows different Ne, cols reps
Ne<-c(50,100,200,500,1000)
for(k in 1:5){
	for(j in 1:10){
		Gms<-rbinom(n=Ne[k],size=2,prob=p_ms)
		G<-matrix(rbinom(n=Ne[k]*L,size=2,p),nrow=L,ncol=Ne[k])
		ld<-cor(Gms,t(G))^2
		LDmat[k,j]<-sum(ld[1,] >0.1,na.rm=TRUE)
	}
}

save(list=ls(),file="LDsims.rdat")
library(scales)
pdf("fig_sim.pdf",width=5,height=5)
par(mar=c(5,5,1,1))
boxplot(t(log10(LDmat)),ylim=c(0,5.5),names=Ne,xlab="Effective population size",ylab="Log number of high-LD loci",cex.lab=1.5,col=alpha("cadetblue",.5))
for(i in 1:5){

	if(i<4){points(jitter(rep(i,10),factor=2),log10(LDmat[i,]),pch=21,bg=alpha("cadetblue",.9))}
	if(i>3){points(jitter(rep(i,10),factor=2),LDmat[i,],pch=21,bg=alpha("cadetblue",.9))}
}

dev.off()

