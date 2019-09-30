setwd('/chimeras/')
#written by Anna Liznar in R

one=read.csv('chimera_piRNA_1.tab', header= F, sep= "\t")
one=one['V2']
one$rep_1=TRUE

two=read.csv('chimera_piRNA_2.tab', header= F, sep= "\t")
two=two['V2']
two$rep_2=TRUE

three=read.csv('chimera_piRNA_3.tab', header= F, sep= "\t")
three=three['V2']
three$rep_3=TRUE

#################################################

#################################################

o_tw=merge(three, two, by= "V2", all=TRUE)
all=merge(o_tw, one, by= "V2", all=TRUE)


all[is.na(all)]=FALSE
all=unique(all)
rownames(all)=all$V2
all=all[,-1]
mall=as.matrix(all)

library(eulerr)
fit_all <- euler(mall)
fit_all
fit_per_all<- fit_all$fitted.values*100/sum(fit_all$fitted.values) #percentages
fit_per_all


pdf(file='./venn_123_reps.pdf', width = 5, height=3)
plot(fit_all,  alpha=0.7,counts=TRUE, cex = 2, main= "eCLASH", legend=T)#labels=list(font=0.5),=F)

dev.off()




