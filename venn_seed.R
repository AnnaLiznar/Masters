setwd('/data1/eCLASH/7-Exon-chimeras/')
#library("Biostrings")

one=read.csv('Exon_chimera_piRNA_1.tab', header= F, sep= "\t")
one=one['V2']
one$rep_1=TRUE

two=read.csv('Exon_chimera_piRNA_2.tab', header= F, sep= "\t")
two=two['V2']
two$rep_2=TRUE

three=read.csv('Exon_chimera_piRNA_3.tab', header= F, sep= "\t")
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

#png(file='./venn_123_reps.png', width=700, height = 400, res = 125)
pdf(file='./venn_123_reps.pdf', width = 5, height=3)
plot(fit_all,  alpha=0.7,counts=TRUE, cex = 2, main= "eCLASH", legend=T)#labels=list(font=0.5),=F)

#counts=TRUE
#fill=c('green', 'red', 'blue'),
dev.off()


##exxxtraccct unqiue outside and inside 

all_Xins= as.data.frame(fit_sm3_Xins)
all_Xins$seq=rownames(all_Xins)
#all= sapply(all, function(x) gsub("TRUE", 1, x))
setwd('/data1/AnLi_FM_MA/data/output_sm3_KD/grep_miRNA/names_norm/thresh_allreps/unique_names_for_GoTerms/')
protected_outside=subset(all_Xins, all_Xins$sm3_Xins==FALSE)
protected_outside=subset(protected_outside, protected_outside$protected==FALSE)
protected_outside=protected_outside[,-c(1:4)]
write.table(protected_outside, "./protected_outside_sm3_kd_Xins.txt", 
            col.names = F, row.names = F, quote = F, sep= '\t')

protected_inside=subset(all_Xins, all_Xins$sm3_Xins==TRUE)
protected_inside=subset(protected_inside, protected_inside$protected==TRUE)
write.table(protected_inside$seq, "./protected_inside_sm3_kd_Xins.txt", 
            col.names = F, row.names = F, quote = F, sep= '\t')
###################


