library(ggplot2)
library(data.table)
library(parallel)
library(Biostrings)
library(Rmisc)

setwd("/common/WORK/mfabijanic/Sponges/CURRENT")
doforgenome <- function(gname, gpath){
	Esu4 <- readDNAStringSet(gpath)
	imena <- c("A","C","G","T")
	sve <- mclapply(1:length(Esu4), function(indeks){
		x <- strsplit(as.character(Esu4[[indeks]]),"")[[1]]
		x <- x[x%in%imena]
		bwt <- data.table(ispred=x, iza=dplyr::lead(x))
		bwt[,get("imena"):=lapply(imena,function(x)mean(na.omit(iza)==x)),ispred]
		bwt <- unique(bwt[,c("ispred",imena), with=FALSE][order(ispred)])
		bwt[,seqid:=indeks]
		bwt
	}, mc.cores=22)
	svesve <- do.call("rbind",sve)
	sve <- melt(svesve,id.vars=c("ispred","seqid"), variable.name="iza", value.name="percentage")
	ggplot(sve, aes(iza, percentage))+stat_boxplot(geom="errorbar")+facet_grid(ispred~.)+theme_light()+geom_boxplot(outlier.color=NA)+xlab("Nucleotide after")+geom_hline(yintercept=0.25,col="red")+coord_flip()+ggtitle(gname)
}
eunapius <- doforgenome("Eunapius", "EsuV4polished.fasta")
suberites <- doforgenome("Suberites", "Sub_repilonwtdbg.fasta")
ephydatia <- doforgenome("Ephydatia", "Emu_Illumina_out.padded.fasta")

pdf("assembledGenomes.pdf")
multiplot(eunapius, suberites, ephydatia, ncol=3)
dev.off()
