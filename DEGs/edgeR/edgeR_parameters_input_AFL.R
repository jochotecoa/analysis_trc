######################################################
### DEseq2 for normalization and DE identification ###
######################################################

library(DESeq2)
require("pheatmap")
require("ggplot2")
require("DESeq2")
require("edgeR")
library("lattice")

###################################################################################
###################################################################################
# PARAMETERS TO SET MANUALLY                            		 
				
basedir <- "/ngs-data-2/data/CEFIC/RNAseq_pipeline/validationDATA/Trimmed_reads_fastp/"
outputdirBASE <- paste(basedir, "edgeR/Options/", sep="")
dir.create(outputdirBASE)
inputs<-c("RSEM", "Cufflinks", "featurecounts")

#EdgeR options 
FITTYPE<-c("NegBinom", "GLM_LRT", "GLM_QLF")
	for (Y in 1:length(FITTYPE)) {
NORM<-c("TMM","TMMwsp","RLE","upperquartile")
	for (Z in 1:length(NORM)) {
Z=1
PARAMETERS= paste(FITTYPE[Y], NORM[Z],sep="_")
#Default= "NegBinom" + "TMM"

edgeRoptions<-PARAMETERS

Total_Count_cutoff<- c(1000000)
Gene_Count_cutoff<- c(10) #on average X reads per sample detected
COMP<-"AFL"
FDR<-0.01



###################################################################################
#DEFINE FUNCTIONS
############################################
Make_PCA<-function(data, NORM_TYPE){
	setwd(plotdir)
	#introduce the input file with transposing column and rows
	pc <- prcomp(t(data)) 
	PCvariance <- pc$sdev^2/sum(pc$sdev^2)
	PC1var<-signif((PCvariance[1]*100), 2)
	PC2var<-signif((PCvariance[2]*100), 2)
	#introduce the similar file with the same sample names with colors and transpose it
	Colors<-matrix(data=NA, ncol=ncol(data), nrow=1)
	colnames(Colors)<-colnames(data)
	pch<-NULL
	for (k in 1:ncol(Colors)) {
		if (substring(colnames(Colors)[k],1,3) == "Con") {Colors[1,k] <- "blue" 
			pch<- c(pch, 21)}
		if (substring(colnames(Colors)[k],1,3) == COMP) {Colors[1,k] <- "red"
			pch<- c(pch, 22)}
	}
	c<-t(Colors)
	#PCA plot of the data with colors
	plot( pc$x[ , 1:2 ], col=1, bg=c , cex=1.5, pch=pch, main=paste("PCA", COMP, NORM_TYPE, sep=" "), xlab=paste0("PC1: ", PC1var, "% variance"), ylab=paste0("PC2: ", PC2var, "% variance"))
	savePlot(filename = paste(COMP, NORM_TYPE, "PCA.png", sep="_"), type = "png")

	legend("bottomleft", legend=c(COMP, VEHI), col=c("red", "blue"), pch=c(15, 16))
	savePlot(filename = paste(COMP, NORM_TYPE, "PCA_legend.png", sep="_"), type = "png")
	#to add labels to the plot
	text(pc$x[ , 1], pc$x[ ,2 ], sapply(strsplit(rownames(c), split='_', fixed=TRUE),function(x)(x[6])), cex = 0.6)
	savePlot(filename = paste(COMP, NORM_TYPE, "PCA_sampNAME.png", sep="_"), type = "png")

} # Make_PCA function defined

###################################################################################

Run_edgeR<-function(sampleData, Samp4compare, Cont4compare, DESIGN){

## Create DESeqDesign file (for annotation):
	edgeRDesign <- matrix(data=NA, ncol=5, nrow= ncol(sampleData))
	colnames(edgeRDesign) <- c("Compound","Machine", "Run", "Comp_Mach", "Comp_Run")
	rownames(edgeRDesign) <- subKeyCompound[,"NAME"]
	edgeRDesign[,"Compound"] <- subKeyCompound[,"Abbr"]
	edgeRDesign[,"Machine"] <- subKeyCompound[,"RS_platform"]
	edgeRDesign[,"Run"] <- subKeyCompound[,"Run"]
	edgeRDesign[,"Comp_Mach"] <- paste(subKeyCompound[,"Abbr"], subKeyCompound[,"RS_platform"], sep="_")
	edgeRDesign[,"Comp_Run"] <- paste(subKeyCompound[,"Abbr"], edgeRDesign[,"Run"], sep="_")

	comparison1<- Cont4compare
	comparison2<- Samp4compare
	DESIGN<- DESIGN

	overview1<-matrix(data=NA, ncol=2, nrow=length(comparison1))
	colnames(overview1)<-c("NORM_TYPE", "DEGs")

	for (x in 1:length(comparison1)){	## for all comparisons
		condition1<- comparison1[x]	    		
		condition2<- comparison2[x]  
		print(paste(condition2, " vs ", condition1))	     		

		DE_Design <- matrix(data=NA, ncol=2)
		DE_Design <- edgeRDesign [c(grep(condition1,edgeRDesign[,DESIGN[x]]), grep(condition2,edgeRDesign[,DESIGN[x]])),]
		numeric_design <-matrix(data=NA, ncol=1, nrow=nrow(DE_Design))
		row.names(numeric_design)<-row.names(DE_Design)
		numeric_design[grep(condition1,DE_Design[,DESIGN[x]]),1] <- 1
		numeric_design[grep(condition2,DE_Design[,DESIGN[x]]),1] <- 2
		samples <- sampleData[, rownames(DE_Design) ]
		
	
		group<- DE_Design[,DESIGN[x]]
		DGElist<-DGEList(counts=samples, group=group)
		keep<- filterByExpr(DGElist) #filter out low expressed genes
		Filtered<-DGElist[keep, , keep.lib.sizes=FALSE]
		DGElist<-DGEList(counts=Filtered, group=group) #recalculate counts

		normalization<- NORM_TYPE<-paste0(inputs[i], edgeRoptions, "_FDR",FDR)
		
		#Normalizing		
		DGElist<-calcNormFactors(DGElist, method=NORM[Z]) #default method is TMM
	if (FITTYPE[Y]== "NegBinom") {
		#Dispersion (DEFAULT: Negative binomial model, estimating biological coefficient of variation)
		DGElist<-estimateCommonDisp(DGElist)
		DGElist<-estimateTagwiseDisp(DGElist)
		#DEGs identification
		DEGs<-exactTest(DGElist)		
	}
		
	if (substring(FITTYPE[Y], 1,3) == "GLM") {
		#Dispersion (Generalized linear models, fitting a log-linear model)
		DGElist<-estimateGLMCommonDisp(DGElist, numeric_design)
		DGElist<-estimateGLMTrendedDisp(DGElist, numeric_design)
		DGElist<-estimateGLMTagwiseDisp(DGElist, numeric_design)
		if (substring(FITTYPE[Y], 5,7) == "LRT") {
			#likelihood  ratiotest
			fit<-glmFit(DGElist, model.matrix(~(factor(numeric_design[,1]))))
			DEGs<-glmLRT(fit, coef=2) #2vs1, 1 =control
		}
		if (substring(FITTYPE[Y], 5,7) == "QLF") {
			#quasi-likelihood F-test
			fit<-glmQLFit(DGElist, model.matrix(~(factor(numeric_design[,1]))))
			DEGs<-glmQLFTest(fit, coef=2) #2vs1, 1 =control
		}			
	}
		
		res<-DEGs$table
		padj<- p.adjust(res[,"PValue"],"fdr") #FDR not calculated within edgeR, needs to be added manually
		res<-cbind (res, padj)
		norm_data<- DGElist$counts*DGElist$samples$norm.factors

		setwd(outputdir)
		FileName<-paste(NORM_TYPE, DESIGN[x], condition2,"vs",condition1, sep="_")
		capture.output(summary(res), file = paste(FileName, "summary_ALLgenes", "table.txt",sep="_"))
		DEsamples <-  subset(res,res$padj < FDR)
		capture.output(summary(DEsamples), file = paste(FileName, "summary_DEgenes", "table.txt",sep="_"))
		DEsamples <- DEsamples[order(DEsamples$padj),]
		filename =  paste(FileName,"DEG", "table.txt",sep="_")
		write.table(DEsamples,file=filename, sep="\t", quote=FALSE)
		DEsamples<- read.delim(file = filename, header = TRUE, row.names=NULL, sep = "\t")
		overview1[x,"NORM_TYPE"]<-paste(NORM_TYPE, condition2,"vs",condition1, sep="_")
		overview1[x,"DEGs"]<-nrow(DEsamples)

		setwd(outputdir)

		write.table(overview1,paste0("overview_", NORM_TYPE), sep="\t", quote=FALSE)
		Make_PCA(norm_data, paste(NORM_TYPE, condition2,"vs",condition1, sep="_"))

		setwd(plotdir)
		plotBCV(DGElist, main = paste0("BCV plot"))
		savePlot(filename = paste(FileName, "_BCV_plot.png", sep="_"), type = "png")
		setwd(outputdir)

	} #Close for condition
} # Run_edgeR function defined

#####################################################################################################################
#
# ALL FUNCTIONS DEFINED
# 
#####################################################################################################################







### Input = RSEM
i=1
outputdirINTERM <- paste(outputdirBASE, inputs[i], "/", sep="")
dir.create(outputdirINTERM)
outputdir <- paste(outputdirINTERM, edgeRoptions, "/", sep="")
dir.create(outputdir)
plotdir<- paste(outputdir, "/plots/", sep="")
dir.create(plotdir)
sampledir <- paste(basedir, inputs[i], "/", sep="")

###################################################################################
###################################################################################

setwd(sampledir)

# Load samplekey 
FullKey<- read.delim("/ngs-data-2/data/CEFIC/samplekeyAFLplatform.csv", stringsAsFactors=FALSE, sep=",", header=TRUE,  quote="\"", row.names=NULL)
ID_SamplesAll<-FullKey[,"RS_run"]

#Merge sample files 
A<- read.delim(paste0(ID_SamplesAll[1], ".genes.results"), stringsAsFactors=FALSE, header=TRUE,  quote="\"", row.names=NULL)
A<- A[,c("gene_id", "expected_count")]
B<- read.delim(paste0(ID_SamplesAll[2], ".genes.results"), stringsAsFactors=FALSE, header=TRUE,  quote="\"", row.names=NULL)
B<- B[,c("gene_id", "expected_count")]
merge<- merge(A, B, by="gene_id", all=TRUE, suffixes=c(ID_SamplesAll[1], ID_SamplesAll[2])) 

for (S in 3:length(ID_SamplesAll)) {
	C<- read.delim(paste0(ID_SamplesAll[S], ".genes.results"), stringsAsFactors=FALSE, header=TRUE,  quote="\"", row.names=NULL)
print(ID_SamplesAll[S])
	C<- C[,c("gene_id", "expected_count")]
	merge<- merge(merge, C, by="gene_id", all=TRUE, suffixes=c(ID_SamplesAll[S-1], ID_SamplesAll[S])) 
}
row.names(merge)<-merge[,1]
merge<-merge[,2:ncol(merge)]
colnames(merge)<-FullKey[,"NAME"]
DATA_SamplesAll<-merge


# Subset samplekey RNAseq Analysis
VEHI<- FullKey[grep(COMP,FullKey[,"Abbr"]), "vehicle_route_designation"][1]
subKeyCompound<- FullKey[c(grep(COMP,FullKey[,"Abbr"]), grep(VEHI,FullKey[,"Abbr"])),] 
ID_SamplesCompound<-subKeyCompound[,"NAME"]

# Subset sampledata RNAseq Analysis
sampleData<- DATA_SamplesAll[,c(ID_SamplesCompound)]
# replace NA 
sampleData[ is.na(sampleData) ] <- 0 
# remove samples with low readcount (below set Total_Count_cutoff)
sampleData<- sampleData[,(colSums(sampleData)>Total_Count_cutoff)]
# remove genes with low readcount (below set cutoff)
#sampleData<- sampleData[(rowSums(sampleData)>Gene_Count_cutoff*ncol(sampleData)),]



##########################
# amount expressed genes #
##########################
overview_amount_expressed_genes<-matrix(data=NA, ncol=3, nrow=1)
colnames(overview_amount_expressed_genes)<-c("Compound", "mean", "sd")

#if (FIRST COMPOUND OF THE LIST) {
subset_VEHI<-sampleData[,grep(VEHI,colnames(sampleData))]
expressed<-NULL
for (x in 1:ncol(subset_VEHI)) {expressed<-c(expressed, sum(subset_VEHI[,x]>0))}
overview_amount_expressed_genes<-rbind(overview_amount_expressed_genes, matrix(data=c(VEHI, mean(expressed), sd(expressed)), ncol=3, nrow=1))
overview_amount_expressed_genes<-overview_amount_expressed_genes[2,]
#}	
subset_COMP<-sampleData[,grep(COMP,colnames(sampleData))]
expressed<-NULL
for (x in 1:ncol(subset_COMP)) {expressed<-c(expressed, sum(subset_COMP[,x]>0))}
overview_amount_expressed_genes<-rbind(overview_amount_expressed_genes, matrix(data=c(COMP, mean(expressed), sd(expressed)), ncol=3, nrow=1))

setwd(outputdir)	
write.table(overview_amount_expressed_genes,paste0(inputs[i], edgeRoptions, "_overview_amount_expressed_genes.csv"), sep=",", quote=FALSE, col.names=T, row.names=F)	
###########


#################
# Plot raw data #
#################

Make_PCA(sampleData, paste0(inputs[i], edgeRoptions, "_raw")) #(data, NORM_TYPE)


##########
# edgeR #
##########

Cont4compare<- c(VEHI, "HiScanSQ", "AFL_B", "NN_OG_B", "NN_OG_A", "NN_OG_B","NN_OG_B", "NN_OG_A")
Samp4compare<- c(COMP, "HiSeq2000", "AFL_A", "NN_OG_A", "AFL_A", "AFL_B", "AFL_A", "AFL_B")
DESIGN<- c("Compound", "Machine", "Comp_Run", "Comp_Run","Comp_Run","Comp_Run","Comp_Run","Comp_Run")
Run_edgeR(sampleData, Samp4compare, Cont4compare, DESIGN)

#Done with RSEM
################################################################################################################################












### Input = Cufflinks
i=2
#cufflinks failed: empty read_count table

#Done with Cufflinks
################################################################################################################################














### Input = featurecounts
i=3
outputdirINTERM <- paste(outputdirBASE, inputs[i], "/", sep="")
dir.create(outputdirINTERM)
outputdir <- paste(outputdirINTERM, edgeRoptions, "/", sep="")
dir.create(outputdir)
plotdir<- paste(outputdir, "/plots/", sep="")
dir.create(plotdir)
sampledir <- paste(basedir, inputs[i], "/", sep="")

###################################################################################
###################################################################################

setwd(sampledir)

# Load samplekey (#AFL_Rep1B deleted, quality dip R2 + outlyer PCA)
allCompKey<-read.delim("/ngs-data-2/data/CEFIC/samplekeyALLplatform.csv", stringsAsFactors=FALSE, sep=",", header=TRUE,  quote="\"", row.names=NULL)
FullKey<- read.delim("/ngs-data-2/data/CEFIC/samplekeyAFLplatform.csv", stringsAsFactors=FALSE, sep=",", header=TRUE,  quote="\"", row.names=NULL)
ID_SamplesCompound<-FullKey[,"NAME"]

#Load sample file 
DATA_SamplesAll<- read.delim("FeatureCounts_basic.txt", stringsAsFactors=FALSE, sep="\t", header=FALSE,  quote="\"", row.names=NULL)
rownames(DATA_SamplesAll)<-DATA_SamplesAll[,1]
DATA_SamplesAll<- DATA_SamplesAll[c(2:nrow(DATA_SamplesAll)),c(7:ncol(DATA_SamplesAll))]

#set colnames
colnames<-NULL
nameTAB<-allCompKey[,c("NAME", "RS_run")]
rownames(nameTAB)<-nameTAB[,"RS_run"]
for (col in 1:ncol(DATA_SamplesAll)) {
	colnames<-c(colnames, nameTAB[substring(DATA_SamplesAll[1,col],0,10),1])	
}
colnames(DATA_SamplesAll)<-colnames
DATA_SamplesAll<-DATA_SamplesAll[2:nrow(DATA_SamplesAll),]

# Subset samplekey RNAseq Analysis
VEHI<- FullKey[grep(COMP,FullKey[,"Abbr"]), "vehicle_route_designation"][1]
subKeyCompound<- FullKey[c(grep(COMP,FullKey[,"Abbr"]), grep(VEHI,FullKey[,"Abbr"])),] 
ID_SamplesCompound<-subKeyCompound[,"NAME"]

# Subset sampledata RNAseq Analysis
sampleData<- DATA_SamplesAll[,c(ID_SamplesCompound)]
write.table(sampleData, "temp.csv", sep=",", quote=FALSE, col.names=T, row.names=T)
sampleData<- read.delim("temp.csv", stringsAsFactors=FALSE, sep=",", header=TRUE,  quote="\"", row.names=1)
colnames(sampleData)<-ID_SamplesCompound

# replace NA 
sampleData[ is.na(sampleData) ] <- 0 
# remove samples with low readcount (below set Total_Count_cutoff)
sampleData<- sampleData[,(colSums(sampleData)>Total_Count_cutoff)]
# remove genes with low readcount (below set cutoff)
#sampleData<- sampleData[(rowSums(sampleData)>Gene_Count_cutoff*ncol(sampleData)),]




##########################
# amount expressed genes #
##########################
overview_amount_expressed_genes<-matrix(data=NA, ncol=3, nrow=1)
colnames(overview_amount_expressed_genes)<-c("Compound", "mean", "sd")

#if (FIRST COMPOUND OF THE LIST) {
subset_VEHI<-sampleData[,grep(VEHI,colnames(sampleData))]
expressed<-NULL
for (x in 1:ncol(subset_VEHI)) {expressed<-c(expressed, sum(subset_VEHI[,x]>0))}
overview_amount_expressed_genes<-rbind(overview_amount_expressed_genes, matrix(data=c(VEHI, mean(expressed), sd(expressed)), ncol=3, nrow=1))
overview_amount_expressed_genes<-overview_amount_expressed_genes[2,]
#}	
subset_COMP<-sampleData[,grep(COMP,colnames(sampleData))]
expressed<-NULL
for (x in 1:ncol(subset_COMP)) {expressed<-c(expressed, sum(subset_COMP[,x]>0))}
overview_amount_expressed_genes<-rbind(overview_amount_expressed_genes, matrix(data=c(COMP, mean(expressed), sd(expressed)), ncol=3, nrow=1))

setwd(outputdir)	
write.table(overview_amount_expressed_genes,paste0(inputs[i], edgeRoptions, "_overview_amount_expressed_genes.csv"), sep=",", quote=FALSE, col.names=T, row.names=F)	
###########


#################
# Plot raw data #
#################

Make_PCA(sampleData, paste0(inputs[i], edgeRoptions, "_raw")) #(data, NORM_TYPE)


##########
# egdeR #
##########

Cont4compare<- c(VEHI, "HiScanSQ", "AFL_B", "NN_OG_B", "NN_OG_A", "NN_OG_B","NN_OG_B", "NN_OG_A")
Samp4compare<- c(COMP, "HiSeq2000", "AFL_A", "NN_OG_A", "AFL_A", "AFL_B", "AFL_A", "AFL_B")
DESIGN<- c("Compound", "Machine", "Comp_Run", "Comp_Run","Comp_Run","Comp_Run","Comp_Run","Comp_Run")
Run_edgeR(sampleData, Samp4compare, Cont4compare, DESIGN)


#Done with featurecounts
################################################################################################################################







}} #for edgeRoptions
########################################################################################################################



	
print("Script completed !!!")        

