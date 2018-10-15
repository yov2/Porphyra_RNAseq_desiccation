source("https://bioconductor.org/biocLite.R")
biocLite("sva")
library(sva)
biocLite("pamr")
library(pamr)
library(edgeR)
library(limma)
library(DESeq2)
biocLite("IHW")
library(IHW)

##input count matrixs and sample information
setwd("/Users/helen/Desktop/bioinformatic/script_rna_pipeline/5_7")
A = read.csv("htseq_star_1.csv",header = TRUE )
row_name_B=A[,1]
B=matrix(0, nrow=dim(A)[1],ncol=24)
B=A[,-1]
row.names(B)=row_name_B
head(B)

C= read.csv("pheno_1.csv",header = TRUE)
row_name_D=C[,1]
coldata=matrix(0, nrow=dim(C)[1],ncol=24)
coldata=C[,-1]
row.names(coldata)=row_name_D
head(coldata)

##filter the reads 
filter = apply(B, 1, function(x) length(x[x>3]) >=10)
filtered = B[filter,]
head(filtered)
dat0=as.matrix(filtered)

# to check if coldata and dat0 match
all(rownames(coldata) %in% colnames(dat0))
all(rownames(coldata) == colnames(dat0))

#read into DESeq2
dds = DESeqDataSetFromMatrix(countData = dat0,
                             colData = coldata,
                             design = ~ condition)
dds
rld = rlog(dds,blind = FALSE)

##head(assay(rld),3)
rld_dataframe = as.matrix(assay(rld))

#compute null and full model
mod = model.matrix(~as.factor(condition), data= coldata)
mod0=model.matrix(~1,data=coldata)
##mod0 = cbind(mod1[,1])

# use default method gave 3 surrogate variables, while using leek method gave 20 surrogate variables.
n.sv = num.sv(rld_dataframe,mod)
n.sv
## n.sv = 2

## the highest n.sv worked is 2
svseq = sva(rld_dataframe, mod, mod0, n.sv=2)
# adjust for surrogate variables
ddssva = dds
ddssva$SV1 = svseq$sv[,1]
ddssva$SV2 = svseq$sv[,2]
design(ddssva) = ~ SV1 + SV2 +condition
ddssva_DE = DESeq(ddssva)

##perform the same test after adjusting for surrogate variables
modSv=cbind(mod,svseq$sv)
mod0Sv = cbind (mod0,svseq$sv)

##input the filtered table matrix into limma package and perform liner fit
dge = DGEList(counts=filtered)
group = as.factor(c(rep("Fresh", 6), rep("Dehydration",6), rep("Desiccation",6), rep("Rehydration",6)))
dge$samples$group = group

#normalize library size
dge_Nor = calcNormFactors(dge, method = "TMM")

##counts are converted to logCPM values 
logCPM = cpm (dge_Nor, log=TRUE, prior.count = 0.25)
plotMDS(logCPM, col = as.numeric(group),cex=0.5)

##install.packages("RColorBrewer")
##library("RColorBrewer")
##pca
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(logCPM, labels=group, col=col.group)
title(main="A. Sample groups")

##perform PCA using original data (PCA plot from read count matrix from RNA-Seq/https://www.biostars.org/p/282685/#282691)
project.pca = prcomp(t(logCPM))
summary(project.pca)
#Determine the proportion of variance of each component
#Proportion of variance equals (PC stdev^2) / (sum all PCs stdev^2)
project.pca.proportionvariances = ((project.pca$sdev^2)/(sum(project.pca$sdev^2)))*100
#bar plot shows the proportion of each variance of component
barplot(project.pca.proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(project.pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))

##bi-plots using original data
##par(mar=c(4,4,4,4), mfrow=c(1,3), cex=1.0, cex.main=0.8, cex.axis=0.8)
par(mar=c(4,4,4,4),cex=1.0, cex.main=0.8, cex.axis=0.8)
#Plots scatter plot for PC 1 and 2
plot(project.pca$x, type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"))
points(project.pca$x, col=as.numeric(group), pch=16, cex=1)
#Plots scatter plot for PC 3 and 4
plot(project.pca$x[,3], project.pca$x[,4], type="n", main="Principal components analysis bi-plot", xlab=paste("PC3, ", round(project.pca.proportionvariances[3], 2), "%"), ylab=paste("PC4, ", round(project.pca.proportionvariances[4], 2), "%"))
points(project.pca$x[,3], project.pca$x[,4], col=as.numeric(group), pch=16, cex=1)

#extract the covariates from modSv
cov=modSv[,5:6]
logCPM_nobatch = removeBatchEffect(logCPM, covariates = cov)
##pca
#par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(logCPM_nobatch, labels=group, col=col.group)
title(main="A. Sample groups")

##pca for after removing covariates
project.pca.nobatch = prcomp(t(logCPM_nobatch))
project.pca.proportionvariances.nobatch = ((project.pca.nobatch$sdev^2)/(sum(project.pca.nobatch$sdev^2)))*100
par(mar=c(4,4,4,4),cex=1.0, cex.main=0.8, cex.axis=0.8)
#Plots scatter plot for PC 1 and 2
plot(project.pca.nobatch$x, type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(project.pca.proportionvariances.nobatch[1], 2), "%"), ylab=paste("PC2, ", round(project.pca.proportionvariances.nobatch[2], 2), "%"))
points(project.pca.nobatch$x, col=as.numeric(group), pch=16, cex=0.1)
label=c('F1','F2','F3','F4','F5','F6','H1','H2','H3','H4','H5','H6','D1','D2','D3','D4','D5','D6','R1','R2','R3','R4','R5','R6')
text(project.pca.nobatch$x,labels=label,col=as.numeric(group))




#Plots scatter plot for PC 3 and 4
plot(project.pca.nobatch$x[,3], project.pca.nobatch$x[,4], type="n", main="Principal components analysis bi-plot", xlab=paste("PC3, ", round(project.pca.proportionvariances.nobatch[3], 2), "%"), ylab=paste("PC4, ", round(project.pca.proportionvariances.nobatch[4], 2), "%"))
points(project.pca.nobatch$x[,3], project.pca.nobatch$x[,4], col=as.numeric(group), pch=16, cex=1)


##desiccation vs fresh
ressva_Desiccation_fresh = results(ddssva_DE, contrast = c("condition", "Desiccation", "Fresh"),filterFun=ihw, alpha=0.05)
summary(ressva_Desiccation_fresh,na.rm=TRUE)
sum(ressva_Desiccation_fresh$padj < 0.05, na.rm=TRUE)
ressva_Desiccation_fresh_Sig = subset (ressva_Desiccation_fresh, padj <0.05)
ressva_Desiccation_fresh_Sig_order = ressva_Desiccation_fresh_Sig[order(ressva_Desiccation_fresh_Sig$log2FoldChange),]
ressva_Desiccation_fresh_Sig_order_up = ressva_Desiccation_fresh_Sig[order(ressva_Desiccation_fresh_Sig$log2FoldChange, decreasing = TRUE),]
write.csv(as.data.frame(ressva_Desiccation_fresh_Sig_order), file = "./desiccation_vs_fresh.csv")

##dehydrtion vs fresh
ressva_Dehydration_fresh = results(ddssva_DE, contrast = c("condition", "Dehydration", "Fresh"),filterFun=ihw, alpha=0.05)
summary(ressva_Dehydration_fresh,na.rm=TRUE)
sum(ressva_Dehydration_fresh$padj < 0.05, na.rm=TRUE)
ressva_Dehydration_fresh_Sig = subset (ressva_Dehydration_fresh, padj <0.05)
ressva_Dehydration_fresh_order = ressva_Dehydration_fresh_Sig[order(ressva_Dehydration_fresh_Sig$log2FoldChange),]
write.csv(as.data.frame(ressva_Dehydration_fresh_order), file = "./dehydration_vs_fresh.csv")

##rehydration vs fresh
ressva_Rehydration_fresh = results(ddssva_DE, contrast = c("condition", "Rehydration", "Fresh"),filterFun=ihw, alpha=0.05)
summary(ressva_Rehydration_fresh,na.rm=TRUE)
sum(ressva_Rehydration_fresh$padj < 0.05, na.rm=TRUE)
ressva_Rehydration_fresh_Sig = subset (ressva_Rehydration_fresh, padj <0.05)
ressva_Rehydration_fresh_order = ressva_Rehydration_fresh_Sig[order(ressva_Rehydration_fresh_Sig$log2FoldChange),]
write.csv(as.data.frame(ressva_Rehydration_fresh_order), file = "./Rehydration_vs_fresh.csv")

##desiccation vs dehydration
ressva_Desiccation_Dehydration = results(ddssva_DE, contrast = c("condition", "Desiccation", "Dehydration"),filterFun=ihw, alpha=0.05)
summary(ressva_Desiccation_Dehydration,na.rm=TRUE)
sum(ressva_Desiccation_Dehydration$padj < 0.05, na.rm=TRUE)
ressva_Desiccation_Dehydration_Sig = subset (ressva_Desiccation_Dehydration, padj <0.05)
ressva_Desiccation_Dehydration_order = ressva_Desiccation_Dehydration_Sig[order(ressva_Desiccation_Dehydration_Sig$log2FoldChange),]
write.csv(as.data.frame(ressva_Desiccation_Dehydration_order), file = "./desiccation_vs_dehydration.csv")

##rehydration vs dehydration
ressva_Rehydration_dehydration = results(ddssva_DE, contrast = c("condition", "Rehydration", "Dehydration"),filterFun=ihw, alpha=0.05)
summary(ressva_Rehydration_dehydration,na.rm=TRUE)
sum(ressva_Rehydration_dehydration$padj < 0.05, na.rm=TRUE)
ressva_Rehydration_dehydration_Sig = subset (ressva_Rehydration_dehydration, padj <0.05)
ressva_Rehydration_dehydration_order = ressva_Rehydration_dehydration_Sig[order(ressva_Rehydration_dehydration_Sig$log2FoldChange),]
write.csv(as.data.frame(ressva_Rehydration_dehydration_order), file = "./Rehydration_vs_dehydration.csv")

##rehydration vs desiccation
ressva_Rehydration_desiccation = results(ddssva_DE, contrast = c("condition", "Rehydration", "Desiccation"),filterFun=ihw, alpha=0.05)
summary(ressva_Rehydration_desiccation,na.rm=TRUE)
sum(ressva_Rehydration_desiccation$padj < 0.05, na.rm=TRUE)
ressva_Rehydration_desiccation_Sig = subset (ressva_Rehydration_desiccation, padj <0.05)
ressva_Rehydration_desiccation_order = ressva_Rehydration_desiccation_Sig[order(ressva_Rehydration_desiccation_Sig$log2FoldChange),]
write.csv(as.data.frame(ressva_Rehydration_desiccation_order), file = "./Rehydration_vs_desiccation.csv")
