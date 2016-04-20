### R code from vignette source 'FunciSNP_vignette.Rnw'

###################################################
### code chunk number 1: packages
###################################################
options(width=80);
library(FunciSNP);
package.version("FunciSNP");


###################################################
### code chunk number 2: FunciSNP_vignette.Rnw:204-210
###################################################
## Full path to the example GWAS SNP regions file for Glioblastoma 
#  (collected from SNPedia on Jan 2012)
glioma.snp <- file.path(system.file('extdata', package='FunciSNP'), 
dir(system.file('extdata',package='FunciSNP'), pattern='.snp$'));
gsnp <- read.delim(file=glioma.snp,sep=" ",header=FALSE);
gsnp;


###################################################
### code chunk number 3: FunciSNP_vignette.Rnw:228-256
###################################################
## Full path to the example biological features BED files 
#  derived from the ENCODE project for Glioblastoma U-87 cell lines.
glioma.bio <- system.file('extdata',package='FunciSNP');
#user supplied biofeatures
as.matrix(list.files(glioma.bio, pattern='.bed$'));
#FunciSNP builtin biofeatures
as.matrix(list.files(paste(glioma.bio, "/builtInFeatures", sep=""),
            pattern='.bed$'));
nrsf.filename <- list.files(glioma.bio, pattern='.bed$')[1];
pol2.filename <- list.files(glioma.bio, pattern='.bed$')[2];
Ctcf <- ctcf_only
Dnase1 <- encode_dnase1_only
Dnase1Ctcf <- encode_dnase1_with_ctcf
Faire <- encode_faire
Promoters <- known_gene_promoters
Nrsf <- read.delim(file=paste(glioma.bio, nrsf.filename,sep="/"), sep="\t",
header=FALSE);
PolII <- read.delim(file=paste(glioma.bio, pol2.filename,sep="/"), sep="\t",
header=FALSE);
dim(Nrsf);
dim(PolII);
dim(Ctcf);
dim(Dnase1);
dim(Dnase1Ctcf);
dim(Faire);
dim(Promoters);
## Example of what the BED format looks like:
head(Nrsf);


###################################################
### code chunk number 4: FunciSNP_vignette.Rnw:270-279
###################################################
## FunciSNP analysis, extracts correlated SNPs from the 
## 1000 genomes db ("ncbi" or "ebi") and finds overlaps between 
## correlated SNP and biological features and then 
## calculates LD (Rsquare, Dprime, distance, p-value).
## Depending on number of CPUs and internet connection, this step may take 
## some time. Please consider using a unix machine to access multiple cores.

# glioma <- getFSNPs(snp.regions.file=glioma.snp, bio.features.loc = glioma.bio,
# bio.features.TSS=FALSE);


###################################################
### code chunk number 5: FunciSNP_vignette.Rnw:285-287
###################################################
data(glioma);
class(glioma);


###################################################
### code chunk number 6: FunciSNP_vignette.Rnw:296-297
###################################################
glioma;


###################################################
### code chunk number 7: FunciSNP_vignette.Rnw:311-312
###################################################
summary(glioma);


###################################################
### code chunk number 8: FunciSNP_vignette.Rnw:339-349
###################################################
glioma.anno <- FunciSNPAnnotateSummary(glioma);
class(glioma.anno);
gl.anno <- glioma.anno;
## remove rownames for this example section.
rownames(gl.anno) <- c(1:length(rownames(gl.anno)))
dim(gl.anno);
names(gl.anno);
head(gl.anno[, c(1:18,20:28)]);
summary(gl.anno[, c(1:18,20:28)]);
rm(gl.anno);


###################################################
### code chunk number 9: FunciSNP_vignette.Rnw:384-385
###################################################
FunciSNPtable(glioma.anno, rsq=0.44);


###################################################
### code chunk number 10: FunciSNP_vignette.Rnw:393-394
###################################################
FunciSNPtable(glioma.anno, rsq=0.44, geneSum=TRUE);


###################################################
### code chunk number 11: FunciSNP_vignette.Rnw:405-406
###################################################
FunciSNPsummaryOverlaps(glioma.anno)


###################################################
### code chunk number 12: FunciSNP_vignette.Rnw:412-413
###################################################
FunciSNPsummaryOverlaps(glioma.anno, rsq=0.44)


###################################################
### code chunk number 13: FunciSNP_vignette.Rnw:424-430
###################################################
rs6010620 <- FunciSNPidsFromSummary(glioma.anno, tagsnpid="rs6010620",
num.features=2, rsq=0.44);
#summary(rs6010620);
dim(rs6010620);
class(rs6010620);
## See FunciSNPbed to visualize this data in a genome browser.


###################################################
### code chunk number 14: FunciSNP_vignette.Rnw:452-455
###################################################
pdf("glioma_dist.pdf")
FunciSNPplot(glioma.anno)
dev.off()


###################################################
### code chunk number 15: FunciSNP_vignette.Rnw:487-489
###################################################
FunciSNPplot(glioma.anno, splitbysnp=TRUE)
ggsave("glioma_dist_bysnp.pdf")


###################################################
### code chunk number 16: FunciSNP_vignette.Rnw:506-509
###################################################
pdf("glioma_heatmap.pdf")
FunciSNPplot(glioma.anno, heatmap=TRUE, rsq = 0.1)
dev.off()


###################################################
### code chunk number 17: FunciSNP_vignette.Rnw:530-532
###################################################
## Following will output a series of plots for each biofeature at rsq=0.5
FunciSNPplot(glioma.anno, tagSummary=TRUE, rsq=0.5)


###################################################
### code chunk number 18: FunciSNP_vignette.Rnw:582-585
###################################################
pdf("glioma_genomic_sum_rcut.pdf")                                               
FunciSNPplot(glioma.anno, rsq=0.5, genomicSum=TRUE, save=FALSE)                  
dev.off()                                                                        


###################################################
### code chunk number 19: FunciSNP_vignette.Rnw:627-630
###################################################
## will output to current working directory.
FunciSNPbed(glioma.anno, rsq=0.22);
# FunciSNPbed(rs6010620, rsq=0.5);


###################################################
### code chunk number 20: FunciSNP_vignette.Rnw:665-666
###################################################
toLatex(sessionInfo())                                                           


