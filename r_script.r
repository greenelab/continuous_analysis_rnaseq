# run r work
base_dir <- "/drone/src/github.com/greenelab/continuous_analysis_rnaseq"
sample_id <- dir(file.path(base_dir,"kallisto_output"))
sample_id
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "kallisto_output", id))
kal_dirs

s2c <- read.table(file.path(base_dir, "samples.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample, condition)
s2c

s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)

# this analysis is from: https://benchtobioinformatics.wordpress.com/2015/07/10/using-kallisto-for-gene-expression-analysis-of-published-rnaseq-data/
# perform PCA
source("https://bioconductor.org/biocLite.R")
library("edgeR")
library('biomaRt')
library("genefilter")

library("ggplot2")
library("stringr")
library("plyr")
library("RColorBrewer")
for (i in 1:length(kal_dirs)){
	print(kal_dirs[i])
	tmp = read.table(file = paste0(kal_dirs[i],"/abundance.tsv"), header = T)
	assign(kal_dirs[i], tmp)
}
sample_list = mget(kal_dirs)
ListNames <- Map(function(x, i) setNames(x, ifelse(names(x) %in% "target_id", names(x), sprintf("%s.%d", names(x), i))), sample_list, seq_along(sample_list))
# merge full kalisto table
# have to use Reduce() function as merge() will only merge two data.frames
full_results <- Reduce(function(...) merge(..., by='target_id', all = T), ListNames)
 
# subset estimate count values (correcting from TPM values)
count_vals <- full_results[, grep("est_counts", names(full_results))]
 
# now we have to assign column and row IDs to tpm table
row.names(count_vals) <- full_results$target_id
 
# setting column names is a little trickier as each file has "output_SRR..." as IDs
# we will use regrex with stringr and dplyr packages to get down to the SRR file IDs
splitNames <- str_split(string=kal_dirs, pattern = "_", n = 2)
df <- ldply(splitNames, data.frame)
FinalNames <- df[df$X..i.. != "output",]
 
FinalNames
# [1] SRR1654626 SRR1654628 SRR1654633 SRR1654636 SRR1654637 SRR1654639 SRR1654641 SRR1654643
# 9 Levels: output SRR1654626 SRR1654628 SRR1654633 SRR1654636 SRR1654637 ... SRR1654643
# actual IDs
ActualNames <- c("mN5",'mN7','mP1','mP4','mP5','mT3','mT5') #,'mT8')
Groups <- data.frame(Subtypes = c("mN","mN",'mP','mP','mP','mT','mT')) #,'mT'))
 
Names_count <- data.frame(FinalNames, ActualNames, Groups)
Names_count
 
colnames(count_vals) <- ActualNames
head(count_vals)

des <- factor(Groups$Subtypes)
design <- model.matrix(~0+des)
colnames(design) <- levels(des)
design

# filter out lowly expressed genes
A <- rowSums(count_vals)
isexpr <- A > 1
table(isexpr)
# isexpr
# FALSE TRUE 
# 21367 66831 
count_vals <- count_vals[isexpr, ]
dim(count_vals)
y <- DGEList(counts = count_vals)
 
# lets do some exploratory data analysis of this RNAseq dataset
# first with Principle componenet analysis like figure 5A in Boj et al., 
# and compare whether TPM normalization can replicate their data 
 
log_count <- log2(y$counts + 1)
 
pcCount <- prcomp(t(log_count))
 
# 2D PCA plot with ggplot2
scoresCount <- data.frame(pcCount$x, Names_count) 
 
brewer.pal(3,"Set2")
my_palette <- brewer.pal(3,"Set2") # "#66C2A5" "#FC8D62" "#8DA0CB"

pcaTpm = (pcaCount <- ggplot(scoresCount, aes(x=PC1, y=PC2)) +
 geom_point(size = 4, aes(col=factor(scoresCount$Subtypes))) +
 ggtitle("Principal Components\nUsing Kallisto Estimated Counts") +
 geom_text(aes(label=scoresCount$ActualNames),hjust=0.5, vjust=1) + 
 theme_minimal())
ggsave(filename=file.path(base_dir, "results/PCA.png"), plot=pcaTpm, width=10, height=8)

# volcano plots
exp <- voom(counts = y$counts, design = design, plot = T)
fit <- lmFit(exp, design = design)
 
# set up contrast matrix of conditions
cont.matrix <- makeContrasts("mT-mN",'mT-mP','mN-mP', levels = design)
fit <- contrasts.fit(fit, cont.matrix)
fit <- eBayes(fit)
topTable(fit, adjust="BH", number = Inf, p=0.05)
summary(decideTests(fit))
 
# compare the mT vs mN and mP vs mN contrasts and sort by LogFC
options(digits = 3)
 
# comparing normal organoids to Tumor organoids 
res_mTvsmN <- (topTable(fit, coef = 1, adjust = 'BH', n = Inf, p = 1)[,-2])
res_mTvsmN <- res_mTvsmN [order(-res_mTvsmN$logFC),]
head(res_mTvsmN)
nrow(res_mTvsmN)

res_mTvsmN$threshold = as.factor(abs(res_mTvsmN$logFC) > 2 & res_mTvsmN$P.Value < 0.05/(dim(res_mTvsmN)[1]))
# Make a basic volcano plot
g = ggplot(data=res_mTvsmN, aes(x=logFC, y=-log10(P.Value), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(c(-10, 10)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value")
ggsave(filename=file.path(base_dir, "results/res_mTvsmN.png"), plot=g, width=10, height=10)

# comparing normal organoids to PanIN organoids
res_mNvsmP <- (topTable(fit, coef = 3, adjust='BH', n = Inf, p = 1)[,-2])
res_mNvsmP <- res_mNvsmP[order(-res_mNvsmP$logFC),]
head(res_mNvsmP)
nrow(res_mNvsmP)

res_mNvsmP$threshold = as.factor(abs(res_mNvsmP$logFC) > 2 & res_mNvsmP$P.Value < 0.05/(dim(res_mNvsmP)[1]))
# Make a basic volcano plot
g = ggplot(data=res_mNvsmP, aes(x=logFC, y=-log10(P.Value), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(c(-10, 10)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value")
ggsave(filename=file.path(base_dir, "results/res_mNvsmP.png"), plot=g, width=10, height=10)


# Sleuth is a new tool to peform differential expresion with diffrential expresion, so we'll do that as well
s2c <- read.table(file.path(base_dir, "samples.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample, condition)
s2c

s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)

# library('sleuth')
# so <- sleuth_prep(s2c, ~ condition)
# so <- sleuth_fit(so)
# so <- sleuth_fit(so, ~1, 'reduced')
# so <- sleuth_lrt(so, 'reduced', 'full')
# models(so)
# results_table <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
# write.table(results_table, file='results/raw_sleuth_results.txt', sep="\t")

# # include gene names
# mart <- biomaRt::useMart(biomart = 'ensembl', dataset = "mmusculus_gene_ensembl")
# t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
#     "external_gene_name"), mart = mart)
# t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
#   ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
# so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)
# so <- sleuth_fit(so)
# so <- sleuth_fit(so, ~1, 'reduced')
# so <- sleuth_lrt(so, 'reduced', 'full')
# results_table <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
# write.table(results_table, file='results/named_sleuth_results.txt', sep="\t")