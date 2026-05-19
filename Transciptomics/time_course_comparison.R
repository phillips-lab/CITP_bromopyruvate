
library(edgeR)
library(openxlsx)
library(stringr)

# Gene annotations.
genes <- readRDS("r_data/gene_annotations_Cel.Rdata")

# Load the feature counts.
fc.bromop <- read.table("BROMOP_featureCounts.txt", header=TRUE)
rownames(fc.bromop) <- fc.bromop$Geneid
fc.bromop <- as.matrix(fc.bromop[c(1:40)])

stopifnot(all(row.names(fc.bromop)==genes$gene_id))

fc.col <- colnames(fc.bromop)
grp <- c()
age <- c()
for(j in 1:length(fc.col)) {
  col.s <- str_split(fc.col[j], "_")[[1]]
  if(col.s[3]=="BROMOP") { grp <- append(grp,"T") } else { grp <- append(grp,"C") }
  age <- append(age, col.s[2])
}

# Grouping factor.
group <- factor(grp)
age <- as.numeric(age)

# Create the data object.
y <- DGEList(counts=fc.bromop, group=group, genes=genes$gene_name)

# Filter out lowly expressed genes.
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Normalize the library sizes.
y <- calcNormFactors(y)

# Save MDS plot.
labels <- paste0(age, group)
colors <- sapply(group, function(x) ifelse(x=="C","navy","darkred"), USE.NAMES=FALSE)

pdf(file=paste0("plots/mds/mds_time_course_comparison.pdf"), width=6, height=6)
plotMDS(y, labels=labels, col=colors)
dev.off()

# Create design matrix.
design <- model.matrix(~group*age)

show(design)

# Estimate common dispersion and tagwise dispersions.
y <- estimateDisp(y, design)

# Perform quasi-likelihood F-tests.
fit <- glmQLFit(y, design, robust=TRUE)

qlf <- glmQLFTest(fit)

tags <- as.data.frame(topTags(qlf, n=nrow(qlf$table)))
tags <- cbind(row.names(tags), tags)

colnames(tags) <- c("geneWbid","gene","logFC","logCPM","F","pValue","FDR")

wb <- createWorkbook()
addWorksheet(wb=wb, sheetName="time_course")
modifyBaseFont(wb, fontSize=16)
writeDataTable(wb=wb, sheet=1, x=tags)
saveWorkbook(wb, paste0("reports/time_course_comparison.xlsx"), overwrite=TRUE)
