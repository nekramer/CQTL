contrast_shrink_and_format <- function(dds, contrast, gse){
de_genes_shrink <- lfcShrink(dds, contrast = contrast, type = "ashr", format = "GRanges") %>%
plyranges::names_to_column("gene_id")
de_genes_shrink <-
inner_join(x = as.data.frame(de_genes_shrink),
y = as.data.frame(rowData(gse)) %>%
dplyr::select(c("gene_id", "symbol", "tx_ids")),
by = "gene_id") %>%
makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
keepStandardChromosomes(pruning.mode = "coarse") %>%
as.data.frame()
return(de_genes_shrink)
}
donorSamplesheet <- read_csv("data/donorSamplesheet.csv") %>%
mutate(Race = replace_na(Race, "Unknown"))

fit_de_model <- function(dds, filter = 0.05){
## Filter out lowly expressed genes
# 10 reads in at least filter% of samples
keep <- rowSums(counts(dds) >= 10) >= ceiling(nrow(colData(gse))*filter)
dds <- dds[keep,]
## Fit model
dds_fit <- DESeq(dds)
return(dds_fit)
}

###### Sex differences in FNF
gse_fnf <- gse[, gse$Condition == "FNF"]

## Join gse with additional donor samplesheet info
colData(gse) <- as(left_join(as.data.frame(colData(gse)),
donorSamplesheet[,c("Donor", "Sex", "Age", "Race")],
by = "Donor"), "DataFrame")
## Put Age into categories
gse$Age_group <- cut(gse$Age, breaks = seq(29, 91, 10),
labels = c("30-39", "40-49", "50-59", "60-69", "70-79", "80-89"))
## Convert to factors
colData(gse)[] <- lapply(colData(gse), factor)
## Add colnames based on sample name/Sex
colnames(gse) <- paste0(colData(gse)[,"names"], "_", colData(gse)[,"Sex"])

fnf_desex <- shrink_and_format(dds_sex_fnf,
coef = "Sex_M_vs_F",
gse = gse_fnf)

load("data/dds_sex_fnf.rda")
results(dds_sex_fnf)




get_gene_age_Counts <- function(gene, dds){
# Get coef matrix and design matrix for fitted spline
coef_mat <- coef(dds)
design_mat <- model.matrix(design(dds), colData(dds))
geneCounts <- plotCounts(dds, gene = gene,
intgroup = "Age",
normalized = TRUE,
returnData = TRUE) %>%
remove_rownames() %>%
mutate(gene_id = gene,
logmu = design_mat %*% coef_mat[gene,])
return(geneCounts)
}