# RAAK counts ------------
raak_genes <- read_csv("../DE/data/RAAK/RAAK_genes.csv")
# Import non- normalized data and convert to DESeq?
# filter for OA/preserved
sample_names <- read_delim("../DE/data/RAAK/GSE57218_series_matrix.txt.gz", 
                           delim = "!", skip = 37, n_max = 1, col_names = "samples") |> 
  pull(X2) |> 
  str_split("\t") |> 
  unlist()

sample_names <- grep("cartilage", sample_names, value = TRUE)
sample_names <- gsub("\\\"", "", sample_names)

col_sample_names <- paste0(rep(sample_names, each = 2), c("_expression", "_detectionP"))

non_norm_raak <- read_delim("../DE/data/RAAK/GSE57218_Non-normalized_data.txt.gz",
                            col_names = c("ProbeID", "ID_REF", col_sample_names), skip = 1) |> 
  # Get rid of healthy sample columns
  dplyr::select(-contains("Healthy")) |> 
  # Just pull out expression columns
  dplyr::select(ProbeID, ID_REF, ends_with("expression"))

# Convert Illumina IDs to ENSEMBL gene ids and symbols
illumina_genes <- AnnotationDbi::select(illuminaHumanv4.db, keys = non_norm_raak$ID_REF, 
                                        columns = c("SYMBOL"),
                                        keytype = "PROBEID")
non_norm_raak <- non_norm_raak |> 
  left_join(illumina_genes, by = join_by("ID_REF" == "PROBEID"), relationship = "many-to-many")

# Define colData
colData <- non_norm_raak |> 
  dplyr::select(starts_with("cartilage")) |> 
  colnames() |> 
  as_tibble_col(column_name = "Sample") |> 
  separate_wider_delim(cols = "Sample", delim = "_", names = c(NA, NA, "Sample_no", "Disease_state", NA)) |> 
  mutate(Disease_state = factor(Disease_state, levels = c("Preserved", "OA")))

# Make count matrix
counts <- non_norm_raak |> 
  dplyr::select("SYMBOL", starts_with("cartilage")) |> 
  distinct(SYMBOL, .keep_all = TRUE) |> 
  filter(!is.na(SYMBOL)) |> 
  column_to_rownames(var = "SYMBOL") |> 
  as.matrix() 
counts_integer <- apply(counts, 2, as.integer)
rownames(counts_integer) <- rownames(counts)

# Make DDS
raak_dds <- DESeqDataSetFromMatrix(counts_integer, colData, design = ~Sample_no + Disease_state) 
raak_dds <- DESeq(raak_dds)


DESeq2::plotCounts(raak_dds, gene = "PAPPA", intgroup = "Disease_state", 
                   normalized = TRUE,
                   returnData = TRUE)