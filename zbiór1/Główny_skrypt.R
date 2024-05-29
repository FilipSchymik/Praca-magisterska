############################ BIBLIOTEKI ########################################
library("readr")
library("DESeq2")
library("dplyr")
library("reshape2")
library("tidyr")
library("RColorBrewer")
library('ggplot2')
library("dplyr")
library("EnhancedVolcano")
library("clusterProfiler")
library("org.Hs.eg.db")
library("edgeR")
library(msigdbr)

############################ METADANE ##########################################

# podajemy który to dataset
dataset <- "dataset_1"

# ustalamy katalog roboczy
dir <- "~/Desktop/Analiza_2.0" 

# ustalamy katalog do zapisywania macierzy zliczeń
matrix_dir <- "~/Desktop/znormalizowane_macierze_zliczeń"

# ustalamy katalog do zapisywania wykresów
plot_dir <- "~/Desktop/master_plots"

# wczytujemy plik z informacjami o próbkach
metadata <- read.csv(file.path(dir, "experiment_metadata.csv"))

############################ WCZYTANIE DANYCH SALMON ###########################

# wczytujemy zliczenia z salmona
salmon_counts <- read.table(
  file = file.path(dir, "nfcore/Salmon_rnaseq/salmon/salmon.merged.gene_counts.tsv"),
  sep = "\t", 
  header = TRUE)

# usuwamy gene name
salmon_counts$gene_name <- NULL

# ustawiamy kolumnę gene_id jako nazwy wierszy
rownames(salmon_counts) <- salmon_counts$gene_id
salmon_counts$gene_id <- NULL

############################ ANALIZA RÓŻNICOWA SALMON ##########################

# tworzymy obiekt DESeqDataSet na podstawie kontrola / depresja dla salmon
dds_salmon <- DESeqDataSetFromMatrix(countData = round(salmon_counts),
                              colData = metadata,
                              design = ~ Group)

# ocenimy średnią głębokość sekwencjonowania
counts_for_depth <- counts(dds_salmon, normalized=FALSE)
average_depth <- mean(colSums(counts_for_depth))
print(average_depth)

# Wyliczymy CPM
filter_counts_salmon <- counts(dds_salmon, normalized=FALSE)
cpm_values_salmon <- cpm(filter_counts_salmon)

# Filtracja- 0.1 CPM w 60% próbek
filt <- rowSums(cpm_values_salmon >= 0.1) >= (0.6 * ncol(filter_counts_salmon))
dds_salmon <- dds_salmon[filt,]

# Analiza różnicowa
dds_salmon <- DESeq(dds_salmon)

# Wykres dyspersji
plotDispEsts(dds_salmon)

############################ WCZYTANIE DANYCH GEO ##############################

# wczytujemy zliczenia z geo
geo_counts <- read.table(
  file = file.path(dir, "geo_counts.tsv"),
  sep = "\t", 
  header = TRUE)

# ustawiamy kolumnę gene_id jako nazwy wierszy
rownames(geo_counts) <- geo_counts$gene_id
geo_counts$gene_id <- NULL

# znajdowanie wierszy specjalnych
special_rows <- grep("^__", rownames(geo_counts), value = TRUE)

# usuwamy wiersze specjalne
geo_counts <- geo_counts[!rownames(geo_counts) %in% special_rows, ]

############################ ANALIZA RÓŻNICOWA GEO #############################

# tworzymy obiekt DESeqDataSet na podstawie kontrola / depresja dla geo
dds_geo <- DESeqDataSetFromMatrix(geo_counts, colData = metadata, design = ~ Group)

# Wyliczymy CPM
filter_counts_geo <- counts(dds_geo, normalized=FALSE)
cpm_values_geo <- cpm(filter_counts_geo)

# Filtracja- 0.1 CPM w 60% próbek
filt <- rowSums(cpm_values_geo >= 0.1) >= (0.6 * ncol(filter_counts_geo))
dds_geo <- dds_geo[filt,]

# Analiza różnicowa
dds_geo <- DESeq(dds_geo)

############################ ANALIZA GENE ID ###################################

# wczytujemy znormalizowane zliczenia z salmona
zliczenia_salmon <- counts(dds_salmon, normalized=TRUE)

# zapiszemy znormalizowane zliczenia z salmona na komputerze 
write.csv(as.data.frame(zliczenia_salmon), file = file.path(
  matrix_dir, dataset, "znormalizowane_zliczenia_salmon.csv"))

# wczytujemy znormalizowane zliczenia z geo
zliczenia_geo <- counts(dds_geo, normalized=TRUE)

# zapiszemy znormalizowane zliczenia z geo na komputerze 
write.csv(as.data.frame(zliczenia_geo), file = file.path(
  matrix_dir, dataset, "znormalizowane_zliczenia_geo.csv"))

# znajdujemy wspólne gene_id obecne i w geo i w salmon
geny_wspólne <- intersect(rownames(zliczenia_salmon), rownames(zliczenia_geo))

# znajdujemy wszystkie geny obecne w obydwu zestawach
geny_wszystkie <- union(rownames(zliczenia_salmon), rownames(zliczenia_geo))

# znajdujemy geny wskazane tylko przez salmon
unikatowe_salmon <- setdiff(rownames(zliczenia_salmon), rownames(zliczenia_geo))

# znajdujemy geny wskazane tylko przez geo
unikatowe_geo <- setdiff(rownames(zliczenia_geo), rownames(zliczenia_salmon))

# obliczamy jaki ułamek wszystkich genów stanowią geny wspólne
wspólne_wobec_wszystkich <- length(geny_wspólne) / length(geny_wszystkie)

# obliczamy jaki ułamek unikatowych genów stanowią geny unikatowe z salmon
unikatowe_salmon_wobec_wszystkich_unikatowych <- length(unikatowe_salmon) / 
  (length(unikatowe_salmon) + length(unikatowe_geo))

# obliczamy jaki ułamek unikatowych genów stanowią geny unikatowe z geo
unikatowe_geo_wobec_wszystkich_unikatowych <- length(unikatowe_geo) / 
  (length(unikatowe_salmon) + length(unikatowe_geo))

# definiujemy kolumne opisującą co wyliczamy
kolumna_proporcji <- c("geny w salmon", "geny w geo", "suma genów", 
                       "geny wspólne", "geny unikatowe dla salmon",
                       "geny unikatowe dla geo", "wspóle wobec wszystkich",
                       "unikatowe salmon wobec wszystkich unikatowych",
                       "unikatowe geo wobec wszystkich unikatowych")

# definiujemy kolumnę z wynikami
kolumna_wartości <- c (length(rownames(zliczenia_salmon)), 
                       length(rownames(zliczenia_geo)),
                       length(geny_wszystkie), length(geny_wspólne), 
                       length(unikatowe_salmon), 
                       length(unikatowe_geo), wspólne_wobec_wszystkich,
                       unikatowe_salmon_wobec_wszystkich_unikatowych,
                       unikatowe_geo_wobec_wszystkich_unikatowych)

# zapisujemy jako ramka danych
macierz_proporcji <- data.frame(kolumna_proporcji, kolumna_wartości)

# zapisujemy na dysku
write.csv(macierz_proporcji, file.path(plot_dir, dataset, "macierz_proporcji.csv"), row.names = FALSE)

############################ KORELACJA ZLICZEŃ #################################

# dodamy macierze dla średnich zliczeń Salmon i GEO tylko z wspólnymi genami. 
common_counts_salmon <- log2(rowMeans(zliczenia_salmon[geny_wspólne, , drop=FALSE], na.rm = TRUE) + 1)
common_counts_geo <- log2(rowMeans(zliczenia_geo[geny_wspólne, , drop=FALSE], na.rm = TRUE) + 1)

# połączymy je w jeden df
data_for_correlation <- data.frame(
  salmon = common_counts_salmon,
  geo = common_counts_geo
)

# Wykres z logarytmicznie transformowanymi zliczeniami dla genów wspólnych(znoramlizowanymi przez DESeq2)

png(file.path(plot_dir, dataset, "spearman_corr.png"), width=800, height=600)
par(mar = c(5, 5, 2, 2), mgp=c(2.2, 0.5, 0))
plot(data_for_correlation$salmon, data_for_correlation$geo, 
     xlab = expression(log[2](N[s])),  
     ylab = expression(log[2](N[g])),     
     xlim = c(0, 20), ylim = c(0, 20),
     cex.lab = 1.4, cex.axis = 1.4)   # Zwiększa rozmiar etykiet osi
abline(lm(data_for_correlation$geo ~ data_for_correlation$salmon), col = "red") # Dodaje linię regresji
# Zwiększa rozmiar tekstu dla wyniku korelacji Spearmana i dodaje kolor
text(x = 6, y = 18, labels = paste("Korelacja Spearmana =", round(
  cor(data_for_correlation$salmon, data_for_correlation$geo, 
      method = "spearman"), 2)), col = "blue", cex = 1.5)
dev.off()

############################ BOX PLOTS #########################################

# wczytujemy NIEznormalizowane zliczenia z salmona
unnormalized_zliczenia_salmon <- counts(dds_salmon, normalized=FALSE)

# wczytujemy NIEznormalizowane zliczenia z geo
unnormalized_zliczenia_geo <- counts(dds_geo, normalized=FALSE)

# nazwy plików
file_names <- c("boxplot_unnormalized_salmon.png", 
                "boxplot_unnormalized_geo.png", 
                "boxplot_normalized_salmon.png",
                "boxplot_normalized_geo.png")

# Lista danych do wykresów
data_list <- list(
  log2(unnormalized_zliczenia_salmon + 1),
  log2(unnormalized_zliczenia_geo + 1),
  log2(zliczenia_salmon + 1),
  log2(zliczenia_geo + 1)
)

# Pętla przez wszystkie dane i nazwy plików
for (i in 1:length(data_list)) {
  png(file.path(plot_dir, dataset, file_names[i]), width=800, height=600)  # Otwarcie pliku PNG
  par(cex.axis=0.7)
  boxplot(data_list[[i]], las=2, ylim=c(0, 26))  # Generowanie boxplota
  dev.off()  # Zamknięcie pliku PNG
}

########################### ANALIZA RÓŻNICOWA ##################################

# wczytamy wyniki z DESeq2 dla salmon
res_salmon <- results(dds_salmon, alpha = 0.05)

# i zapiszemy na komputerze
write.csv(res_salmon, file.path(matrix_dir, dataset, "salmon_results_DESeq2.csv"),
          row.names = TRUE)

# wczytamy wyniki z DESeq2 dla geo
res_geo <- results(dds_geo, alpha = 0.05)

# i zapiszemy na komputerze
write.csv(res_geo, file.path(matrix_dir, dataset, "geo_results_DESeq2.csv"),
          row.names = TRUE)

# wykres wulkaniczny dla salmon z użyciem standardowego p-value
EnhancedVolcano(res_salmon,
                lab = rownames(res_salmon),
                x = 'log2FoldChange',
                y = 'pvalue')

########################### WZBOGACENIE GENÓW ##################################

#### dla Salmona ####
# we want the log2 fold change 
original_gene_list_salmon <- res_salmon$log2FoldChange
names(original_gene_list_salmon) <- rownames(res_salmon)

# omit any NA values 
original_gene_list_salmon <-na.omit(original_gene_list_salmon)

# sort the list in decreasing order (required for clusterProfiler)
gene_list_salmon = sort(original_gene_list_salmon, decreasing = TRUE)

gse <- gseGO(geneList=gene_list_salmon, 
             ont ="BP", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "BH")

require(DOSE)

png(file.path(plot_dir, dataset, "biological_process_salmon.png"), width=800, height=600)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()

#### dla GEO ####
# we want the log2 fold change 
original_gene_list_geo <- res_geo$log2FoldChange
names(original_gene_list_geo) <- rownames(res_geo)

# omit any NA values 
gene_list_geo <-na.omit(original_gene_list_geo)

# sort the list in decreasing order (required for clusterProfiler)
gene_list_geo = sort(gene_list_geo, decreasing = TRUE)

gse_geo <- gseGO(geneList=gene_list_geo, 
                 ont ="BP", 
                 keyType = "ENSEMBL", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = org.Hs.eg.db, 
                 pAdjustMethod = "BH")

require(DOSE)

png(file.path(plot_dir, dataset, "biological_process_geo.png"), width=800, height=600)
dotplot(gse_geo, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()