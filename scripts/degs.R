
# Libraries
require(pacman)
pacman::p_load(
  tidyverse,
  DESeq2,
  rgl,
  BiocParallel,
  reshape2,
  scatterplot3d,
  cowplot,
  patchwork,
  gridExtra,
  ggrepel,
  org.Hs.eg.db,
  clusterProfiler,
  AnnotationDbi
)

# read count_table
rawCounts <- read.table("./bulk/data/mapped_trimmed/count_table.txt", header=T)[,-c(2:6)] %>%
  rename_with(~str_remove(., '_R1_001_trimmed.fq.gz.sam_sorted.bam')) 

names(rawCounts) <- gsub('P4_coinf_16_S12','P4_PR8RSV_S12',colnames(rawCounts))

# Sample information
Run <- names(rawCounts)[2:length(names(rawCounts))] %>% sort()
treatment <-rep(c("Mock", "Flu", "Flu-RSV", "RSV"), 4)
replicate <- rep(c(1:4), 4) %>% sort()
#incubate_on_ice <- rep(c("Yes", "No"), 8) %>% sort(decreasing=T)

sampleData <- data.frame(Run=Run,virusType=treatment, replicateID=as.factor(replicate))
rownames(sampleData) <- sampleData$Run
sampleData <- sampleData[,-c(1)]
head(sampleData)

                        
# Convert count data to a matrix of appropriate form for DEseq2
geneID <- rawCounts$Geneid
sampleIndex <- grepl("P\\d+", colnames(rawCounts))
rawCounts <- as.matrix(rawCounts[,sampleIndex])
rownames(rawCounts) <- geneID
head(rawCounts)

rm(list=setdiff(ls(), c("rawCounts",  "sampleData")))

# Put the columns of the count data in the same order as rows names of the sample mapping, 
  #then make sure it worked
rawCounts <- rawCounts[,unique(rownames(sampleData))]
all(colnames(rawCounts) == rownames(sampleData))

# Order the tissue types so that it is sensible and make sure the control sample 
  #is first: normal sample -> primary tumor -> metastatic tumor
sampleData$virusType <- factor(sampleData$virusType, 
                               levels=c("Mock", "Flu", "RSV", "Flu-RSV"))

# Create the DEseq2DataSet object
deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData,
                                     design= ~ virusType + replicateID)

# Register the number of cores to use
register(MulticoreParam(5))

#Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data <- DESeq(deseq2Data,parallel = TRUE)

# #remove noisy genes (those with rowsums < 50)
# keep <- rowSums(counts(deseq2Data)) >= 50
# deseq2Data <- deseq2Data[keep,]

summary(deseq2Data)

# Change Ensembl IDs to gene ids
gene_symbols <- mapIds(org.Hs.eg.db, rownames(counts(deseq2Data)), "SYMBOL", "ENSEMBL")
gene_symbols2 <- mapIds(org.Hs.eg.db, rownames(counts(deseq2Data)), "ENTREZID", "ENSEMBL")
  # Update DESeq2 object with gene symbols
rowData(deseq2Data)$gene <- as.data.frame(gene_symbols)[,1]
rowData(deseq2Data)$entrez_id <- as.data.frame(gene_symbols2)[,1]

# Write normalized counts to file
write.csv(counts(deseq2Data,normalized=TRUE), file="./bulk/data/deseq_normalised_counts.csv")

#make 2D PCA plot
# 'regularized log' transformation
rld<-rlog(deseq2Data) #rlog is less sensitive to size factors, can be an issue when size factors vary widely
# varianceStabilizingTransformation
rld2<-varianceStabilizingTransformation(deseq2Data)
# distsRL <- dist(t(assay(rld)))
# mat <- as.matrix(distsRL)
# rownames(mat) <- colnames(mat) <- with(colData(deseq2Data), paste(shortname, sep=" : "))
# plotPCA(rld, intgroup=c("virusType", "replicateID"),returnData=TRUE)
# plotPCA(rld, intgroup=c("virusType", "replicateID"),returnData=F) +
#   theme_bw()

pcaData <- plotPCA(rld2, intgroup=c("virusType", "replicateID"),returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca<-ggplot(pcaData, aes(PC1, PC2, color=virusType)) +
  geom_point(size=3) +
  # xlim(-100, 100)+
  # ylim(-100, 100)+
  #geom_text(aes(label=rownames(pcaData)))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw() +
  coord_fixed()

ggsave2(
  filename="./bulk/output/pca.png",
  plot = ggplot2::last_plot(),
  width = 7,
  height = 4,
  units = c("in"),
  dpi = 300)

head(pcaData)
resultsNames(deseq2Data)

#results(dds, contrast=c("condition","C","B")) meaning genes with logFC > 0 are overexpressed in C.

flu <- results(deseq2Data, contrast=c("virusType","Flu","Mock"))
rsv <- results(deseq2Data, contrast=c("virusType","RSV","Mock"))
coinf <- results(deseq2Data, contrast=c("virusType","Flu-RSV","Mock"))

#data for IPA
write.csv(as.data.frame(rsv), file="./bulk/data/ipa_rsv.csv")
write.csv(as.data.frame(flu), file="./bulk/data/ipa_flu.csv")
write.csv(as.data.frame(coinf), file="./bulk/data/ipa_coinf.csv")


# Coerce to a data frame
coerce_to_df <- function(df){
  as.data.frame(df) %>%
    dplyr::mutate(significant = case_when(
      #padj < 0.05 & !between(log2FoldChange, -1, 1) ~ "Significant",
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ 'NA'))
}
rsv_df <- coerce_to_df(rsv)
flu_df <- coerce_to_df(flu)
coinf_df <- coerce_to_df(coinf)

gene_id_dict<- as.data.frame(cbind(rownames(deseq2Data), rowData(deseq2Data)$gene, rowData(deseq2Data)$entrez_id)) %>%
  dplyr::mutate(V2 = coalesce(V2,V1))
  
names(gene_id_dict) <- c("ensembl", "gene_id", "entrez_id")

# Volcano plots #########
plot_volcano <- function(mydata){
  data_sig <- mydata %>%
    drop_na()  %>%
    dplyr::filter(significant!="NA") %>%
    slice_min(padj, n=10)
  data_sig <- merge(data_sig, gene_id_dict, by.x=0, by.y="ensembl")
  title=paste(deparse(substitute(mydata)))
  ggplot() + 
    geom_point(data = mydata, aes(x=log2FoldChange, y=-log10(padj), color=significant)) + 
    labs(title=title, x="Log2Fold", y="-Log10 p")+ 
    geom_vline(xintercept=1, linetype="dashed", color = "grey", linewidth=0.5) +
    geom_vline(xintercept=-1, linetype="dashed", color = "grey", linewidth=0.5) +
    geom_hline(yintercept=-log10(0.05),linetype="dashed", color = "grey", linewidth=0.8) +
    theme_bw() + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.position="none" #none
    ) +
    xlim(c(-20, 20)) + ylim(c(0, 300)) + 
    scale_colour_manual(
      values = c("Upregulated" = "tan2","Downregulated" = "steelblue2","NA" = "grey"),
      limits = c("Upregulated", "Downregulated"))+
    ggrepel::geom_label_repel(data = data_sig,
                              mapping = aes(x = log2FoldChange, y=-log10(padj), label = gene_id))
  
}  

a = plot_volcano(rsv_df)
b = plot_volcano(flu_df) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

c = plot_volcano(coinf_df) + 
  theme(legend.position="right",
        legend.text=element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) 

all_volcanos <- a+b+c
all_volcanos

ggsave2(
  filename="./bulk/output/volcanos.png",
  plot =all_volcanos,
  width = 15,
  height = 6,
  units = c("in"),
  dpi = 300)



# DEGs per condition #######
extract_DEGS <- function(mydata){
  mydata_degs= rownames(mydata %>%
                          drop_na() %>%
                          dplyr::filter(significant != 'NA'))
  new_df = data.frame(ensembl=unlist(mydata_degs, use.names = FALSE)) 
  new_df$virus<-deparse(substitute(mydata))
  new_df <- new_df %>% 
    dplyr::inner_join(., gene_id_dict, by="ensembl")
  return(new_df)
}

head(gene_id_dict)

tmp_a = extract_DEGS(rsv_df) 
tmp_b = extract_DEGS(flu_df) 
tmp_c = extract_DEGS(coinf_df)

tmp_d = rbind(tmp_a, tmp_b)
deg_table = rbind(tmp_d, tmp_c)

table(deg_table$virus)
######
nrow(coinf_df %>%
  drop_na() %>%
  dplyr::filter(significant != 'NA') %>%
  dplyr::filter(significant == 'Upregulated'))/ length(coinf_gene_ids)



# check dispersal estimates in data
    #ie plot per-gene dispersion estimates together with the fitted mean-dispersion relationship.
pdf(file = "./bulk/output/dispersal_estimates.pdf",  width = 7, height = 4)
plotDispEsts(deseq2Data)
dev.off()


# Venn diagram ######
# data

# list1<-c(tmp_a,tmp_b, tmp_c)
# list2<-list1[c(1,3,5)]
# names(list2) <- c('RSV', 'IAV', 'Coinf')

list1<- list(RSV = tmp_a[,1], IAV = tmp_b[,1], Coinf = tmp_c[,1])

library(ggvenn)
venn<-ggvenn(
  list1, 
  show_percentage = TRUE,
  fill_color = c("#0073C2FF", 'bisque', "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 5
)

venn

ggsave2(
  filename="./bulk/output/ggvenn.png",
  width = 5,
  height = 5,
  units = c("in"),
  dpi = 300)



# set diff
#setdiff(A,B)
## DEGs unique to each condition (not shared with another condition) ####
unique_degs_per_condition <- lapply(1:length(list1), function(n) setdiff(list1[[n]], 
                                                                         unlist(list1[-n])))
names(unique_degs_per_condition) <- c('rsv', 'flu', 'coinf')


## 3D Scatter plot ##### 
    #(how does expression of a given gene compare across conditions?)

# Get DEGs
degs_rsv <- extract_DEGS(rsv_df)[,1]
degs_flu <- extract_DEGS(flu_df)[,1]
degs_coinf <- extract_DEGS(coinf_df)[,1]

## Get shared DEGs #######
shared_degs <- Reduce(intersect, list(degs_rsv,degs_flu,degs_coinf))

#filter_shared_degs
filter_shared_degs <- function(mydata, col_name){
  output= mydata %>%
    dplyr::mutate(ensembl=rownames(mydata)) %>%
    dplyr::filter(rownames(mydata) %in% shared_degs) %>%
    dplyr::select(log2FoldChange, ensembl) %>%
    dplyr::rename(!!sym(col_name) := log2FoldChange)
  
  return(output)
  }

rsv_tmp <- filter_shared_degs(rsv_df, "rsv")
flu_tmp <- filter_shared_degs(flu_df, "flu")
coinf_tmp <- filter_shared_degs(coinf_df, "coinf")

# merge
library(plotly)

scatter_df <- rsv_tmp %>%
  inner_join(., flu_tmp, by ='ensembl') %>%
  inner_join(., coinf_tmp, by ='ensembl')%>%
  inner_join(., gene_id_dict, by ='ensembl') %>%
  dplyr::rename(
    RSV=rsv,
    IAV=flu,
    Coinfected=coinf)

#plot
fig <- plot_ly(scatter_df, x = ~RSV, y = ~IAV, z = ~Coinfected, color = I('steelblue3'),
               text = ~paste('Gene:', gene_id))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'L2F RSV'),
                                   yaxis = list(title = 'L2F IAV'),
                                   zaxis = list(title = 'L2F Coinfected')))

fig

  #export to html  
    
## 2D option (bubble) #####
### facet the L2f changes

plot_scatter<- function(data, x, y){
  # Identify genes positively expressed in one condition and negative expressed in the other
  mylist<- sign(data[[deparse(substitute(x))]]) == sign(data[[deparse(substitute(y))]])
  data$annotation <- mylist
  #return(data)
  
  ggplot(data=data, aes(x={{x}}, y={{y}})) + # colour by GO?
    geom_point(aes(size = {{x}}, color = {{x}}),alpha=0.7) +
    scale_color_viridis_c(guide = "legend") +
    #scale_size_continuous(range = c(-5, 15))
    theme_bw() +
    theme(legend.position="none") +
    geom_text_repel(data=data %>% filter(annotation=='FALSE'),
                    aes(x={{x}}, y={{y}},label=gene_id), size=3) +
    geom_vline(xintercept=0) +
    geom_hline(yintercept=0) +
    xlim(-3, 15) +
    ylim(-4, 14)
  }

# Spearman's correlation
a<-cor(scatter_df$IAV,scatter_df$RSV, method = "spearman") %>% round(., 2)
b<-cor(scatter_df$Coinfected,scatter_df$IAV, method = "spearman") %>% round(., 2)
c<-cor(scatter_df$Coinfected,scatter_df$RSV, method = "spearman") %>% round(., 2)

a <- plot_scatter(scatter_df, IAV, RSV) + annotate("text", x=12, y=2, label= paste0("Corr=", a))
b <- plot_scatter(scatter_df, IAV, Coinfected) + annotate("text", x=12, y=2, label= paste0("Corr=", b))
c <- plot_scatter(scatter_df, Coinfected, RSV) + annotate("text", x=12, y=2, label= paste0("Corr=", c))



bubbles<- a/b/c
bubbles


ggsave2(
  filename="./bulk/output/scatter.png",
  width = 8,
  height = 9,
  units = c("in"),
  dpi = 300)

# HeatMap #####
# # use Flu as the reference for significant genes
# df <- na.omit(as.data.frame(flu_df))
# 
# # top genes (or use p-value in place of basemean)
# df.top <- df[ (df$baseMean > 50) & (abs(df$log2FoldChange) > 1),]
# df.top <- df.top[order(df.top$log2FoldChange, decreasing = TRUE),]
# df.top
# 
# #get normalized count data from dds object
# rlog_out <- rlog(deseq2Data, blind=FALSE) 
# mat<-assay(rlog_out)[rownames(df.top), rownames(sampleData)] #sig genes x samples
# colnames(mat) <- rownames(sampleData)
# base_mean <- rowMeans(mat)
# mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
# colnames(mat.scaled)<-colnames(mat)
# 
# # keep top 25 and bottom 25
# num_keep <- 25
# #1 to num_keep len-num_keep to len
# rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)) )
# 
# #getting log2 value for each gene we are keeping
# l2_val <- as.matrix(df.top[rows_keep,]$log2FoldChange) 
# colnames(l2_val)<-"logFC"
# #getting mean value for each gene we are keeping
# mean <- as.matrix(df.top[rows_keep,]$baseMean) 
# colnames(mean)<-"AveExpr"
# 
# #BiocManager::install("ComplexHeatmap")
# library(ComplexHeatmap)
# library(circlize)
# 
# #maps values between b/w/r for min and max l2 values
# col_logFC <- colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue", "white", "red")) 
# #maps between 0% quantile, and 75% quantile of mean values --- 0, 25, 50, 75, 100
# col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))
# 
# 
# ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), 
#                                                height = unit(2, "cm")))
# h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = F, 
#               column_labels = colnames(mat.scaled), name="Z-score",
#               cluster_columns = T)
# 
# pdf(file = "./output/top50genes_flu_heatmap.pdf",  width = 9, height = 9)
# h1
# dev.off()
# 
# 
# 
# h2 <- Heatmap(l2_val, row_labels = df.top$symbol[rows_keep], 
#               cluster_rows = F, name="logFC", top_annotation = ha, col = col_logFC,
#               cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
#                 grid.text(round(l2_val[i, j],2), x, y)
#               })
# h3 <- Heatmap(mean, row_labels = df.top$symbol[rows_keep], 
#               cluster_rows = F, name = "AveExpr", col=col_AveExpr,
#               cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
#                 grid.text(round(mean[i, j],2), x, y)
#               })
# h<-h1+h2+h3
# h
# 
# # all genes 
# # use Flu as the reference for significant genes
# df <- as.data.frame(flu_df)
# #get normalized count data from dds object
# rlog_out <- rlog(deseq2Data, blind=FALSE) 
# mat<-assay(rlog_out)[rownames(df), rownames(sampleData)] #sig genes x samples
# colnames(mat) <- rownames(sampleData)
# # order by flu
# 
# 
# 
# mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
# colnames(mat.scaled)<-colnames(mat)
# Heatmap(mat.scaled, cluster_rows = F, 
#         column_labels = colnames(mat.scaled), name="Z-score",
#         km=T,# Apply k-means clustering on rows
#               #heatmap will be split by rows according to the k-means clustering
#         row_km_repeats=1000,
#         show_row_names=F, # hide gene names as too many
#         #cluster_columns = T
#         )

# Shared genes (Venn intersect)
shared_degs
rlog_out <- rlog(deseq2Data, blind=FALSE) 
  #filter for genes after transformation 
mat<-assay(rlog_out)[shared_degs, rownames(sampleData)]
colnames(mat) <- rownames(sampleData)
base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled)<-colnames(mat)
#getting log2 value for each gene we are keeping
# a<-as.data.frame(base_mean) 
# a$gene<- rownames(a)
# a[order(a$base_mean, decreasing = TRUE),]
# gene_order_rowmeans<-a$gene  
# 
# heat_map_order_by_flu<-flu_df[order(flu_df$log2FoldChange, decreasing = TRUE),]
# l2_val <- as.matrix(heat_map_order_by_flu[shared_degs,]$log2FoldChange) 
# colnames(l2_val)<-"logFC"
# col_logFC <- colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue", "white", "red")) #maps values between b/w/r for min and max l2 values

intersect_heat<-Heatmap(mat.scaled[shared_degs,], cluster_rows = T, 
                        #row_order = gene_order_rowmeans,
                        column_labels = colnames(mat.scaled), name="Z-score",
                        cluster_columns = T, show_row_names=F,)

intersect_heat

pdf(file = "./output/intersect_genes_heatmap.pdf",  width = 9, height = 9)
intersect_heat
dev.off()


# # using geom_tile instead of heatmap
# tile_data<- as.data.frame(mat.scaled) %>%
#   dplyr::mutate(geneID=rownames(mat.scaled))%>%
#   pivot_longer(!geneID, names_to = "sample", values_to = "Z_score")
# 
# sample_order <- c('P1_mock_S1', 'P2_mock_S2', 'P3_mock_S9', 'P4_mock_10_S9', 'P1_RSV_S5', 
#                   'P2_RSV_S6','P3_RSV_S11', 'P4_RSV_14_S11', 'P1_PR8RSV_S7', 'P2_PR8RSV_S8', 
#                   'P3_PR8RSV_S12', 'P4_PR8RSV_S12', 'P1_PR8_S3', 'P2_PR8_S4', 'P3_PR8_S10',
#                   'P4_PR8_12_S10')
# 
# coinf_tmp <- filter_shared_degs(coinf_df, "coinf")
# gene_order<- rownames(coinf_tmp[order(coinf_tmp$coinf, decreasing = TRUE),])
# 
# ggplot(tile_data, aes(x=factor(sample, level=sample_order), 
#                       y=factor(geneID, level=gene_order), fill=Z_score)) +
#   geom_tile()+
#   scale_fill_gradient2(low = "#075AFF",
#                        mid = "#FFFFCC",
#                        high = "#FF0000")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         axis.text.y = element_blank())




# heatmap using L2FC

      # tmp_rsv<- rsv_df %>%
      #   dplyr::select(log2FoldChange)
      # tmp_flu<- flu_df %>%
      #   dplyr::select(log2FoldChange) 
      # tmp_coinf<- coinf_df %>%
      #   dplyr::select(log2FoldChange) %>%
      #   dplyr::mutate(geneID=rownames(coinf_df))
      # tmp_rsv_flu<- merge(tmp_rsv, tmp_flu, by=0)
      # names(tmp_rsv_flu) <- c('geneID', 'rsv', 'flu')
      # tile_data_l2fc <- merge(tmp_rsv_flu, tmp_coinf, by='geneID')
      # 
      # names(tile_data_l2fc) <- c('geneID', 'rsv', 'flu', 'coinf')
      # names(tile_data_l2fc) <- c('geneID', 'cond', 'L2F')
      # 
      # head(tile_data_l2fc)
        # 
        # tile_data_l2fc <-  tile_data_l2fc %>%
        #   pivot_longer(!geneID, names_to = "sample", values_to = "L2F")

head(scatter_df)

gene_order <- scatter_df[order(scatter_df$Coinfected),][,2]

tile_data<- gather(scatter_df, condition, L2F, c(RSV, IAV, Coinfected), factor_key=TRUE)

range(tile_data$L2F)

tmp <- coinf_tmp %>%
  dplyr::inner_join(., gene_id_dict, by="ensembl") 

gene_order<-tmp[order(tmp$coinf, decreasing = TRUE),][['gene_id']]

intersect_l2f<- ggplot(tile_data, aes(x=factor(condition, 
                               level=c('RSV', 'IAV', 'Coinfected')), 
                      y=factor(gene_id, level=gene_order), fill=L2F)) +
  geom_tile()+
  #scale_fill_gradientn(colours = c("cyan", "white", "green","brown1", "red"))
  scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue')+
  theme(axis.text.y = element_blank())+
  ylab("Genes")+ xlab("Condition")
  


intersect_l2f


pdf(file = "./bulk/output/l2f_intersect_genes_geom_tile.pdf",  width = 9, height = 9)
intersect_l2f
dev.off()


# GSEA gene set enrichment analysis (GSEA) #######
    # BiocManager::install("org.Hs.eg.db") 
    # BiocManager::install("clusterProfiler")
    # BiocManager::install("AnnotationDbi")



# order data
get_genes<- function(mydata){
  # # filter out noisy genes
  # mydata<- mydata[mydata$baseMean > 50,]
  # order using stat. You can use L2FC
  mydata <- mydata[order(-mydata$log2FoldChange),] %>%
    dplyr::filter(significant !="NA")
  # get gene list
  gene_list <- mydata$log2FoldChange
  names(gene_list) <- rownames(mydata)
  return(gene_list)
}

rsv_gene_ids <- names(get_genes(rsv_df))
flu_gene_ids <- names(get_genes(flu_df))
coinf_gene_ids <- names(get_genes(coinf_df))
shared_degs

# # Covert geneIDs to ensemblIDs and/or EntrezIds using the org.Hs.eg.db package
# # Specify the human gene IDs that you want to convert
# convert_geneIDs <- function(gene_ids){
#   ensembl_ids <- mapIds(org.Hs.eg.db, 
#                             keys = gene_ids, 
#                             column = "ENSEMBL", 
#                             keytype = "SYMBOL") %>% as.data.frame()
#   entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "ENTREZID",
#                        keytype = "SYMBOL") %>% as.data.frame()
#   combined_df <- merge(ensembl_ids, entrez_ids, by=0, all=TRUE)
#   names(combined_df) <- c("gene_ids","ensembl_ids", "entrez_ids")
#   return(combined_df)
#   }
# 
# rsv_genes_df<- convert_geneIDs(rsv_gene_ids)


      # gse <- gseGO(geneList=sort(rsv_genes_df[,2], decreasing = T),
      #              ont = "BP",
      #              keyType = "ENSEMBL",
      #              OrgDb = "org.Hs.eg.db")
      # 
      # 
      # library(enrichR)
      # # Perform gene set enrichment analysis using Enrichr
      # result <- enrichR::enrichr(genes = sort(rsv_genes_df[,3], decreasing = T), databases ="GO_Biological_Process_2018")

# None of these two approaches map :(

library(enrichplot)
library(ggupset)

# GroupGO: Functional Profile of a gene set at specific GO level. 
    # Given a vector of genes, this function will return the GO profile at a specific level.

rsv_ggo <- groupGO(gene     = rsv_gene_ids,
               OrgDb    = org.Hs.eg.db,
               keyType =  "ENSEMBL",         # "ENTREZID",
               ont      = "BP",
               readable = TRUE)
flu_ggo <- groupGO(gene     = flu_gene_ids,
                   OrgDb    = org.Hs.eg.db,
                   keyType =  "ENSEMBL",         # "ENTREZID",
                   ont      = "BP",
                   readable = TRUE)
coinf_ggo <- groupGO(gene     = coinf_gene_ids,
                   OrgDb    = org.Hs.eg.db,
                   keyType =  "ENSEMBL",         # "ENTREZID",
                   ont      = "BP",
                   readable = TRUE)
View(as.data.frame(flu_ggo))


# GO Ontology ###########
# enrichGO: GO Enrichment Analysis of a gene set. 
    #Given a vector of genes, this function will return the enrichment GO categories after FDR control.

# Note: biological_process != pathway. 
 # Specifically it doesnt represent any of the dynamics or dependencies that would be required to describe a pathway.
 # High-level processes such as ‘cell death’ [GO:0008219] can have both subtypes, such as ‘apoptotic process’ [GO:0006915], and subprocesses, such as ‘apoptotic chromosome condensation’ [GO:0030263]

rsv_go_enrich <- enrichGO(gene = rsv_gene_ids, 
                       OrgDb = "org.Hs.eg.db", 
                       keyType = "ENSEMBL", # or "ENTREZID" 
                       ont = "BP") #BP, MF, and CC represent Biological Process, Molecular Function, and Cellular Component groups of GO

flu_go_enrich <- enrichGO(gene = flu_gene_ids, 
                          OrgDb = "org.Hs.eg.db", 
                          keyType = "ENSEMBL", # or "ENTREZID" 
                          ont = "BP")

coinf_go_enrich <- enrichGO(gene = coinf_gene_ids, 
                          OrgDb = "org.Hs.eg.db", 
                          keyType = "ENSEMBL", # or "ENTREZID" 
                          ont = "BP")

intersect_go_enrich <- enrichGO(gene = shared_degs, 
                            OrgDb = "org.Hs.eg.db", 
                            keyType = "ENSEMBL", # or "ENTREZID" 
                            ont = "BP")


View(as.data.frame(rsv_go_enrich) ) 

barplot(rsv_go_enrich, showCategory = 20)
upsetplot(rsv_go_enrich)
emapplot(rsv_go_enrich)

pdf(file = "./output/rsv_GO_enrich_upsetplot.pdf",  width = 14, height = 9)
upsetplot(rsv_go_enrich)
dev.off()

pdf(file = "./output/flu_GO_enrich_upsetplot.pdf",  width = 14, height = 9)
upsetplot(flu_go_enrich)
dev.off()

original_gene_list <- df$log2FoldChange
# omit any NA values 
gene_list<-na.omit(rsv_df$log2FoldChange)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
cnetplot(rsv_go_enrich, categorySize="pvalue", foldChange=gene_list)

# REACTOME #####
# BiocManager::install("ReactomePA")
library(ReactomePA)

flu_df




