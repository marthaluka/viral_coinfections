
pacman::p_load(viridis,
               enrichplot,
               DOSE,
               ggnewscale)


# Summarise GO Ontologies and pathway analyses
rsv_go_enrich
flu_go_enrich
coinf_go_enrich
intersect_go_enrich

# Number of GO ontologies per condition
nrow(as.data.frame(rsv_go_enrich)) 
nrow(as.data.frame(flu_go_enrich)) 
nrow(as.data.frame(coinf_go_enrich)) 
nrow(as.data.frame(intersect_go_enrich))
#
nrow(as.data.frame(intersect_go_enrich) %>%
  dplyr::filter(ID %in% shared_GOs))


# Merge the datasets 
merge_dfs <- function(){
  tmp <- rbind(as.data.frame(rsv_go_enrich) %>%
                 dplyr::mutate(virus = "RSV"), as.data.frame(flu_go_enrich) %>%
                 dplyr::mutate(virus = "IAV")) %>%
    rbind(., as.data.frame(intersect_go_enrich)%>%
            dplyr::mutate(virus = "Intersect_genes"))
    
  merged_df <- rbind(tmp, as.data.frame(coinf_go_enrich)%>%
                       dplyr::mutate(virus = "Coinf"))
  
  dataset<-separate_wider_delim(merged_df, GeneRatio, delim ="/", names =c("count", "tot_genes")) %>%
    separate_wider_delim(., BgRatio, delim ="/", names =c("db_count", "db_tot_genes")) %>%
    dplyr::mutate(
      geneRatio = as.numeric(count)/as.numeric(tot_genes),
      bgRatio = as.numeric(db_count)/as.numeric(db_tot_genes)
    )
  return(dataset)
}


# plot most common ontologies
plot_common_ontologies <- function(virusType, colName){
  dataset <- merge_dfs() %>%
    dplyr::filter(virus == virusType) %>%
    mutate(Description=replace(Description, Description=="adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains", "adaptive immune response"))

  plot_data <- as.data.frame(dataset) %>% 
    dplyr::slice_max(!! sym(colName), n=30) # we plot top 30
  
  plot_data$Description <- factor(plot_data$Description, levels = plot_data$Description[order(plot_data$Count)])

  # bar_chat <- ggplot(plot_data)+
  #   geom_bar(aes(x=Description, y=plot_data[,colName]), fill="darkorchid",  stat = "identity")+ 
  #   coord_flip() +
  #   theme_bw()
  
  bubble_plot<-ggplot(plot_data)+
    geom_point(aes(x=Description, y=plot_data[,colName], color=p.adjust,  size = Count))+ 
    coord_flip() +
    theme_bw() +
    ylab("Gene ratio")+
    ggtitle(virusType)
  return(bubble_plot)
  
  # if (colName == "Count"){
  #   return(bar_chat)
  # } else {
  #   return(bubble_plot)
  # }
}


a <- plot_common_ontologies("RSV", "geneRatio")
b <- plot_common_ontologies("IAV", "geneRatio")
c <- plot_common_ontologies("Coinf", "geneRatio")
d <- plot_common_ontologies("Intersect_genes", "geneRatio")
(a+b)/(c+d)

ggsave2(
  filename="./bulk/output/go_ontologies.png",
  plot = ggplot2::last_plot(),
  width = 13,
  height = 10,
  units = c("in"),
  dpi = 300)



# Plot overlapping GOs in the three conditions in one bubble plot
list2<- list(RSV = as.data.frame(rsv_go_enrich)[,1], IAV = as.data.frame(flu_go_enrich)[,1], 
             Coinf = as.data.frame(coinf_go_enrich)[,1])

ggvenn(
  list2, 
  show_percentage = TRUE,
  fill_color = c("#0073C2FF", 'bisque', "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 5
)


## GOs unique to each condition (not shared with another condition) ####
unique_GOs_per_condition <- lapply(1:length(list2), function(n) setdiff(list2[[n]], 
                                                                         unlist(list2[-n])))
names(unique_GOs_per_condition) <- c('rsv', 'flu', 'coinf')

## Get shared DEGs #######
shared_GOs <- Reduce(intersect, list2)

merged_df %>%
  dplyr::filter(ID %in% shared_GOs) %>%
  ggplot()+
  geom_point(aes(x=Description, y=geneRatio, color=virus,  size = Count))+ 
  coord_flip() +
  theme_bw() +
  ylab("Gene ratio") 

###########
coinf_gene_ids
flu_gene_ids
rsv_gene_ids


# Gene concept networks

gene_network <- function(gene_ids, condition_df){
  
  # convert from ensembl to gene IDs using a pre-made dictionary
  entrez_ids <- gene_id_dict %>%
    dplyr::filter(ensembl %in% gene_ids) %>%
    drop_na()
  
  # Enrichment analysis based on the DisGeNET 
  edo <- enrichDGN(entrez_ids[,3])
  
  # create a geneList (geneIDs and L2FCs)
   #github.com/YuLab-SMU/DOSE/wiki/how-to-prepare-your-own-geneList
  ## feature 1: numeric vector
  tmp<-condition_df %>% 
    dplyr::mutate(ensembl = rownames(.)) %>%
    inner_join(., gene_id_dict, by="ensembl") %>%
    drop_na()
  geneList = tmp[,"log2FoldChange"]
  ## feature 2: named vector
  names(geneList) = as.character(tmp[,"entrez_id"])
  ## feature 3: decreasing order
  geneList = sort(geneList, decreasing = TRUE)
  
  # convert gene ID to Symbol
  edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
  
  #plot
  cnet <- cnetplot(edox, foldChange=geneList2,
                node_label="all", cex_label_gene =.5,
                cex_label_category=.7)
  
  return(cnet)
  
}

a<- gene_network(rsv_gene_ids, rsv_df)
b<- gene_network(flu_gene_ids, flu_df)
c<- gene_network(coinf_gene_ids, coinf_df)

# p3 <- cnetplot(edox, foldChange=geneList2, circular = TRUE, colorEdge = TRUE) 
cowplot::plot_grid(a, b, c, ncol=2, labels=LETTERS[1:3], rel_widths=c(.8, 1, 1))


(a+b) / c

ggsave2(
  filename="./bulk/output/gene_networks.png",
  plot = ggplot2::last_plot(),
  width = 10,
  height = 8,
  units = c("in"),
  dpi = 300)



# Heatmap-like functional classification
p1 <- heatplot(edox, showCategory=10)
p2 <- heatplot(edox, foldChange=geneList2, 
               showCategory=15)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])


# 
# edox2 <- pairwise_termsim(edox)
# p1 <- treeplot(edox2)
# p2 <- treeplot(edox2, hclust_method = "average")
# aplot::plot_list(p1, p2, tag_levels='A')


# Enrichment Map
edo_ema <- pairwise_termsim(edo)
p1 <- emapplot(edo_ema)



# Merge the datasets 
tmp <- rbind(as.data.frame(rsv_go_enrich) %>%
  dplyr::mutate(virus = "RSV"), as.data.frame(flu_go_enrich) %>%
  dplyr::mutate(virus = "IAV"))
merged_df <- rbind(tmp, as.data.frame(coinf_go_enrich)%>%
  dplyr::mutate(virus = "Coinf"))


# Disease Gene Set Enrichment Analysis

## feature 1: numeric vector
tmp<-rsv_df %>% 
  dplyr::mutate(ensembl = rownames(.)) %>%
  inner_join(., gene_id_dict, by="ensembl") %>%
  drop_na()
geneList2 = tmp[,"log2FoldChange"]
## feature 2: named vector
names(geneList2) = as.character(tmp[,"entrez_id"])
## feature 3: decreasing orde
geneList2 = sort(geneList2, decreasing = TRUE)


edo2 <- gseDO(geneList2)
dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")

coinf<-ridgeplot(edo2) + ggtitle("Coinf")

ridgeplots<-(rsv+flu)/coinf
ggsave2(
  filename="./bulk/output/ridgeplots.png",
  plot = ggplot2::last_plot(),
  width = 16,
  height = 9,
  units = c("in"),
  dpi = 300)


# rank plot

df <- merged_df %>% 
  dplyr::select(ID, Description, geneRatio, virus)%>%
  dplyr::group_by(virus)%>%
  arrange(desc(geneRatio))%>% 
  mutate(ranking = row_number())%>%
  dplyr::ungroup()%>%
  as.data.frame()%>%
  dplyr::filter(ranking<50)


ggplot(data = df, aes(x = virus, y = ranking, group = Description)) +
  geom_line(aes(color = Description, alpha = 1), linewidth = 2) +
  geom_point(aes(color = Description, alpha = 1), size = 4) +
  scale_y_reverse(breaks = 1:nrow(df))+
  # geom_text(data = df %>% filter(virus == "Coinf"),
  #           aes(label = Description, x = 0.5) , hjust = .85, fontface = "bold", color = "#888888", size = 4)+
  theme(legend.position = "none")+
  theme_bw()



