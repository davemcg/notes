# not necessarily runnable code
# more for me to find what i did in the past 
# so I don't have to re-google a bazillion things

#################################################################################
# heatmap 1
# non-default things:
## 1. many categorical columns programatically defined
## 2. as well as their colors
## 3. column order is done myself
## 4. remove legend colors
library(tidyverse)
library(ComplexHeatmap)
### https://bmcophthalmol.biomedcentral.com/articles/10.1186/s12886-021-01830-9#MOESM4

shin_seo_table <- readxl::read_xlsx('2886_2021_1830_MOESM4_ESM.xlsx', skip = 2)

unique_genes <- shin_seo_table$`Nearby/containing Gene` %>% unique() %>% str_split(., "\\/|;") %>% unlist() %>% unique()
### grab go terms for genes in the table
unique_genes_go_terms <-  AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys=c(unique_genes), columns = c('SYMBOL','GENENAME','GO'), keytype = "SYMBOL")
### now get the names of each go term
go_term_names <- list()
for (i in unique_genes_go_terms$GO %>% unique()){
  if (!is.na(i)){
    go_term_names[[i]] <- GO.db::GOTERM[[i]]@Term
  }
}
# align gene names and go terms, filter down to select terms
annotations <- unique_genes_go_terms %>% 
  left_join(go_term_names %>% unlist() %>% enframe(value = 'Description'), by =c("GO" = "name")) %>% 
  filter(grepl("pigmentation|complement|mitoch|lysosom|autophagy|oxidative|angiog", Description)) %>% 
  ungroup() %>% 
  group_by(SYMBOL, Description) %>% 
  summarise(GOterms = paste0(unique(GO), collapse =', ')) %>% 
  mutate(simplified_description = case_when(grepl("mito", Description)~ 'mitochondria',
                                            grepl("lyso", Description)~ 'lysosome',
                                            grepl("angiog", Description)~ 'angiogenesis',
                                            grepl("comple", Description)~ 'complement',
                                            grepl("oxidati", Description)~ 'oxidative stress',
                                            grepl("pigmen", Description)~ 'pigmentation'))
### create anno table for HeatmapAnnotation
go_table <- shin_seo_table %>% 
  rowwise() %>% 
  mutate(genes = (str_split(`Nearby/containing Gene`, "\\/|;"))) %>% unnest(genes) %>%  
  left_join(annotations %>% select(SYMBOL, simplified_description) %>% unique(),by = c("genes" = "SYMBOL")) %>% 
  mutate(id = paste0(`SNP ID`, ' - ', genes, ' - ', simplified_description)) %>% 
  filter(!is.na(simplified_description)) %>% 
  select(`SNP ID`, simplified_description) %>% 
  unique() %>% 
  mutate(count = 1) %>% 
  pivot_wider(values_from = count, names_from = simplified_description, values_fill = 0) %>% 
  column_to_rownames('SNP ID')

### value table for the heatmap itself
log10_table <- shin_seo_table %>% 
  filter(`SNP ID` %in% row.names(go_table)) %>% 
  mutate(id = paste0(`SNP ID`, ' - ', gsub("\\/|;", ", ", `Nearby/containing Gene`))) %>% 
  column_to_rownames('id') %>% 
  mutate(across(contains("log10"),as.numeric))
log10_table[is.na(log10_table)] <- 0
snps_with_sig_qvalue <- log10_table %>% 
  select(`SNP ID`, contains('log10')) %>% pivot_longer(-`SNP ID`) %>% 
  mutate(qvalue = value**-10) %>% filter(qvalue < 0.05)

#### sanity check
table(row.names(go_table) == log10_table$`SNP ID`)

log10_table$rmeans <- log10_table %>% select(contains("log10")) %>% rowMeans()
### custom row order
row_order <- cbind(go_table, log10_table) %>% 
  arrange(complement, angiogenesis, mitochondria, lysosome, pigmentation, `oxidative stress`, rmeans)
### anno function for white/black
col_fun = circlize::colorRamp2(c(0,1), c("white","black"))
col_list <- list()
### build list structure for colors
for (i in c('complement', 'angiogenesis', 'mitochondria', 'lysosome', 'pigmentation', 'oxidative stress')){
  col_list[[i]] <- col_fun
}
### customize row order and remove legend
ha_col <- HeatmapAnnotation(df =go_table[row.names(row_order),
                                         c('complement', 'angiogenesis', 'mitochondria', 'lysosome', 'pigmentation', 'oxidative stress')],
                            col = col_list,
                            show_legend = FALSE)
### enhance colors a bit (some outliers)
hm_col <-  circlize::colorRamp2(c(-100, 0, 100), c("blue", "white", "red"))
Heatmap(log10_table[row.names(row_order),] %>% select(contains("log10")) %>% as.matrix() %>% t(),
        top_annotation = ha_col, 
        name = 'log10 qvalue', 
        col = hm_col, 
        cluster_columns = FALSE,
        column_names_max_height = unit(20,'cm'))
##############################################################################
