library(tidyverse)

cor_data <- read_csv("data/cor_result.csv")
cor_data %>%
  filter(!is.na(cor))# %>%
  ggplot(aes(x=cor)) +
  geom_histogram() +
#  scale_y_log10() +
  xlab("PCC") +
  theme_linedraw() +
  geom_vline(xintercept = cor_data %>% arrange(-cor) %>% filter(row_number() <= 4) %>% pull(cor)) +
  theme(axis.text.x = element_text(size =30), axis.text.y = element_text(size=30), aspect.ratio = 1/3)

cor_data %>%
  filter(!is.na(cor)) %>%
  ggplot(aes(y=cor)) +
  geom_boxplot() +
  #  scale_y_log10() +
  ylab("PCC")

single_cell <- read_tsv("data/rna_single_cell_type.tsv")
colnames(single_cell) <- c("Gene", "Gene_name", "Cell_type", "nTPM")

top_5_genes_in_mouse <- c("ARRB1", "RNF13", "CPNE9", "TNNC1")
cell_type_list <- c("Rod photoreceptor cells", "Bipolar cells", "Muller glia cells", "Horizontal cells", "
Cone photoreceptor cells")
extracted <- single_cell %>%
  filter(Gene_name %in% top_5_genes_in_mouse) %>%
  filter(Cell_type %in% cell_type_list) %>%
  full_join(tibble(Cell_type=c("Cone photorecepter cells"), Gene_name=top_5_genes_in_mouse, nTPM=0))


heatmap_celltype <- extracted %>%
  group_by(Gene_name) %>% mutate(zscore=scale(nTPM)) %>% ungroup() %>% arrange(Cell_type) %>%
  ggplot(aes(x=Gene_name, y=factor(Cell_type, levels = rev(unique(Cell_type))), fill = zscore)) + geom_tile() +
  geom_text(aes(label = nTPM)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  xlab("Gene name") +
  ylab("Cell type") +
  theme_bw() +
  theme(
    axis.text.x = element_text(face = "italic"),  # x軸ラベルをイタリックに設定
    axis.text = element_text(color = "black")
  )


  ggsave(plot = heatmap_celltype, filename="heatmap.svg", height = 6, width = 5)
