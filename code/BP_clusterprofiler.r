#### R Code for Gene Enrichment Analysis ###

library(clusterProfiler)
library(ggplot2)

### Read input file ####
de <- read.csv("BP_list", header = TRUE)
ego <- enrichGO(de$Entrezid, OrgDb = "org.Hs.eg.db", ont="BP", readable = TRUE, keyType = "SYMBOL") ### Use ont="BP" for biological processes, and ont="MF" for molecular functions
ego2 <- simplify(ego, cutoff=0.8, by="p.adjust", select_fun=min, measure = "Wang")
dotplot(ego2, x = "GeneRatio", color = "p.adjust", showCategory = 20, font.size = 10, label_format = 50)
ggsave("Output.png", units="in", width=7, height=6, dpi=600)
write.csv(ego2, file = "BP_complete_Output")
