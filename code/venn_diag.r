#### R code ###
#### Venn Diagram ###

library("ggVennDiagram")
library("ggplot2")library("ggVennDiagram")
library("ggplot2")
data<-read.csv("input.csv", header = TRUE)   ### Comma separated file of genes in each subtype

p <- ggVennDiagram(data[1:4], label = "count")                  ## Count will shown only common numbers and not percentage
p + scale_fill_distiller(palette = "Reds", direction = 1)

ggsave("Image.png", units="in", width=6, height=4, dpi=600, bg = 'white')
dev.off()
