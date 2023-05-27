###### Code for generating Barplots #######

library(ggplot2)
df = read.csv("input.csv", header = TRUE)
p<-ggplot(df, aes(y = reorder(CellType, +Fraction), x = Fraction)) + geom_col(aes(fill = Fraction), width = 0.5)
pp <- p + ggtitle("LumB Rep") + xlab("Log Observe/Expected Fraction(%)") + ylab("CellType") + theme(plot.title = element_text(hjust = 0.5))

pp + theme(
     plot.title = element_text(color="black", size=14, face="bold"),
     axis.text.x = element_text(color="black", size=10, face="bold"),
     axis.text.y = element_text(color="black", size=10, face="bold"),
     axis.title.x = element_text(color="black", size=14, face="bold"),
     axis.title.y = element_text(color="black", size=14, face="bold"))

ggsave("Image.png", units="in", width=6, height=4, dpi=600)
dev.off()

