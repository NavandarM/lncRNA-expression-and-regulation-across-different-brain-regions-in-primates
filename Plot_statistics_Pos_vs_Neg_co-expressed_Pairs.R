library("ggplot2")

Statistics <- read.delim("DF_for_positive_Negative_correlations.txt")

ggplot(Statistics, aes(x=Organism, y= Fractions, fill=Type)) + geom_bar(width=0.5, stat = "identity") + coord_polar("y", start=0) + theme_bw() + scale_fill_manual(values = c("#E69F00", "#009E73"))