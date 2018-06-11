library(ggplot2)
library(gridExtra)
library(grid)

inputs = commandArgs(trailingOnly=T)

# Check if there are 2 arguments given. If not, return a help message.
if (length(inputs)!=2) {
  stop("USAGE: Rscript plot_sfs.R /path/to/the/sfs_file/ /path/to/the/output.png_file", call.=FALSE)
}

sfs = read.table(inputs[1], header = T)

png(inputs[2], width = 11, height = 8, res = 300, units = "in")
ggplot(data=sfs, aes(x=bin, y=count)) + 
  geom_bar(stat = "identity", position=position_dodge(), colour="black", fill=alpha("grey60", 0.7)) +
  theme(legend.title=element_blank(), text = element_text(size=14), plot.title = element_text(hjust = 0.5, face="bold")) + 
  labs(title="Site frequency spectrum", x="Count", y="Number of variants") +
  scale_x_continuous(breaks = sfs$bin)
dev.off()

