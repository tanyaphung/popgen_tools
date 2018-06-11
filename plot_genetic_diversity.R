library(ggplot2)
library(gridExtra)
library(grid)

inputs = commandArgs(trailingOnly=T)

# Check if there are 2 arguments given. If not, return a help message.
if (length(inputs)!=2) {
  stop("USAGE: Rscript plot_sfs.R /path/to/the/pi_file/ /path/to/the/output.png_file", call.=FALSE)
}

diversity = read.table(inputs[1], header = T)

png(inputs[2], width = 11, height = 8, res = 300, units = "in")

ggplot(data=diversity, aes(x=start, y=pi_per_site)) + 
  geom_point() +
  theme(legend.title=element_blank(), text = element_text(size=14), plot.title = element_text(hjust = 0.5, face="bold")) + 
  labs(title="Genetic diversity", y=expression(paste(pi)), x="Coordinates")

dev.off()