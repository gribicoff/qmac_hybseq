args <- commandArgs(trailingOnly=T)

# args: plink .mds file, output filename stem

library(ggplot2)

m <- read.table(args[1],header=T)
spp <- rep(NA, length(m[,1]))
spp[grep("mac", m[,1])] <- "macrocarpa"
spp[grep("alb", m[,1])] <- "alba"
spp[grep("bic", m[,1])] <- "bicolor"
spp[grep("ste", m[,1])] <- "stellata"
spp[grep("mue", m[,1])] <- "muehlenbergii"

pdf(paste0(args[2],"_macpops_mds.pdf"))
b <- ggplot(m, aes(C1, C2, col = spp)) + geom_point(shape=1,size=2,fill=NA)
b <- b + scale_colour_manual(values = c("red", "blue", "brown", "green", "yellow"))
b <- b + theme_light()
b + xlab("MDS1") + ylab("MDS2") + ggtitle(paste0("MDS Plot - "), args) + theme(plot.title = element_text(hjust = 0.5))
dev.off()
