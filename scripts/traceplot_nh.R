args <- commandArgs(trailingOnly=T)

library(ggplot2)
library(reshape2)

# args <- c("traceplotfile.txt",output file name)

# read in trace file

trace.df <- read.table(args[1],header=T)
trace.df$runnum <- as.factor(trace.df$runnum)
trace.df_long <- melt(trace.df,id.vars=c("itnum","runnum"))

# create trace plot

pdf(args[2])
ggplot(trace.df_long,aes(x=itnum,y=value,col=runnum,linetype=variable)) +
  geom_line(alpha=0.5) +
  ggtitle(args[2])
dev.off()
