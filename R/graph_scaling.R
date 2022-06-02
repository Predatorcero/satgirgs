library(ggplot2)
library(reshape2)

tbl <- read.csv("../graph_scaling.csv")

ggplot(tbl, aes(t, clustering, group=t)) + geom_boxplot() + xlab('t') + ylab('Closed 4-cycle probability')
