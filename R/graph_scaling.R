library(ggplot2)
library(reshape2)
library(dplyr)

tbl <- read.csv("graph_scaling.csv")

tbl <- tbl %>% group_by(n, t) %>% summarise_all(mean)


plot <- ggplot(tbl, aes(n, t, fill=closedProbability)) + geom_tile()
ggsave("uniform-geometric-sat-scaling-n-t-closedProbability.pdf", width=5, height=3)

plot <- ggplot(tbl, aes(n, t, fill=fourCycles)) + geom_tile()
ggsave("uniform-geometric-sat-scaling-n-t-fourCycles.pdf", width=5, height=3)