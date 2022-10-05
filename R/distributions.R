library(ggplot2)
library(reshape2)
library(dplyr)

tbl <- read.csv("distributions-d.csv")

# tbl <- tbl %>%
#    group_by(n, t) %>%
#    summarise_all(mean)


plot <- ggplot(tbl, aes(closedProbability, color = factor(t))) +
    geom_density() + facet_grid(rows=vars(d))
ggsave("distribution-prob.pdf", width = 5, height = 10)


plot <- ggplot(tbl, aes(fourCycles, color = factor(t))) +
    geom_density() + facet_grid(rows=vars(d))
ggsave("distribution-cycles.pdf", width = 5, height = 10)
