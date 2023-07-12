# Checking RL2 results.
#
# Jull 11

library(pomp)
library(spatPomp)
library(tidyverse)

results <- readRDS("model3/run_level_2/h3_revision.rds")
# results3 <- readRDS("model3/run_level_3/h3_revision.rds")

s1_traces <- melt(results$search1$traces) %>%
  rename(L1 = `.L1`) %>%
  mutate(iteration = as.numeric(iteration))

s2_traces <- melt(results$search2$traces) %>%
  rename(L1 = `.L1`) %>%
  # filter(as.numeric(iteration) != 0) %>%
  mutate(iteration = as.numeric(iteration) + 100)

s3_traces <- melt(results$search3$traces) %>%
  rename(L1 = `.L1`) %>%
  # filter(as.numeric(iteration) != 0) %>%
  mutate(iteration = as.numeric(iteration) + 200)

unit_specific <- c("betaB", "foi_add", "aHur", "hHur")

shared <- c(paste0(c(
  "mu_B", "XthetaA", "thetaI", "lambdaR", "r", "std_W",
  "epsilon", "k"
), "1"), 'loglik')

s1_topn <- order(-results$search1$logLiks$logLik)[1:9]
s2_topn <- order(-results$search2$logLiks$logLik)[1:6]
s3_top10 <- order(-results$search3$logLiks$logLik)[1:10]

s2_ancestor <- function(id) {
  s1_topn[((id - 1) %/% 4) + 1]
}

s3_ancestor2 <- function(id) {
  s2_topn[((id - 1) %/% 6) + 1]
}

s3_ancestor1 <- function(id) {
  s2_ancestor(s3_ancestor2(id))
}

s2_traces %>%
  mutate(ancest_1 = s2_ancestor(L1)) -> s2_traces

s3_traces %>%
  mutate(ancest_2 = s3_ancestor2(L1),
         ancest_1 = s3_ancestor1(L1)) -> s3_traces

merged_results <- bind_rows(
  filter(s1_traces, L1 %in% unique(s3_traces$ancest_1)) %>%
    mutate(L2 = L1),
  filter(s2_traces, L1 %in% unique(s3_traces$ancest_2)) %>%
    mutate(L2 = ancest_1) %>%
    select(-ancest_1) %>%
    mutate(L1 = L1 + 72),
  s3_traces %>%
    mutate(L2 = ancest_1) %>%
    select(-ancest_1, -ancest_2) %>%
    mutate(L1 = L1 + 72 + 36)
)

# filter(s1_traces, variable %in% shared) %>%
merged_results %>%
  mutate(top_10 = L1 %in% (s3_top10 + 72 + 36)) %>%
  filter(variable %in% shared) %>%
  filter(iteration > 5) %>%
  ggplot(aes(x = iteration, y = value)) +
  geom_line(aes(group = L1, color = top_10)) +
  facet_wrap(~variable, scales = 'free') +
  theme_bw()






