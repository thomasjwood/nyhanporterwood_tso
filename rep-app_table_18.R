library(plyr)
library(tidyverse)
library(magrittr)

t4 <- "https://github.com/thomasjwood/nyhanporterwood_tso/raw/main/t4.rds" %>% 
  url %>% 
  gzcon %>% 
  readRDS

t4 %>% 
  group_by(
    cond, any_attrit
    ) %>% 
  tally %>% 
  spread(
    any_attrit, n
    )

table(
  t4$cond, t4$any_attrit
  ) %>% 
  chisq.test