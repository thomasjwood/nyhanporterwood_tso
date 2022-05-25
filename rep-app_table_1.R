library(plyr)
library(tidyverse)
library(magrittr)

t1 <- "https://github.com/thomasjwood/nyhanporterwood_tso/raw/main/t1.rds" %>% 
  url %>% 
  gzcon %>% 
  readRDS

t2 <- t1 %>% 
  select(
    cond_7, wave, workerid, education:race, -cont_pk
  ) %>% 
  mutate(
    ideol = ideol %>% 
      mapvalues(
        1:7,
        c("liberal", "moderate", "conservative") %>% 
          rep(c(3, 1,3))
      )
  ) %>% 
  gather(
    cat, val, education:race
  ) %>% 
  mutate(
    val = val %>% 
      factor(
        c("18-34", "35-44", "45-54", "55-64", "65+",
          
          "college", "no college", 
          
          "Female", "Male", "Other",
          
          "conservative", "moderate", "liberal", 
          
          "democrat", "republican", "Non-white", "White")
      )
  ) %>% 
  group_by(
    cond_7, wave, cat, val
  ) %>% 
  tally %>%
  mutate(
    n = n %>% 
      divide_by(
        n %>% sum
      ) %>% 
      multiply_by(100) %>% 
      round
  ) %>% 
  unite(
    "wc", c(wave, cat)
  ) %>% 
  spread(
    cond_7, n
  ) %>% 
  na.omit %>% 
  filter(
    wc %>% 
      str_detect("wave 1")
  )


t2 %>% 
  bind_rows(
    t1 %>% 
      select(
        cond_7, wave, workerid
      ) %>% 
      group_by(
        cond_7, wave
      ) %>% 
      tally %>% 
      spread(
        cond_7, n
      ) %>% 
      slice(1) %>% 
      mutate(
        wc = "",
        val = "n"
      ) %>% 
      select(-wave)
  ) %>% 
  slice(
    c(65, 1:64)
  )