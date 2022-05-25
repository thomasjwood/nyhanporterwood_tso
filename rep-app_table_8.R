# setting up modeliing list columns ---------------------------------------

library(plyr)
library(factoextra)
library(psych)
library(glmnet)
library(glmnetUtils)
library(broom)
library(estimatr)
library(tidyverse)
library(magrittr)

t1 <- "https://github.com/thomasjwood/nyhanporterwood_tso/raw/main/t1.rds" %>% 
  url %>% 
  gzcon %>% 
  readRDS


# add perceived issue importance

l1 <- list(
  # dvs
  list(
    "dv_cc_happen",
    "dv_cc_cause",
    "dv_cc_govdo",
    "dv_cc_enrgymix", 
    "dv_scitrust",
    "dv_cc_import"
  ),
  # waves
  list(
    str_c("wave ", 2:4),
    str_c("wave ", 2:4),
    str_c("wave ", 2:4),
    str_c("wave ", 2:4),
    str_c("wave ", 4),
    str_c("wave ", 4)
  )
)


mdf <- l1 %>% 
  pmap_dfr(
    function(i, j)
      
      expand_grid(
        dvs = i,
        wave = j %>% unlist
      ) %>% 
      arrange(
        dvs, wave
      )
  )  %>% 
  left_join(
    t1 %>% 
      group_by(wave) %>% 
      nest
  )


mdf$l_m <- mdf %>% 
  select(
    dvs, data
  ) %>% 
  pmap(
    possibly(
      function(
    dvs, data
      ){
        
        set.seed(1234)
        # 
        # dvs <- mdf$dvs[[62]]
        # data <- mdf$data[[62]]
        
        dvs %>% 
          str_c(
            " ~ pt_cc_cause  + pt_cc_enrgymix + import_cc + education + age + gender + ideol + cont_pk + race + cont_polint + ind_scitrust + region"
          ) %>% 
          as.formula %>% 
          cv.glmnet(
            alpha = 1,
            data = data,
            family = "gaussian"
          )
      }, 
    NULL
    )
  )

mdf$frms <- mdf$dvs %>% 
  map2(
    mdf$l_m %>% 
      map(
        ~coef(.) %>%
          tidy %>%
          slice(-1) %>% 
          use_series(row)
      ),
    \(i, j)
    
    str_c(
      i,
      " ~ cond_7 * pid + ",
      j %>% 
        unlist %>% 
        str_c(collapse = " + ")
    ) %>% 
      as.formula
  )

mdf$r_m <- mdf$frms %>% 
  map2(
    mdf$data,
    \(i, j)
    
    lm_robust(
      i,
      se_type = "HC2", 
      data =  j
    )
  )


mdf %>% 
  slice(1:12) %>%
  use_series(r_m) %>% 
  set_attr(
    "names", 
    mdf %>% 
      slice(1:12) %>% 
      use_series(wave) %>% 
      str_replace("wave ", "W")
  ) %>% 
  modelsummary::modelsummary(
    stars = c(
      "*" = .05,
      "**" = .01,
      "***" = .005
      # Significance threshold: p < 0.05 (two-sided; we will also report if p<.01 and p<.005)
    ),
    output = "gt",
    coef_map = c(
      t1$cond_7 %>%
        levels,
      t1$cond_7 %>%
        levels %>% 
        str_c(" * Republican")
    ) %>% 
      str_to_title %>% 
      set_attr(
        "names",
        c(
          t1$cond_7 %>%
            levels %>%
            str_c("cond_7", .),
          t1$cond_7 %>%
            levels %>%
            str_c("cond_7", . ,":pidrepublican")
        )
      )
  ) %>% 
  gt::tab_spanner(label = 'CC is happening', columns = 2:4) %>% 
  gt::tab_spanner(label = 'CC has human cause', columns = 5:7) %>% 
  gt::tab_spanner(label = 'Favor govt action', columns = 8:10) %>% 
  gt::tab_spanner(label = 'Favor renewable energy', columns = 11:13)