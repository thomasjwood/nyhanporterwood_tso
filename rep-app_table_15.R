library(plyr)
library(multcomp)
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


t1 %<>% 
  select(
    -cont_pk
  ) %>% 
  left_join(
    t1 %>% 
      filter(
        wave == "wave 1"
      ) %>% 
      select(workerid, cont_pk) %>% 
      mutate(
        cont_pk = cont_pk %>% 
          sjmisc::dicho() %>% 
          mapvalues(
            c(0:1),
            c("low",
              "high")
          )
      )
  )


# setting up modeliing list columns ---------------------------------------

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
        
        dvs %>%
          str_c(
            " ~ pt_cc_cause  + pt_cc_enrgymix + import_cc + education + age + gender + ideol +  race + cont_polint + ind_scitrust + region"
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


mdf$l_coefs <- mdf$l_m %>% 
  map(
    ~coef(.) %>% 
      tidy %>% 
      slice(-1) %>% 
      use_series(row)
  )



mdf$frms <- mdf$dvs %>% 
  map2(
    mdf$l_coefs,
    \(i, j)
    
    str_c(
      i, 
      " ~ cond_7 * cont_pk + ",
      j %>% 
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
  slice(13:14) %>%
  use_series(r_m) %>% 
  set_attr(
    "names", 
    mdf %>% 
      slice(13:14) %>% 
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
        str_c(" * High Pol_Knowl")
    ) %>% 
      str_to_title %>%
      str_replace("_k", ".K") %>% 
      set_attr(
        "names",
        c(
          t1$cond_7 %>%
            levels %>%
            str_c("cond_7", .),
          t1$cond_7 %>%
            levels %>%
            str_c("cond_7", . ,":cont_pkhigh")
        )
      )
    
  ) %>% 
  gt::tab_spanner(label = 'Trust scientists', columns = 2) %>% 
  gt::tab_spanner(label = 'CC is important', columns = 3)