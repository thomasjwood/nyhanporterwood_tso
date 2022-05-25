library(plyr)
library(multcomp)
library(factoextra)
library(psych)
library(glmnet)
library(glmnetUtils)
library(broom)
library(estimatr)
library(emmeans)
library(patchwork)
library(tidyverse)
library(magrittr)


t1 <- "https://github.com/thomasjwood/nyhanporterwood_tso/raw/main/t1.rds" %>% 
  url %>% 
  gzcon %>% 
  readRDS


l1 <- list(
  # effect names
  list(
    "h1 & h2 main models",
    "rq4 policy models",
    "sci affect"
  ),
  # relevant waves 
  list(
    str_c("wave ", 2:4),
    str_c("wave ", 2:4),
    str_c("wave ", 4)
  ),
  # conditional indicator
  # list(
  list(
    c("cond_within", "cond_across")
  ) %>%
    rep(times = 3),
  # dvs
  list(
    c("dv_cc_happen", "dv_cc_cause", "dv_cc_import"),
    c("dv_cc_govdo", "dv_cc_enrgymix"),
    "dv_scitrust"
  )
)

mdf <- l1 %>% 
  pmap_dfr(
    function(i, j, k, l)
      
      expand_grid(
        modgrp = i,
        wave = j,
        cond = k,
        dvs = l %>% 
          unlist
      ) %>% 
      arrange(
        modgrp, dvs
      )
  ) %>% 
  bind_rows(
    # partisan difference
    l1 %>% 
      pmap_dfr(
        function(i, j, k, l)
          
          expand_grid(
            modgrp = "rq5 partisan differences",
            wave = j,
            cond = k %>% 
              str_c(" * pid"),
            dvs = l %>% 
              unlist
          ) %>% 
          arrange(
            modgrp, dvs
          )
      ),
    l1 %>% 
      pmap_dfr(
        function(i, j, k, l)
          
          expand_grid(
            modgrp = "rq6 anthro prebelief",
            wave = j,
            cond = k %>% 
              str_c(" * ind_reject_anthro"),
            dvs = l %>% 
              unlist
          ) %>% 
          arrange(
            modgrp, dvs
          )
      )
  ) %>% 
  left_join(
    t1 %>% 
      group_by(wave) %>% 
      nest
  )



# removing waves 2 and 3 from dv_cc_import --------------------------------

mdf %<>%  
  filter(
    (
      dvs %>% 
        equals("dv_cc_import") &
        wave %>% 
        is_in(
          c("wave 2",
            "wave 3")
        )
    ) %>%
      not
  )


# keep only cond_within

mdf %<>% 
  filter(
    wave %>% 
      is_in(
        "wave " %>% 
          str_c(2:3)
      ) &
      cond %>% 
      str_detect("cond_within")
  )


# estimating lasso models -------------------------------------------------

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
            " ~ pt_cc_cause  + pt_cc_enrgymix + import_cc + education + age + gender + pid + ideol + cont_pk + race + cont_polint + ind_scitrust + region"
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
      slice(-1)
  )



# extracting coefficients from lasso models -------------------------------


mdf$l_coefs <- mdf %>% 
  select(modgrp, l_coefs) %>% 
  pmap(
    function(modgrp, l_coefs)
      
      if(
        modgrp %>% 
        str_detect("partisan")
      ) {
        l_coefs %>%
          filter(
            row %>% 
              str_detect("pid") %>% 
              not
          )
      } else {
        l_coefs
      } 
  )


# regression models ------------------------------------------------

mdf$r_m <- mdf %>% 
  select(
    dvs, cond, data, l_coefs
  ) %>% 
  pmap(
    possibly(
      function(dvs, cond, data, l_coefs)
        
        dvs %>% 
        str_c(
          " ~ ",
          cond,
          " + ",
          l_coefs$row %>% 
            unlist %>%
            str_replace("piddemocrat", "pid") %>% 
            str_replace("raceWhite", "race") %>% 
            str_remove("pidrepublican") %>% 
            str_remove("raceNon-white") %>% 
            stringi::stri_remove_empty() %>% 
            str_c(collapse = " + ")
        ) %>% 
        as.formula %>% 
        lm_robust(
          se_type = "HC2", data =  data
        ), 
      NULL
    )
  )




# emmeans objects from robust models --------------------------------------


mdf$emm <- mdf %>% 
  select(
    cond, dvs, data, r_m
  ) %>% 
  pmap(
    function(cond, dvs, data, r_m){
      
      str_c(
        "trt.vs.ctrl ~ ",
        cond %>% 
          str_replace("\\s\\*\\s", " | ")
      ) %>% 
        as.formula %>% 
        emmeans(
          object = r_m,
          adjust = "none"
        )
    }
  )

mdf$conts <- mdf %>% 
  select(
    cond, dvs, emm
  ) %>% 
  pmap(
    function(cond, dvs, emm){
      
      emm %>% 
        extract2("contrasts") %>% 
        tidy %>% 
        mutate(
          cond = cond,
          dvs = dvs
        )
    }
  )

tl_1_1 <- mdf %>% 
  filter(
    dvs %>% 
      is_in(
        c("dv_cc_cause", "dv_cc_happen")
      ) 
    # &
    # cond %>% 
    #   str_detect("cond_within\\s\\*")
  ) %>% 
  arrange(dvs)

tl_1_1$r_m %>% 
  set_attr(
    "names", 
    tl_1_1$wave %>% 
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
    coef_map = tl_1_1$data[[1]]$cond_within %>% 
      levels %>%
      str_to_title %>% 
      c(
        "Partisan x Republican",
        "Science x Republican",
        "Opinion x Republican",
        "Partisan x Reject Anthro CC", 
        "Science x Reject Anthro CC", 
        "Opinion x Reject Anthro CC"
      ) %>% 
      set_attr(
        "names",
        tl_1_1$data[[1]]$cond_within %>% 
          levels %>% 
          str_c("cond_within", .) %>% 
          c("cond_withinpartisan:pidrepublican",
            "cond_withinscience:pidrepublican",
            "cond_withinopinion:pidrepublican",
            "cond_withinpartisan:ind_reject_anthroReject anthropogenic cc", 
            "cond_withinscience:ind_reject_anthroReject anthropogenic cc", 
            "cond_withinopinion:ind_reject_anthroReject anthropogenic cc")
      )
  ) %>% 
  gt::tab_spanner(label = 'Climate change has human cause', 2:7) %>% 
  gt::tab_spanner(label = 'Climate change is happening', columns = 8:13) 
