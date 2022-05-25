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
library(gt)

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
  # str_c("cond_", c(3, 7)),
  # str_c("cond_", c(3, 7))
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




# removing wave 4 from cond_within ----------------------------------------

# mdf %<>% 
#   filter(
#     (
#       cond %>% 
#         str_detect("cond_within") &
#         wave %>% 
#         equals("wave 4")
#     ) %>%
#       not
#   )



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


# keep only cond_acrpss

mdf %<>%
  filter(
    # wave %>%
    #   is_in(
    #     "wave " %>%
    #       str_c(2:3)
    #   ) &
      cond %>%
      str_detect("cond_across")
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
        # 
        # dvs <- mdf$dvs[[62]]
        # data <- mdf$data[[62]]
        
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
        
        # dvs <- mdf$dvs[[20]]
        # cond <- mdf$cond[[20]]
        # data <- mdf$data[[20]]
        # l_coefs <- mdf$l_coefs[[20]]
        
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
      
      # cond <- mdf$cond[[25]]
      # dvs <- mdf$dvs[[25]]
      # data <- mdf$data[[25]]
      # r_m <- mdf$r_m[[25]]
      
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
      
      # cond <- mdf$cond[[25]]
      # dvs <- mdf$dvs[[25]]
      # emm <- mdf$emm[[25]]
      
      
      emm %>% 
        extract2("contrasts") %>% 
        tidy %>% 
        mutate(
          cond = cond,
          dvs = dvs
        )
    }
  )



t4 <- mdf %>% 
  filter(
    modgrp %>% 
      equals("rq4 policy models")
  )

t4$lc_1 <- t4$r_m %>% 
  map(
    function(i){
      
      i %>% 
        glht(
          linfct = c(
            "`cond_acrossscience partisan` - `cond_acrossscience placebo` = 0"
          )
        )
    }
  )

t4$lc_2 <- t4$r_m %>% 
  map(
    function(i){
      
      i %>% 
        glht(
          linfct = c(
            "`cond_acrossscience opinion` - `cond_acrossscience placebo` = 0"
          )
        )
    }
  )

t4$lc_3 <- t4$r_m %>% 
  map(
    function(i){
      
      i %>% 
        glht(
          linfct = c(
            "`cond_acrossscience science` - `cond_acrosspartisan partisan` = 0"
          )
        )
    }
  )


t5 <- t4 %>% 
  select(
    wave, dvs, lc_1:lc_3
  ) %>% 
  pivot_longer(
    starts_with("lc_")
  ) %>% 
  pmap_dfr(
    function(wave, dvs, name, value)
      
      value %>% 
      summary(
        test = adjusted("none")
      ) %>% 
      tidy %>% 
      mutate(
        type = name,
        dvs = dvs,
        wave = wave
      )
  ) %>% 
  mutate(
    lab = estimate %>% 
      round(3) %>% 
      equals(0) %>% 
      ifelse(
        estimate %>% 
          round(5) %>%
          as.character %>% 
          str_replace("0\\.", "\\."),
        estimate %>% 
          round(3) %>%
          as.character %>% 
          str_replace("0\\.", "\\.")
      ) %>% 
      str_c(
        c("***", 
          "**",
          "*",
          "") %>% 
          extract(
            p.value %>% 
              findInterval(
                c(-Inf, .005, .01, .05, Inf)
              )
          ),
        "(",
        std.error %>% 
          round(3) %>% 
          as.character %>% 
          str_sub(2),
        ")"
      ),
    contrast = contrast %>% 
      str_remove_all("cond_across")
  ) %>% 
  unite(
    "dw", 
    dvs, wave
  ) %>% 
  mutate(
    dw = dw %>% 
      fct_inorder
  ) %>% 
  select(
    contrast, dw, lab
  ) %>% 
  spread(dw, lab)




t4$r_m %>% 
  set_attr(
    "names", t4$wave %>% str_replace("wave ", "W")
  ) %>% 
  modelsummary::modelsummary(
    stars = c(
      "*" = .05,
      "**" = .01,
      "***" = .005
      # Significance threshold: p < 0.05 (two-sided; we will also report if p<.01 and p<.005)
    ),
    output = "gt",
    add_rows = tribble(
      ~contrast,
      "Auxiliary quantities"
    ) %>% 
      bind_rows(
        t5 %>% 
          mutate(
            contrast = contrast %>% str_to_title
          )
      ) %>% 
      modify(
        function(i)
          i %>% 
          is.na %>% 
          ifelse(
            "", i
          )
      ) %>% 
      modify(
        function(i)
          i %>% 
          str_remove_all("`")
      ),
    coef_map = t4$data[[1]]$cond_7 %>%
      levels %>%
      set_attr(
        "names",
        t4$data[[1]]$cond_7 %>%
          levels %>%
          str_c("cond_across", .)
      )
  ) %>% 
  gt::tab_spanner(label = 'Favor renewable energy', columns = 2:4) %>% 
  gt::tab_spanner(label = 'Favor Govt action', columns = 5:7)