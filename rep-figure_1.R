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
  list(
    c("cond_within", "cond_across")
  ) %>%
    rep(times = 3),
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


d3 <- mdf %>% 
  filter(
    cond %>% 
      equals("cond_within")
  ) %>% 
  select(
    modgrp:dvs, emm
  ) %>% 
  pmap_dfr(
    \(modgrp, wave, cond, dvs, emm)
    
    emm$contrasts %>% 
      tidy %>% 
      mutate(
        modgrp, wave, cond, dvs
      )
  ) %>% 
  mutate(
    cond = contrast %>% 
      str_remove(" - placebo") %>% 
      str_to_title,
    lab = estimate %>% 
      round(3) %>% 
      str_replace("0.", ".") %>% 
      str_c(
        c("***", "**", "*", "") %>% 
          extract(
            p.value %>%
              findInterval(
                c(-Inf, .005, .01, .05, Inf)
              )
          )
      ),
    dvs = case_when(
      dvs %>% 
        equals("dv_cc_cause") ~ "Climate change has human cause",
      dvs %>% 
        equals("dv_cc_happen") ~ "Climate change is happening",
      dvs %>% 
        equals("dv_cc_enrgymix") ~ "Favor renewable energy",
      dvs %>% 
        equals("dv_cc_govdo") ~ "Favor government redress climate change"
    ) %>% 
      fct_inorder,
    modgrp = modgrp %>% 
      str_detect("h1 ") %>% 
      ifelse(
        "Scientific understanding",
        "Policy beliefs"
      ),
    lo = estimate %>% 
      subtract(
        std.error %>% 
          multiply_by(1.96)
      ),
    hi = estimate %>% 
      add(
        std.error %>% 
          multiply_by(1.96)
      ),
    wave = wave %>% 
      str_to_title
  )

p_1 <- d3 %>% 
  filter(
    dvs %>% 
      is_in(
        dvs %>% 
          levels %>% 
          extract(1:2)
      )
  ) %>% 
  ggplot() +
  geom_vline(
    xintercept = 0,
    size = .25, 
    linetype = "dashed" 
  ) +
  geom_linerange(
    aes(
      xmin = lo, xmax = hi,
      y = cond
    )
  ) +
  geom_point(
    aes(
      estimate, cond
    ),
    size = 13,
    shape = 21,
    fill = "white"
  ) +
  geom_text(
    aes(
      estimate, cond, label = lab
    ),
    size =  2.5
  ) +
  facet_grid(
    wave ~ dvs, scales = "free", space = "free_y"
  ) +
  labs(
    x = "",
    y = "",
    title = "Scientific understanding"
  ) +
  theme(
    strip.text.y = element_blank(),
    plot.title.position = "plot",
    strip.text.y.left = element_text(angle = 0, lineheight = .5),
    strip.placement = "outside",
    legend.position = "none", 
  )

p_2 <- d3 %>% 
  filter(
    dvs %>% 
      is_in(
        dvs %>% 
          levels %>% 
          extract(1:2)
      ) %>% 
      not
  ) %>% 
  ggplot() +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    size = .25, 
  ) +
  geom_linerange(
    aes(
      xmin = lo, xmax = hi,
      y = cond
    )
  ) +
  geom_point(
    aes(
      estimate, cond
    ),
    size = 13,
    shape = 21,
    fill = "white"
  ) +
  geom_text(
    aes(
      estimate, cond, label = lab
    ),
    size =  2.5,
  ) +
  facet_grid(
    wave ~ dvs, scales = "free", space = "free_y"
  ) +
  labs(
    x = "",
    y = "",
    title = "Policy attitudes"
  ) +
  scale_y_discrete(
    breaks = NULL
  ) +
  theme(
    plot.title.position = "plot",
    strip.text.y.right = element_text(angle = 0, lineheight = .5),
    legend.position = "none", 
  )

p_1 | p_2 
