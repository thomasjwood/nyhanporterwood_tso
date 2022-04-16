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


# setting up modeliing list columns ---------------------------------------


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
    str_c("cond_", c("within", "across"))
  ),
  list(
    c("dv_cc_happen", "dv_cc_cause"),
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


mdf %<>% 
  filter(
    (
      cond %>% 
        str_detect("cond_within") &
        wave %>% 
        equals("wave 4")
    ) %>%
      not
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


# robust regression models ------------------------------------------------

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

# plotting main effects ---------------------------------------------------

t2 <- mdf %>% 
  filter(
    (
      modgrp %>%
        str_detect("main models") |
        modgrp %>%
        str_detect("rq4 ")
    ) &
      cond %>% 
      equals("cond_across")
  ) %>% 
  select(
    modgrp, wave, dvs, conts
  ) %>% 
  pmap_dfr(
    function(modgrp, wave, dvs, conts)
      
      conts %>% 
      mutate(
        modgrp = modgrp,
        wave = wave, 
        dvs = dvs
      )
  ) %>% 
  mutate(
    c2  = contrast %>% 
      str_sub(
        ,
        contrast %>% 
          str_locate(" - ") %>% 
          extract(, 1) %>% 
          subtract(1)
      )
  )

t2$c2 %<>% 
  plyr::mapvalues(
    c("science science",
      "science opinion", 
      "science placebo",
      "science partisan",
      "placebo science",
      "partisan partisan"
    ),
    c("Science~~W2%->%Science~~W3",
      "Science~~W2%->%Opinion~~W3", 
      "Science~~W2%->%Placebo~~W3",
      "Science~~W2%->%Partisan~~W3",
      "Placebo~~W2%->%Science~~W3",
      "Partisan~~W2%->%Partisan~~W3")
    
  ) %>% 
  factor(
    c("Science~~W2%->%Science~~W3",
      "Science~~W2%->%Opinion~~W3", 
      "Science~~W2%->%Placebo~~W3",
      "Science~~W2%->%Partisan~~W3",
      "Placebo~~W2%->%Science~~W3",
      "Partisan~~W2%->%Partisan~~W3")
  )

t2$dvs %<>% 
  mapvalues(
    c("dv_cc_cause", "dv_cc_happen", 
      "dv_cc_enrgymix", "dv_cc_govdo"),
    c("atop('    Climate change has    ', 'human cause')",
      "atop('Climate change is', 'happening')",
      "atop('Favor renewable', 'energy')",
      "atop('Favor', '     government action     ')")
  )

t2$lab <- t2$estimate %>% 
  round(2) %>% 
  equals(0) %>% 
  ifelse(
    t2$estimate %>% 
      round(3) %>% 
      as.character %>% 
      str_replace("0\\.", "\\."),
    t2$estimate %>% 
      round(2) %>% 
      as.character %>% 
      str_replace("0\\.", "\\.")
  ) %>% 
  str_c(
    t2$p.value %>% 
      gtools::stars.pval() %>% 
      str_remove_all("\\.") %>% 
      str_trim
  )


p1 <- t2 %>%
  filter(
    modgrp %>% 
      str_detect(" main ")
  ) %>% 
  mutate(
    c2 = c2 %>% 
      str_to_title %>% 
      factor(
        c2 %>% 
          levels %>% 
          str_to_title
      )
  ) %>% 
  ggplot() +
  geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  geom_linerange(
    aes(
      wave, 
      estimate,
      ymin = estimate %>% 
        subtract(
          std.error %>% 
            multiply_by(1.96)
        ),
      ymax = estimate %>% 
        add(
          std.error %>% 
            multiply_by(1.96)
        ),
      
    ),
    size = .25
  ) +
  geom_line(
    aes(
      wave %>% 
        factor %>% 
        as.numeric, 
      estimate
    ),
    size = .25
  ) +
  geom_point(
    aes(
      wave,
      estimate
    ),
    shape = 21,
    fill = "white",
    size = 7.5
  ) +
  geom_text(
    aes(
      wave, estimate, label = lab
    ),
    size = 1.85
  ) +
  facet_grid(
    dvs ~ c2,
    scales = "free_y", 
    switch = "y",
    # labeller(c2 = label_parsed)
    labeller = label_parsed
  ) +
  labs(
    x = "",
    y = "",
    title = "Treatment effects on scientific understanding",
  ) +
  scale_x_discrete(
    labels = 2:4
  ) +
  scale_y_continuous(position = "right") +
  theme(
    plot.title.position = "plot",
    strip.text.y.left = element_text(angle = 0, lineheight = .5),
    strip.placement = "outside",
    legend.position = "none", 
    plot.margin = margin(.15, .15, -.2, .15, unit = "cm")
  )

p2 <-  t2 %>%
  filter(
    modgrp %>% 
      str_detect(" main ", negate = T)
  ) %>% 
  mutate(
    c2 = c2 %>% 
      str_to_title %>% 
      factor(
        c2 %>% 
          levels %>% 
          str_to_title
      )
  ) %>% 
  ggplot() +
  geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  geom_linerange(
    aes(
      wave, 
      estimate,
      ymin = estimate %>% 
        subtract(
          std.error %>% 
            multiply_by(1.96)
        ),
      ymax = estimate %>% 
        add(
          std.error %>% 
            multiply_by(1.96)
        ),
      
    ),
    size = .25
  ) +
  geom_line(
    aes(
      wave %>% 
        factor %>% 
        as.numeric, 
      estimate
    ),
    size = .25
  ) +
  geom_point(
    aes(
      wave,
      estimate
    ),
    shape = 21,
    fill = "white",
    size = 7.5
  ) +
  geom_text(
    aes(
      wave, estimate, label = lab
    ),
    size = 1.85,
  ) +
  scale_y_continuous(position = "right") +
  scale_x_discrete(
    labels = 2:4
  ) +
  facet_grid(
    dvs ~ c2,
    scales = "free_y", 
    switch = "y",
    labeller = label_parsed
    
  ) +
  labs(
    x = "Wave",
    y = "",
    title = "Treatment effects on policy attitudes",
  ) +
  theme(
    plot.title.position = "plot",
    strip.text.y.left= element_text(angle = 0), 
    strip.placement = "outside",
    legend.position = "none", 
    plot.margin = margin(-2.15, .15, .15, .15, unit = "cm")
  )


p1 / p2 

