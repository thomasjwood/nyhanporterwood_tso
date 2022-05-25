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
library(ggtext)


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
  # dvs
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

# remove pid controls from pid subgroups



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



t3 <- mdf %>% 
  filter(
    modgrp %>%
      str_detect("rq5 ") & 
      cond %>%
      str_detect("cond_across ")
  ) %>% 
  filter(
    dvs %>% 
      str_detect("scitrust", negate = T)
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

t3$c2 %<>% 
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


t3$dvs %<>% 
  mapvalues(
    c("dv_cc_cause", "dv_cc_happen", 
      "dv_cc_enrgymix", "dv_cc_govdo"),
    c("atop('    Climate change has    ', 'human cause')",
      "atop('Climate change is', 'happening')",
      "atop('Favor more', 'renewable energy')",
      "atop('Govt should do more to ', 'redress climate change')")
  )

t3$lab <- t3$estimate %>% 
  round(2) %>% 
  equals(0) %>% 
  ifelse(
    t3$estimate %>% 
      round(3) %>% 
      as.character %>% 
      str_replace("0\\.", "\\."),
    t3$estimate %>% 
      round(2) %>% 
      as.character %>% 
      str_replace("0\\.", "\\.")
  ) %>% 
  str_c(
    t3$p.value %>% 
      gtools::stars.pval() %>% 
      str_remove_all("\\.") %>% 
      str_trim
  )

p5 <- t3 %>%
  filter(
    dvs %>% 
      str_detect("Climate change")
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
      color = pid 
    ),
    size = .25,
    position = position_dodge(width = .2)
  ) +
  geom_line(
    aes(
      wave %>% 
        factor %>% 
        as.numeric, 
      estimate,
      color = pid
    ),
    size = .25,
    position = position_dodge(width = .2),
  ) +
  geom_point(
    aes(
      wave,
      estimate,
      color = pid,
      size = p.value %>% 
        is_less_than(.05)
    ),
    shape = 21,
    fill = "white",
    position = position_dodge(width = .2)
  ) +
  geom_text(
    aes(
      wave, estimate, label = lab,
      color = pid, 
    ),
    size = 1.85,
    position = position_dodge(width = .2),
    data = t3 %>%
      filter(
        dvs %>% 
          str_detect("Climate change") &
          p.value < .05
      ) %>% 
      mutate(
        c2 = c2 %>% 
          str_to_title %>% 
          factor(
            c2 %>% 
              levels %>% 
              str_to_title
          )
      )
  ) +
  facet_grid(
    dvs ~ c2,
    scales = "free_y", 
    switch = "y",
    label = label_parsed
  ) +
  labs(
    x = "",
    y = "",
    title = "Treatment effects on scientific understanding",
    subtitle = "Effects reported among <span style = 'color:#2fabe1'>**Democratic  identifiers/leaners**</span> and <span style = 'color:#d2181f'>**Republican  identifiers/leaners**</span> (pre-treatment)"
  ) +
  scale_x_discrete(
    # breaks = NULL 
    labels = 2:4
  ) +
  scale_y_continuous(position = "right") +
  scale_color_manual(
    values = c("#2fabe1",
               "#d2181f")
  ) +
  scale_size_manual(
    values = c(1.5, 8)
  ) +
  theme(
    strip.text.x = element_text(size = 8),
    plot.title.position = "plot", 
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    legend.position = "none",
    plot.subtitle = element_markdown(),
  )


p6 <- t3 %>%
  filter(
    dvs %>% 
      str_detect("Climate change", negate = T)
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
      color = pid 
    ),
    size = .25,
    position = position_dodge(width = .2)
  ) +
  geom_line(
    aes(
      wave %>% 
        factor %>% 
        as.numeric, 
      estimate,
      color = pid
    ),
    size = .25,
    position = position_dodge(width = .2),
  ) +
  geom_point(
    aes(
      wave,
      estimate,
      color = pid,
      size = p.value %>% 
        is_less_than(.05)
    ),
    shape = 21,
    fill = "white",
    position = position_dodge(width = .2)
  ) +
  geom_text(
    aes(
      wave, estimate, label = lab,
      color = pid, 
    ),
    size = 1.85,
    family = "Inter",
    position = position_dodge(width = .2),
    data = t3 %>%
      filter(
        dvs %>% 
          str_detect("Climate change", negate = T) &
          p.value < .05
      ) %>% 
      mutate(
        c2 = c2 %>% 
          str_to_title %>% 
          factor(
            c2 %>% 
              levels %>% 
              str_to_title
          )
      )
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
    subtitle = "Effects reported among <span style = 'color:#2fabe1'>**Democratic identifiers/leaners**</span> and <span style = 'color:#d2181f'>**Republican identifiers/leaners**</span> (pre-treatment)",
    labels = label_parsed
  ) +
  scale_x_discrete(
    labels = 2:4
    # breaks = NULL 
  ) +
  scale_y_continuous(position = "right") +
  scale_color_manual(
    values = c("#2fabe1",
               "#d2181f")
  ) +
  scale_size_manual(
    values = c(1.5, 8)
  ) +
  theme(
    plot.title.position = "plot",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    strip.text.x = element_text(size = 8),
    plot.subtitle = element_markdown(),
    legend.position = "none",
  )

p5 / p6  +
  plot_layout(heights = c(1.05, 1))

