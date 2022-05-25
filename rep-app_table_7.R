library(plyr)
library(multcomp)
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

l_coef <- {
  set.seed(1234)
  
  "dv_cc_import ~ pt_cc_cause  + pt_cc_enrgymix + import_cc + education + age + gender + pid + ideol + cont_pk + race + cont_polint + ind_scitrust + region"  %>% 
    as.formula %>% 
    cv.glmnet(
      alpha = 1,
      data = t1,
      family = "gaussian"
    ) %>% 
    coef %>% 
    tidy %>% 
    slice(-1) %>% 
    use_series(row)
}

rm <- lm_robust(
  str_c(
    "dv_cc_import ~ ",
    "cond_7 + ",
    l_coef %>%       
      str_remove("pidrepublican") %>% 
      str_replace("piddemocrat", "pid") %>% 
      stringi::stri_remove_empty() %>% 
      str_c(collapse = " + ")
  ) %>% 
    as.formula, 
  data = t1
)

rm %>% 
  list %>% 
  set_attr(which = "names", " ") %>% 
  modelsummary::modelsummary(
    stars = c(
      "*" = .05,
      "**" = .01,
      "***" = .005
    ),
    output = "gt",
    coef_map = t1$cond_7 %>% 
      levels %>% 
      str_to_title %>% 
      set_attr(
        "names",
        t1$cond_7 %>% 
          levels %>% 
          str_c("cond_7", .)
      ),
    add_rows = tribble(
      ~contrast,
      "Auxiliary quantities"
    ) %>% 
      bind_rows(
        c("`cond_7science partisan` - `cond_7science placebo` = 0",
          "`cond_7science opinion` - `cond_7science placebo` = 0",
          "`cond_7science science` - `cond_7partisan partisan` = 0") %>% 
          map_dfr(
            \(i)
            rm %>% 
              glht(
                linfct = i
              ) %>% 
              tidy
          ) %>% 
          mutate(
            contrast = contrast %>%
              str_remove_all("cond_7") %>% 
              str_to_title,
            estimate = estimate %>% 
              round(3) %>% 
              as.character %>% 
              str_sub(2) %>% 
              str_c(
                "(",
                std.error %>% 
                  round(3) %>% 
                  as.character %>% 
                  str_sub(2),
                ")"
              )
          ) %>% 
          select(contrast, estimate)
      ) %>% 
      set_attr(which = "position", 13:16)
  )