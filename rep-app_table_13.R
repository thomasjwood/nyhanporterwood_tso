
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


mdf$frms <- mdf$dvs %>% 
  map(
    \(i)
    
    str_c(
      i, " ~ cond_7"
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

mdf$linc <- mdf$r_m %>% 
  map(
    \(i){
      
      c("`cond_7science partisan` - `cond_7science placebo` = 0",
        "`cond_7science opinion` - `cond_7science placebo` = 0",
        "`cond_7science science` - `cond_7partisan partisan` = 0") %>% 
        map_dfr(
          \(j)
          i %>% 
            glht(
              linfct = j
            ) %>% 
            tidy
        )
      
    }
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
    coef_map =  t1$cond_7 %>%
      levels %>% 
      str_to_title %>% 
      set_attr(
        "names",
        c(
          t1$cond_7 %>%
            levels %>%
            str_c("cond_7", .)
        )
      ),
    add_rows = mdf %>% 
      slice(13:14) %>% 
      use_series(linc) %>% 
      bind_rows %>% 
      mutate(
        est = estimate %>% 
          round(4) %>%
          format(scientific = F) %>% 
          as.character %>% 
          str_replace("0.", ".") %>% 
          str_c(
            c("***", "**", "*", "") %>% 
              extract(
                findInterval(
                  adj.p.value,
                  c(-Inf, .005, .01, .05, Inf)
                )
              ),
            "(",
            std.error %>%
              round(3) %>% 
              as.character %>% 
              str_replace("0.", "."),
            ")"
          ),
        mod = 1:2 %>%
          rep(each = 3),
        contrast = contrast %>% 
          str_remove_all("cond_7") %>% 
          str_to_title %>% 
          str_remove_all("`")
      ) %>% 
      select(contrast, mod, est) %>% 
      spread(mod, est) %>% 
      mutate(
        contrast = contrast %>% 
          str_replace_all("Science", "Sci") %>% 
          str_replace_all("Opinion", "Opi") %>% 
          str_replace_all("Placebo", "Pla") %>% 
          str_replace_all("Partisan", "Par") %>% 
          str_replace_all(" - ", "-") %>%
          str_replace_all(" ", ".")
      ) %>% 
      set_attr(
        "position", 13:15
      )
  ) %>% 
  gt::tab_spanner(label = 'Trust scientists', columns = 2) %>% 
  gt::tab_spanner(label = 'CC is important', columns = 3)