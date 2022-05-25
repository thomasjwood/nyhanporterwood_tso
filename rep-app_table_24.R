library(plyr)
library(multcomp)
library(survey)
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

t1 <- "https://github.com/thomasjwood/nyhanporterwood_tso/raw/main/t7.rds" %>% 
  url %>% 
  gzcon %>% 
  readRDS

# setting up modeliing list columns ---------------------------------------

# add perceived issue importance

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

# remove pid controls from pid subgroups



# extracting coefficients from lasso models -------------------------------


mdf$l_coefs <- mdf %>% 
  dplyr::select(modgrp, l_coefs) %>% 
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

# load ipw weights

t_ipw <- "https://github.com/thomasjwood/nyhanporterwood_tso/raw/main/t_ipw.rds" %>% 
  url %>% 
  gzcon %>% 
  readRDS


# merge the weights to data

mdf$data %<>% 
  map(
    \(i)
    i %>% 
      left_join(
        t_ipw %>% 
          rename(
            workerid = MID
          ),
        by = "workerid"
      ) %>% 
      filter(
        ipw %>% 
          is.na %>% 
          not
      )
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

# survey weights

mdf$s_m <- mdf %>% 
  select(
    dvs, cond, data, l_coefs
  ) %>% 
  pmap(
    possibly(
      function(dvs, cond, data, l_coefs){
        
        svd <- svydesign(
          ids = ~1, 
          weights = data$ipw, 
          data = data
        )
        
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
          svyglm(design = svd, family = gaussian())
      }, 
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


# emmeans for svyglm objects ----------------------------------------------
mdf$s_emm <- mdf %>% 
  select(
    cond, dvs, data, s_m 
  ) %>% 
  pmap(
    function(cond, dvs, data, s_m){
      
      str_c(
        "trt.vs.ctrl ~ ",
        cond %>% 
          str_replace("\\s\\*\\s", " | ")
      ) %>% 
        as.formula %>% 
        emmeans(
          object = s_m,
          adjust = "none"
        )
    }
  )

mdf$s_conts <- mdf %>% 
  select(
    cond, dvs, s_emm
  ) %>% 
  pmap(
    function(cond, dvs, s_emm){
      
      s_emm %>% 
        extract2("contrasts") %>% 
        tidy %>% 
        mutate(
          cond = cond,
          dvs = dvs
        )
    }
  )


m_4 <- mdf %>% 
  filter(
    
    cond %>% 
      equals("cond_across") &
      dvs %>% 
      is_in(
        c("dv_cc_enrgymix",
          "dv_cc_govdo")
      )
  )



m_4$lc_1 <- m_4$s_m %>% 
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

m_4$lc_2 <- m_4$s_m %>% 
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

m_4$lc_3 <- m_4$s_m %>% 
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


m_4_1 <- m_4 %>% 
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





m_4$s_m %>% 
  map2_dfr(
    m_4$cond,
    \(i, j){
      
      k <- j %>% 
        str_detect(" * ") %>% 
        ifelse(
          j %>% 
            str_sub(15),
          j
        )
      
      tibble(
        conts = i %>% 
          tidy  %>% 
          slice(-1) %>% 
          filter(
            term %>% 
              str_detect("cond_", negate = T) 
          ) %>% 
          filter(
            term %>% 
              str_detect(k, negate = T)
          ) %>% 
          use_series(term)
      )
    }, 
    .id = "mod"
  ) %>% 
  mutate(
    mod = mod %>% as.numeric,
    inc = 1
  ) %>% 
  spread(
    mod, inc
  )


# 

mdf2 <- expand.grid(
  dvs = c("dv_scitrust", "dv_cc_import"),
  wave = "wave 4",
  cond = "cond_across"
) %>% 
  left_join(
    tibble(
      data = mdf %>% 
        filter(
          wave == "wave 4"
        ) %>% 
        slice(1) %>% 
        extract2("data")
    ), 
    by = character() 
  ) %>% 
  tibble


mdf2$l_m <- mdf2 %>% 
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



mdf2$l_coefs <- mdf2$l_m %>% 
  map(
    ~coef(.) %>% 
      tidy %>% 
      slice(-1)
  )


mdf2$s_m <- mdf2 %>% 
  select(
    dvs, cond, data, l_coefs
  ) %>% 
  pmap(
    possibly(
      function(dvs, cond, data, l_coefs){
        
        svd <- svydesign(
          ids = ~1, 
          weights = data$ipw, 
          data = data
        )
        
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
          svyglm(design = svd, family = gaussian())
      },
      NULL
    )
  )

mdf2$lc_1 <- mdf2$s_m %>% 
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

mdf2$lc_2 <- mdf2$s_m %>% 
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

mdf2$lc_3 <- mdf2$s_m %>% 
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




m_5_1 <- mdf2 %>% 
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


mdf2$s_m %>% 
  set_attr(
    "names", 
    mdf2$wave %>% 
      str_replace("wave ", "W")
  ) %>% 
  modelsummary::modelsummary(
    stars = c(
      "*" = .05,
      "**" = .01,
      "***" = .005
    ),
    output = "gt",  
    add_rows = tribble(
      ~contrast,
      "Auxiliary quantities"
    ) %>% 
      bind_rows(
        m_5_1 %>% 
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
    coef_map = mdf2$data[[1]]$cond_7 %>% 
      levels %>% 
      set_attr(
        "names",
        mdf2$data[[1]]$cond_7 %>% 
          levels %>% 
          str_c("cond_across", .)
      ), 
    gof_omit = c(
      "R2 Adj"
    )
  ) %>%
  gt::tab_spanner(
    label = 'Science trust',
    columns = 2
  ) %>% 
  gt::tab_spanner(
    label = 'Issue importance',
    columns = 3
    )