library(plyr)
library(tidyverse)
library(magrittr)
library(broom)

t5 <- "https://github.com/thomasjwood/nyhanporterwood_tso/raw/main/t5.rds" %>% 
  url %>% 
  gzcon %>% 
  readRDS


mac <- "https://github.com/thomasjwood/nyhanporterwood_tso/raw/main/t_mac.rds" %>% 
  url %>% 
  gzcon %>% 
  readRDS

# t tests

t6 <- t5 %>% 
  select(
    MID, any_attrit, starts_with("num_"), pt_sci_conf
  ) %>% 
  gather(
    val, est, -c(MID, any_attrit)
  ) %>% 
  nest_by(val) %>% 
  mutate(
    cond = "overall"
  ) %>% 
  bind_rows(
    t5 %>% 
      select(
        MID, any_attrit, starts_with("num_"), pt_sci_conf
      ) %>% 
      left_join(mac) %>% 
      gather(
        val, est, -c(MID, any_attrit, cond)
      ) %>% 
      nest_by(val, cond)
  )

# i <- t6$data[[1]]

t6$ttests <- t6$data %>% 
  map(
    \(i)
    t.test(
      est ~ any_attrit, 
      data = i, 
      var.equal = F
    )
  )

t7 <- t6 %>% 
  ungroup %>% 
  select(
    val, ttests, cond
  ) %>% 
  pmap_dfr(
    \(val, ttests, cond)
    
    ttests %>% 
      tidy %>% 
      mutate(val, cond)
  ) %>% 
  select(val, everything()) %>% 
  rename(
    attrit_mean =  estimate1, 
    nonatt_mean = estimate2
  ) %>% 
  mutate(
    adj.p.value = p.value %>% 
      p.adjust(method = "BH")
  )



t7p <- t7 %>% 
  modify_at(
    .at = c("p.value", "adj.p.value"),
    ~round(., 4)
  ) %>% 
  arrange(
    cond, adj.p.value
  ) %>% 
  mutate(
    val = val %>% 
      mapvalues(
        
        c("num_nonwhite", "num_ideol", "pt_sci_conf", "num_pt_cc_enrgymix", 
          "num_pid_democrat", "num_pid_republican", "num_college", "num_pt_cc_import", 
          "num_pt_cc_cause", "num_gender"),
        c("Race-Non white",
          "Ideology-Conservative",
          "Science confidence",
          "Policy - support renewable energy",
          "Partisanship - democrat",
          "Partisanship - republican",
          "Demography - college educated",
          "Policy - CC important",
          "Policy - CC has human cause",
          "Demography - Non male")
      )
  ) %>% 
  select(
    cond, val, attrit_mean, nonatt_mean, estimate, statistic, p.value, adj.p.value
  ) %>% 
  rename(
    "v" = val,
    "Mean (attriters)" = attrit_mean,
    "Mean (non attriters)" = nonatt_mean,
    "Difference" = estimate,
    "t" = statistic,
    "p" = "p.value",
    "adjusted p" = adj.p.value
  ) %>% 
  mutate(
    test_type = "t test"
  ) %>% 
  bind_rows(
    t5 %>% 
      select(
        MID, age, any_attrit
      ) %>% 
      mutate(
        cond = "Overall"
      ) %>% 
      bind_rows(
        t5 %>% 
          select(
            MID, age, any_attrit
          ) %>% 
          left_join(mac)
      ) %>% 
      nest_by(cond) %>% 
      pmap_dfr(
        \(cond, data)
        
        table(
          data$age, 
          data$any_attrit
        ) %>% 
          chisq.test %>% 
          tidy %>% 
          mutate(
            cond = cond %>% 
              str_to_title,
            v = "Age"
          )
      ) %>% 
      rename(
        "t" = statistic,
        "p" = "p.value"
      ) %>% 
      select(-method) %>% 
      mutate(
        test_type = "chisq"
      )
  )


# adjustment of both sets of pvalues

t7p %<>% 
  mutate(
    `adjusted p` = p %>% 
      p.adjust(method = "BH")
  )


t7p$v %<>% 
  fct_reorder(t7p$`adjusted p`) %>% 
  fct_rev


t7p$plab <- c("***", "**", "*", "") %>% 
  extract(
    t7p$`adjusted p` %>% 
      findInterval(
        c(-Inf, .005, .01, .05, Inf)
      )
  )


t7p$cond %<>% 
  str_to_title %>% 
  fct_inorder



tab_att <- t6 %>% 
  ungroup %>% 
  select(
    val, cond, ttests
  ) %>% 
  pmap_dfr(
    \(val, cond, ttests)
    ttests %>% 
      tidy %>% 
      mutate(
        val, cond
      )
  ) %>% 
  rename(
    mu_diff = estimate,
    mu_attrited = estimate1,
    mu_nonattrite = estimate2
  ) %>% 
  mutate(
    across(
      mu_attrited:mu_nonattrite,
      \(i){
        
        j <-i %>% 
          round(., 2) %>%
          as.character %>%
          str_replace("0.", ".")
        
        j %>% 
          str_detect(
            "\\d\\."
          ) %>% 
          ifelse(
            j %>% 
              str_pad(
                width = 4, side = "right", pad = "0"
              ),
            j %>% 
              str_pad(
                width = 3, side = "right", pad = "0"
              )
          )
      }
    ),
    diff_lab = mu_diff %>% 
      round(3) %>% 
      as.character %>% 
      str_replace("0.", ".") %>% 
      str_c(
        c("***", "**", "*", "") %>% 
          extract(
            p.value %>% 
              findInterval(
                c(-Inf, .001, .005, .05, Inf)
              )
          )
      ),
    val = val %>% 
      mapvalues(
        
        c("num_nonwhite", "num_ideol", "pt_sci_conf", "num_pt_cc_enrgymix", 
          "num_pid_democrat", "num_pid_republican", "num_college", "num_pt_cc_import", 
          "num_pt_cc_cause", "num_gender"),
        c("Race-Non white",
          "Ideology-Conservative",
          "Science confidence",
          "Policy - support renewable energy",
          "Partisanship - democrat",
          "Partisanship - republican",
          "Demography - college educated",
          "Policy - CC important",
          "Policy - CC has human cause",
          "Demography - Non male")
      ),
    cond = cond %>% 
      str_to_title %>% 
      factor(
        c("Overall",
          "Science Science",
          "Science Opinion",
          "Science Placebo",
          "Science Partisan",
          "Placebo Science", 
          "Partisan Partisan",
          "Placebo Placebo"
        )
      )
  )

ta2 <- tab_att %>% 
  select(
    val, cond, diff_lab, mu_attrited, mu_nonattrite
  ) %>% 
  gather(
    quant, x, -c(1:2)
  ) %>% 
  unite(
    "cq", cond, quant, sep = "_"
  ) %>% 
  mutate(
    cq = cq %>% 
      factor(
        str_c(
          c("Overall",
            "Science Science",
            "Science Opinion",
            "Science Placebo",
            "Science Partisan",
            "Placebo Science", 
            "Partisan Partisan",
            "Placebo Placebo"
          ) %>% 
            rep(each = 3),
          "_",
          c("mu_attrited", "mu_nonattrite", "diff_lab")
        )
      )
  ) %>% 
  filter(
    cq %>% 
      str_detect("diff_")
  ) %>% 
  spread(cq, x)  

t5 %>% 
  select(MID, any_attrit) %>% 
  mutate(
    cond = "Overall"
    ) %>% 
  bind_rows(
    t5 %>% 
      select(MID, any_attrit) %>% 
      left_join(
        mac
      )
    ) %>% 
  group_by(
    cond, any_attrit
  ) %>% 
  tally %>% 
  mutate(perc = n %>% divide_by(n %>% sum) %>% multiply_by(100) %>% round(1)) %>% 
  filter(any_attrit %>% equals("Attrited")) %>% 
  select(-n) %>% 
  spread(cond, perc)

