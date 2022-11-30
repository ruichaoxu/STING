library(tidyverse)
library(ggplot2)
library(magrittr)
library(stringr)

nov_2022_urine_conc <- read.csv('./Data/20221123_urine/TAK-676-1002_QPS_Urine_2022-11-22.csv',na.string = "n/a")
nov_2022_urine_vol <- Oct_2022_conc %>%
  filter(PCSPEC == "URINE"&PCTESTCD == "VOLUME") %>% #find urine volume
  select(USUBJID, PCORRES) %>%
  mutate(Volume = as.numeric(PCORRES), .keep = "unused")

nov_2022_urine <- nov_2022_urine_conc %>%
  left_join(nov_2022_urine_vol, by = "USUBJID") %>%
  mutate(UrineConc = as.numeric(PCORRES)) %>% 
  filter(!is.na(UrineConc)&!is.na(Volume)) %>%
  mutate(UrineAMT_mg =UrineConc*Volume/1000000) %>%
  left_join(Oct_2022_cohort_std.SDTM, by = "SUBJID") %>%
  mutate(PctRecoveryUrine = UrineAMT_mg/as.numeric(DOSE)*100)

write.csv(nov_2022_urine, './output/20221123_urine/TAK-676-1002_Urine_recov_11222022.csv', row.names = FALSE)

hist(filter(nov_2022_urine,Volume>50)$PctRecoveryUrin, breaks = 10)
mean(nov_2022_urine$PctRecoveryUrine) #9.68%
sd(nov_2022_urine$PctRecoveryUrine)/mean(nov_2022_urine$PctRecoveryUrine) #0.7363885
median(nov_2022_urine$PctRecoveryUrine) #8.60%
range(nov_2022_urine$PctRecoveryUrine) #0.006% - 26.1% 
  