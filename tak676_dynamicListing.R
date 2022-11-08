library(tidyverse)
library(ggplot2)
library(magrittr)
library(stringr)

#import file (crossover subj only)
crossover.Oct_2022_conc <- read.csv('./Data/20221028_SDTM_cohort/PC_TAK-676-1002_Dynamic Listings_25-Oct-2022_crossover.csv',na.string = "n/a")

#select ID and ARM to work with
crossover.Oct_2022_conc.SUBJID.ARM <- crossover.Oct_2022_conc %>%
  select(Unique.Subject.Identifier..USUBJID.,
         Description.of.Planned.Arm..ARM.) %>% 
  distinct()


#import file (all subject)
all.Oct_2022_conc <- read.csv('./Data/20221028_SDTM_cohort/PC_TAK-676-1002_Dynamic Listings_25-Oct-2022_all.csv',na.string = "n/a")

#select ID and ARM to work with
all.Oct_2022_conc.SUBJID.ARM <- all.Oct_2022_conc %>%
  select(Unique.Subject.Identifier..USUBJID.,Description.of.Planned.Arm..ARM.) %>% 
  distinct() %>%
  mutate(SUBJID = str_sub(Unique.Subject.Identifier..USUBJID.,start = 14L), #extract SUBJID info
         DOSE = str_extract(Description.of.Planned.Arm..ARM.,               #extract dose into
                            "([:digit:]\\.[:digit:])"),
         PEMBRO = case_when(
           str_detect(Description.of.Planned.Arm..ARM.,"\\-(.+)\\+") == 1 ~ "N", #assign pembro
           str_detect(Description.of.Planned.Arm..ARM.,"\\+") == 1 ~ "Y",
           
           TRUE ~ "N"
         ))

#add patient tracker info for side-by-side comparison
cohort.SDTM.Oct_2022 <- all.Oct_2022_conc.SUBJID.ARM %>%
  left_join(
    select(filter(Oct_2022_cohort, Visit.Name == "Enrollment (Cycle 1 Day 1)"),
           Subject.Number,
           Pembrolizumab.Visit.Dose,
           Visit.Name), 
    by = c("SUBJID" = "Subject.Number")) %>%
  
  mutate(PEMBRO2 = if_else(is.na(Pembrolizumab.Visit.Dose), "N","Y"))%>%
  rename_with(~paste0(.x,".SDTM"), 1:5)%>%
  rename_with(~paste0(.x,".PtTracker"), 6:8) %>%
  
  mutate(PEMBRO_CONFLICT = if_else(PEMBRO.SDTM == PEMBRO2.PtTracker, "N","Y"))
#export the finding
write.csv(cohort.SDTM.Oct_2022,"./output/20221028_SDTM_cohort/cohort_compare_SDTM_vs_PtTracker_OCT_2022.csv", row.names = FALSE)  

  
