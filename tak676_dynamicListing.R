library(tidyverse)
library(ggplot2)
library(magrittr)
library(stringr)

#import file (crossover subj only)
#crossover.Oct_2022_conc <- read.csv('./Data/20221028_SDTM_cohort/PC_TAK-676-1002_Dynamic Listings_25-Oct-2022_crossover.csv',na.string = "n/a")

#select ID and ARM to work with
#crossover.Oct_2022_conc.SUBJID.ARM <- crossover.Oct_2022_conc %>%
#  select(Unique.Subject.Identifier..USUBJID.,
#         Description.of.Planned.Arm..ARM.) %>% 
#  distinct()


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


  
