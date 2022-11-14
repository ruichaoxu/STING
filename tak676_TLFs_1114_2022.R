library(tidyverse)
library(ggplot2)
library(magrittr)
library(readxl)
library(stringr)

#1.0 Import----

PP.OCT27.2022.source <- read_xlsx("./Data/20221114_TLFs/TAK676-1002 PK parameters_10272022.xlsx")
SDTM.OCT25.2022.source <- read_excel('./Data/20221114_TLFs/TAK-676-1002_Dynamic Listings_25-Oct-2022.xlsm', sheet = "PC")







#2.1 Trim PP.source----

unique(PP.OCT27.2022.source$Parameter)
pp.parameter.selector <- c("Tmax", 
                           "Cmax", 
                           "AUCINF_obs", 
                           "AUClast", 
                           "HL_Lambda_z",
                           "Cl_obs",
                           "Vss_obs")
PP.OCT27.2022.source.trimmed <- PP.OCT27.2022.source %>%
  select(DOSE,PEMBRO,SUBJID, VISIT, Parameter, Units, Estimate)%>%
  filter(Parameter %in% pp.parameter.selector) %>%
  filter((DOSE<=2.5&PEMBRO == "Y") | (DOSE<= 3.5&PEMBRO == "N")) %>%
  mutate(VISIT = str_to_upper(VISIT))







#2.2 Trim SDTM----

SDTM.PC.colnames <- colnames(SDTM.OCT25.2022.source)
SDTM.col.selector <-  SDTM.PC.colnames[c(1,2,3,20,24,25)] #based on communication with SQS

SDTM.OCT25.2022.source.trimmed <- SDTM.OCT25.2022.source %>%
  select(SDTM.col.selector) %>%
  filter(`Specimen Material Type (PCSPEC)` != "URINE") #urine data only contains volume rather than concentration

SDTM.OCT25.2022.source.trimmed.mod <- SDTM.OCT25.2022.source.trimmed %>%
  mutate(SUBJID = str_sub(`Unique Subject Identifier (USUBJID)`, start =14L)) %>% #extract subject ID from USUBJID
  filter(`Visit Number (VISITNUM)` %in% c("1.01","1.08")) %>% #filter out irrelevant visits
  drop_na() %>%
  filter(`Epoch (EPOCH)`!= "SCREENING") %>%
  group_by(SUBJID, `Visit Name (VISIT)`) %>%
  slice(1)







#3.0 Merge PP&SDTM----


PP.OCT27.2022.merged <- PP.OCT27.2022.source.trimmed %>%
  left_join(SDTM.OCT25.2022.source.trimmed.mod,by = c("SUBJID", "VISIT" = "Visit Name (VISIT)"))






#4.0 Finalize----

#PP name

#Replace name based on reference table
#create a reference table, containing old and new names 
#creat vector, assign new name to vecter value, assign old name to vector row names
#use the vector in str_replace_all fucntion
pp.name.reference <- tibble(PHOENIX = pp.parameter.selector, 
                            PPTESTCD  = c("TMAX", "CMAX", "AUCIFO","AUC024","LAMZHL","CLO","VSSO"),
                            PPTEST = c("Time of CMAX",
                                       "Max Concentration",
                                       "AUC Infinity Obs",
                                       "AUC 0-24hr",
                                       "Half-Life",
                                       "Clearance",
                                       "Volume of Distribution"))

lookup.1 <- pp.name.reference$PPTESTCD
names(lookup.1) <- pp.name.reference$PHOENIX

lookup.2 <-pp.name.reference$PPTEST
names(lookup.2) <- pp.name.reference$PHOENIX
#_______________________________________________________________________




PP.OCT27.2022.final <- PP.OCT27.2022.merged %>%
  mutate(USUBJID = `Unique Subject Identifier (USUBJID)`,
         PPGRPID = VISIT,
         
         PPTESTCD = str_replace_all(Parameter,lookup.1),
         PPTEST = str_replace_all(Parameter,lookup.2),
         PPCAT = "TAK676",
         PPSCAT = "TAK676 PK",
         
         PPORRES = Estimate,
         PPORRESU = Units,
         PPSTRESC = Estimate,
         PPSTRESN = Estimate,
         PPSTRESU = Units,
         VISIT = VISIT,
         VISITNUM = `Visit Number (VISITNUM)`,
         EPOCH = `Epoch (EPOCH)`,
         PCDTC = `Date/Time of Specimen Collection (PCDTC)`,
         PPSPEC = `Specimen Material Type (PCSPEC)`,
         
         .keep = "unused") %>%
  select(-PEMBRO,-SUBJID,-DOSE) %>% 
  relocate(VISIT, .before  =VISITNUM)

write.csv(PP.OCT27.2022.final, "./output/20221114_TLFs/PP.csv", row.names = F)
