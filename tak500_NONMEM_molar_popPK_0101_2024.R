library(tidyverse)
library(ggplot2)
library(magrittr)
library(stringr)
library(mrgsolve)
library(linpk)

MW_tak500 <- 146035
MW_tak676 <- 754.48
MW_PAYLOAD <- 1779.57


#0.SPEC----
# dose - nmol/kg (dose reduction incorporated)
# DV - nM


#1.0 CONC----

# tak_500_DV <- tak500_data_1110_2023%>%
#   mutate(PCORRES = as.numeric(PCORRES)) %>%
#   drop_na(PCORRES) %>%
#   mutate(DV = as.character(PCORRES)) %>%
#   
#   mutate(TIME = if_else(
#     PCTPTNUM == -1, 0, PCTPTNUM+1
#   )) %>%
#   
#   mutate(DVID = case_when(
#     PCTEST == "Total Antibody" ~ "1",
#     PCTEST == "Conjugated TAK-676" ~ "2",
#     PCTEST == "Deconjugated TAK-676" ~ "3"
#   )) %>%
#   
#   mutate(EVID = 0) %>%
#   mutate(
#     # CMT = 
#     #   case_when(
#     #     DVID == "1" ~ 1, #CMT1 - totl.mAb in CENTRAL
#     #     DVID == "2" ~ 3, #CMT3 - bald mAb in CENTRAL
#     #     DVID == "3" ~ 99),
#     CMT = 1,
#     AMT = ".",
#     DUR = ".",
#     ID = ".",
#     MDV = 0) %>%
#   
#   select(SUBJID, ID, AMT, DUR,CMT, TIME, DV, DVID,MDV, EVID )

#1.1 CONC_nM----

tak_500_DV_mol <- tak500_data_1110_2023%>%
  mutate(PCORRES = as.numeric(PCORRES)) %>%
  drop_na(PCORRES) %>%
  
  mutate(DVID = case_when(
    PCTEST == "Total Antibody" ~ "1",
    PCTEST == "Conjugated TAK-676" ~ "2",
    PCTEST == "Deconjugated TAK-676" ~ "3"
  )) %>%
  
  filter(DVID != "3") %>% #free DAZO not applicable for modeling now
  
  mutate(DV = as.character(if_else(DVID=="1", PCORRES/145000*1000, PCORRES/MW_tak676*1000))) %>% #transform ng/mL conc. into nM for 2 species
  
  mutate(
    TAD = if_else(PCTPTNUM == -1, 0, PCTPTNUM+1),
    TIME = round(VISITNUM-1)*24*21+PCTPTNUM+1) %>%
  
  
  mutate(EVID = 0) %>%
  mutate(
    # CMT = 
    #   case_when(
    #     DVID == "1" ~ 1, #CMT1 - totl.mAb in CENTRAL
    #     DVID == "2" ~ 3, #CMT3 - bald mAb in CENTRAL
    #     DVID == "3" ~ 99),
    CMT = 1,
    AMT = ".",
    DUR = ".",
    ID = ".",
    MDV = 0,
    TYPE = DVID) %>%
  arrange(SUBJID, TIME) %>%
  
  select(SUBJID, ID, TYPE,AMT, DUR,CMT, TIME, TAD, DV, DVID,MDV, EVID )


#2.0 ID ------

tak_500_DV_ID <- tak_500_DV_mol %>%
  select(SUBJID) %>%
  group_by(SUBJID) %>%
  arrange(SUBJID) %>%
  slice(1)%>%
  ungroup(SUBJID) %>%
  mutate(ID = as.character(row_number(SUBJID)))




#3.0 DOSE ------

##3.1 DOSE old ------
# tak_500_AMT <- tak_500_AUG_dosing %>%
#   
#   mutate(TIME = (Cycle-1)*7*3*24) %>% #hours
#   
#   mutate(AMT = as.character(DOSE.reduced_ug.kg))%>% #ug/mg
#   mutate(DUR = "1") %>%
#   
#   mutate(ID=".",
#          DV = ".",
#          MDV = 1,
#          EVID = 1,
#          CMT = 1,
#          DVID = "0",
#   ) %>%
#   
#   select(SUBJID, ID, AMT, DUR,CMT, TIME, DV, DVID,MDV, EVID)


##3.2 DOSE new molar per kg ------

#total mAb
tak_500_AMT_mol_tMab <- tak_500_AUG_dosing %>%
  
  mutate(TIME = (Cycle-1)*7*3*24) %>% #hours
  
  mutate(AMT = as.character(DOSE.reduced_ug.kg/145000*1000))%>% #nmol/kg
  mutate(DUR = "1") %>%
  
  mutate(ID=".",
         DV = ".",
         MDV = 1,
         EVID = 1,
         CMT = 1,
         DVID = "0",
         TYPE = "1"
  ) %>%
  
  select(SUBJID, ID, TYPE,AMT, DUR,CMT, TIME, DV, DVID,MDV, EVID)

#conjugated DAZO
tak_500_AMT_mol_acDAZO <- tak_500_AUG_dosing %>%
  
  mutate(TIME = (Cycle-1)*7*3*24) %>% #hours
  
  mutate(AMT = as.character(4*DOSE.reduced_ug.kg/145000*1000))%>% #nmol/kg, assuming DAR = 4
  mutate(DUR = "1") %>%
  
  mutate(ID=".",
         DV = ".",
         MDV = 1,
         EVID = 1,
         CMT = 1,
         DVID = "0",
         TYPE = "2"
  ) %>%
  
  select(SUBJID, ID, TYPE,AMT, DUR,CMT, TIME, DV, DVID,MDV, EVID)


#combining 2 species

tak_500_AMT_mol <- tak_500_AMT_mol_tMab %>%
  bind_rows(tak_500_AMT_mol_acDAZO) %>%
  arrange (SUBJID, TIME, TYPE)






##3.3 DOSE molar ------

#total mAb
tak_500_AMT_mol_tMab_2 <- tak_500_AUG_dosing %>%
  
  mutate(TIME = (Cycle-1)*7*3*24) %>% #hours
  
  mutate(AMT = as.character(Dose.actual_ug/145000*1000))%>% #nmol
  mutate(DUR = "1") %>%
  
  mutate(ID=".",
         DV = ".",
         MDV = 1,
         EVID = 1,
         CMT = 1,
         DVID = "0",
         TYPE = "1"
  ) %>%
  
  select(SUBJID, ID, TYPE,AMT, DUR,CMT, TIME, DV, DVID,MDV, EVID)

#conjugated DAZO
tak_500_AMT_mol_acDAZO_2 <- tak_500_AUG_dosing %>%
  
  mutate(TIME = (Cycle-1)*7*3*24) %>% #hours
  
  mutate(AMT = as.character(4*Dose.actual_ug/145000*1000))%>% #nmol/kg, assuming DAR = 4
  mutate(DUR = "1") %>%
  
  mutate(ID=".",
         DV = ".",
         MDV = 1,
         EVID = 1,
         CMT = 1,
         DVID = "0",
         TYPE = "2"
  ) %>%
  
  select(SUBJID, ID, TYPE,AMT, DUR,CMT, TIME, DV, DVID,MDV, EVID)


#combining 2 species

tak_500_AMT_mol_2 <- tak_500_AMT_mol_tMab_2 %>%
  bind_rows(tak_500_AMT_mol_acDAZO_2) %>%
  mutate(TAD = 0) %>%
  arrange (SUBJID, TIME, TYPE)



#4.0 Merge----

tak_500_NONMEM_Jan2024_2 <- tak_500_DV_mol %>%
  bind_rows(tak_500_AMT_mol_2) %>%
  
  select(-ID)%>%
  left_join(tak_500_DV_ID, by = "SUBJID") %>%
  select(-DVID)%>%
  
  mutate(CMT = if_else(TYPE =="1", 1,3))%>% #tot mAb in cmt1, acDAZO in cmt3
  
  arrange(as.numeric(ID), TIME, AMT, TYPE) %>%
  
  
  mutate(LDV = if_else(DV ==".", DV, as.character(log(as.numeric(DV))))) %>%
  mutate(C = ".")%>% 
  
  filter(!is.na(TIME))%>%
  
  
  select(C,SUBJID, ID, TYPE, TAD, TIME, CMT, DV, AMT, DUR, EVID, MDV,LDV)


write.csv(tak_500_NONMEM_Jan2024_2, '../TAK500/output/20240103_NONMEM/tak_500_NONMEM_Jan2024_v2.csv', row.names = FALSE)
