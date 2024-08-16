#aim: Convert popPK dataset NCA dataset suitable for PKanalix

library(tidyverse)

#load data----
data.tak500.NONMEM <- read.csv('../../Monolix/TAK500/Apr2024_popPK/data/tak_500_NONMEM_Jun2024_withBW.csv')




#reconstruct new ID ----
#previous ID was built for popPK, not directly applicable for NCA
data.tak500.NONMEMtoNCA <- data.tak500.NONMEM %>%
  mutate(cycle = floor(as.numeric(TIME)/(24*7*3))+1) %>% #compute treatment cycle from TIME variable
  filter(BQLflag != "1") %>% #No BQL handling for NCA, this was carried over from popPK workflow
  group_by(ID,cycle)%>%
  mutate(ID.cycle = paste(ID,cycle, sep = "."),
         ID.cycle = as.numeric(ID.cycle))%>% #construct ID.cycle variable for subsequent grouping and numbering the sequence
  ungroup()%>%
  group_by(ID.cycle)%>%
  mutate(NCA.ID = cur_group_id()) #Assign new ID
  




#output NCA dataset----
data.tak500.NONMEMtoNCA.tot.mAb <- data.tak500.NONMEMtoNCA %>%
  filter(TYPE == "1")

data.tak500.NONMEMtoNCA.acDAZO <- data.tak500.NONMEMtoNCA %>%
  filter(TYPE == "2")

write.csv(data.tak500.NONMEMtoNCA.tot.mAb, '../../Monolix/TAK500/Apr2024_popPK/data/tak_500_NONMEM_Jun2024_withBW_NCAdataset.totmAb.csv', row.names = F)

write.csv(data.tak500.NONMEMtoNCA.acDAZO, '../../Monolix/TAK500/Apr2024_popPK/data/tak_500_NONMEM_Jun2024_withBW_NCAdataset.acDazo.csv', row.names = F)
