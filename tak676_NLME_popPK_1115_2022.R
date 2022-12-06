library(tidyverse)
library(ggplot2)
library(magrittr)
library(stringr)
library(mrgsolve)
library(datapasta)
library(linpk)

#_____________________----- 
#0. Import data-----
#data was from phoenix dataset
dataAll_1115.2022 <- read.csv("./Data/20221115_NLME_popPK/popPK.conc.OCT.2022.csv")




#_____________________----- 
#1. Data treatment-----


##1.1, rename columns-----
glimpse(dataAll_1115.2022)

dataAll_1115.2022_1.1 <- dataAll_1115.2022 %>%
  group_by(DOSE,SUBJID) %>%
  mutate(
    DV = CONC, 
    AMT = DOSE, 
    OCC = if_else(VISIT == "CYCLE 1 DAY 1",1,2),
    PEMBRO = if_else (PEMBRO == "Y",1,0),
    ID = cur_group_id(),
    .keep = "unused")
write.csv(dataAll_1115.2022_1.1,"./output/20221115_NLME_popPK/xyplot_tak676_all_OCT2022.csv",row.names = FALSE)



##1.2, extract dose info-----
dataAll_1115.2022_doseRow <- dataAll_1115.2022_1.1 %>%
  group_by(SUBJID,OCC) %>%
  slice(1) %>%
  mutate( 
    
    DV = ".",
    AMT = as.character(AMT), #making sure same column type for merging
    TIME = 1,    #force dosing time at 1hr approximating IV infusion
    .keep = "unused"
  )

unique(dataAll_1115.2022_doseRow$TIME)


###1.2.1. Dose as infusion-----
dataAll_1115.2022_doseRow_inf <- dataAll_1115.2022_doseRow %>%
  mutate(TIME = 0,
         DUR = as.character(1))



##1.3. popPK ready worksheet----

dataAll_1115.2022_1.3 <- dataAll_1115.2022_1.1 %>%
  mutate(
    
    DV = as.character(DV),
    AMT = "."
  ) %>%
  union(dataAll_1115.2022_doseRow)%>%
  arrange(ID,OCC,TIME,desc(AMT)) %>%  #Making sure data properly ordered
  mutate(ID2 = paste0(ID,".",OCC))

write.csv(dataAll_1115.2022_1.3,"./output/20221115_NLME_popPK/tak676_all_Oct2022.csv",row.names = FALSE)

###1.3.1.  Infusion dataset-----
dataAll_1115.2022_1.3.1 <- dataAll_1115.2022_1.1 %>%
  mutate(
    DUR = ".", .after = AMT,
    DV = as.character(DV),
    AMT = "."
  ) %>%
  union(dataAll_1115.2022_doseRow_inf)%>%
  arrange(ID,OCC,TIME,desc(AMT)) %>%  #Making sure data properly ordered
  mutate(ID2 = paste0(ID,".",OCC), .after = ID)

write.csv(dataAll_1115.2022_1.3.1,"./output/20221115_NLME_popPK/tak676_all_Oct2022_inf.csv",row.names = FALSE)






###1.3.2.  Infusion dataset de-occasion-----
#Calculate C1D8 time as natural time


dataAll_1115.2022_1.3.2<- dataAll_1115.2022_1.3.1 %>%
  mutate(TIME = TIME + (OCC-1)*24*7)

write.csv(dataAll_1115.2022_1.3.2,"./output/20221115_NLME_popPK/tak676_all_Oct2022_inf_deOCC.csv",row.names = FALSE)











###1.3.3.  Infusion dataset meanConc-----
#Calculate mean concentration among subject

#using earlier dataset (DV as numeric)
dataAll_1115.2022_1.1_meanCON <- dataAll_1115.2022_1.1 %>%
  ungroup()%>%
  group_by(SUBJID,TIME)%>%
  summarise(meanCONC = mean(DV))%>%
  left_join(distinct(select(dataAll_1115.2022_1.1,-TIME,-DV,-OCC)))%>% #extract out dataset for left join (TIME,DV and OCC are not useful info so dropped)
  mutate(DUR = ".")

dataAll_1115.2022_1.1.3<- dataAll_1115.2022_1.1_meanCON  %>%
  mutate(
    
    DV = as.character(meanCONC),
    AMT = ".",
    .keep = "unused"
  ) %>%
  union(select(ungroup(dataAll_1115.2022_doseRow_inf),-OCC))%>%
  arrange(ID,TIME,desc(AMT)) 

write.csv(dataAll_1115.2022_1.1.3,"./output/20221115_NLME_popPK/tak676_all_July2022_inf_meanConc.csv",row.names = FALSE)






#3. POST hoc-----

#model id: deOCC 2cmp pro inf
asianSubj.Aug2022 <- c("07001-003","58001-004","58004-006","58006-003")

##3.0 Data import-----
CLposthoc.Nov2022 <- read.csv("./Data/20221115_NLME_popPK/deOCC _2cmp_pro_inf_postHoc_11152022.csv")

CLposthoc_clean.Nov2022 <- CLposthoc.Nov2022 %>%
  select(-time) %>%
  mutate(V = 1000*V,
         V2 = 1000*V2,
         Cl = 1000*Cl,
         Cl2 = 1000*Cl2)%>% #get the unit right
  left_join(distinct(select(dataAll_1115.2022_1.1,-TIME,-DV,-OCC)), by = c("id" = "ID")) %>%
  #combine subj info
  mutate(AUC_eta = DOSE/Cl*1000)#%>% #1000*mg/(L/hr) = ng/mL/hr
  #left_join(BW_June2022_ave, by = "SUBJID") %>%
  #mutate(Asian = if_else(SUBJID %in% asianSubj.Aug2022, 1,0)) #add in Asian marker as per communication with Jane

##3.1 CL vs Dose-----
ggplot(data = CLposthoc_clean.Nov2022)+
  geom_point(aes(
    x= factor(DOSE, labels = c("0.2", "0.4", "0.8", "1.2", "1.6", "2.0", "2.5","3.5")), 
    y = Cl, 
    color = factor(PEMBRO, labels = c("No","Yes"))),size = 4)+
  
  
  labs(x = "DOSE (mg)", y = "CL (L/hr)", color = "PEMBRO")+
  
  theme_bw()

ggsave("./output/20221115_NLME_popPK/CLeta_vs_DOSE.png", width = 5, height = 4, dpi = "print", bg = "transparent" )  

##3.2 AUCeta vs Dose-----

filter(select(CLposthoc_clean.Nov2022,Cl, SUBJID,DOSE,PEMBRO,AUC_eta), SUBJID %in% asianSubj.Aug2022)
#         Cl    SUBJID DOSE PEMBRO   AUC_eta
# 1 59.75311 58006-003  0.4      1  6.694212
# 2 46.40389 07001-003  1.2      1 25.859903
# 3 54.69503 58001-004  2.0      0 36.566396
# 4 47.69738 58004-006  2.5      0 52.413780

#Mean Cl 53.14
#Median Cl 51.20

ggplot(data = CLposthoc_clean.Nov2022)+
  geom_point(aes(
    x= factor(DOSE, labels = c("0.2", "0.4", "0.8", "1.2", "1.6", "2.0", "2.5","3.5")), 
    y =AUC_eta, 
    color = factor(PEMBRO, labels = c("No","Yes"))),size = 4)+
  
  
  labs(x = "DOSE (mg)", 
       y = expression(AUC[inf]~"(ng/mL*hr)"), 
       color = "PEMBRO")+
  
  theme_bw()

ggsave("./output/20221115_NLME_popPK/AUCeta_vs_DOSE.png", width = 5, height = 4, dpi = "print", bg = "transparent" ) 
