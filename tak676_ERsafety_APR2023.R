library(tidyverse)
library(ggplot2)
library(magrittr)
library(stringr)
library(scales)


#AIM: Logistic regression for Cmax-CRS

#_____________________----- 
#0. Import data-----

Apr_2023_AE <- read.csv('./Data/20230424_ERsafety/AE_TAK-676-1002_Dynamic Listings_03-Apr-2023.csv',na.string = "n/a")
#`Cytokine.Release.Syndrome.Flag..CRSFL.` is a good variable to use

#Apr_2023_PC <- read.csv('./Data/20230404_PK/TAK-676-1002_Dynamic Listings_03-Apr-2023.csv',na.string = "n/a")





#_____________________----- 
#1. Construct Cmax------

cmax.Apr2023.SDTM<- 
  Apr_2023_PC %>%
  mutate(SUBJID = str_sub(Unique.Subject.Identifier..USUBJID., start = 14L),
         TIME = case_when(
           Planned.Time.Point.Number..PCTPTNUM. == -2 ~ 0,
           Planned.Time.Point.Number..PCTPTNUM. == 0 ~ 1,
           Planned.Time.Point.Number..PCTPTNUM. == 0.5 ~ 1.5,
           Planned.Time.Point.Number..PCTPTNUM. == 1 ~ 2,
           Planned.Time.Point.Number..PCTPTNUM. == 2 ~ 3,
           Planned.Time.Point.Number..PCTPTNUM. == 3 ~ 4,
           Planned.Time.Point.Number..PCTPTNUM. == 6 ~ 7,
           Planned.Time.Point.Number..PCTPTNUM. == 24 ~ 25,
           TRUE ~ NA_real_
         ),
         CONC = as.numeric(Numeric.Result.Finding.in.Standard.Units..PCSTRESN.)) %>%
  filter(Specimen.Material.Type..PCSPEC. == "PLASMA") %>%
  filter(Planned.Time.Point.Number..PCTPTNUM. == 0) %>%
  filter(Visit.Number..VISITNUM. %in% c(1.01,1.08,1.15)) %>%
  group_by(SUBJID,Unique.Subject.Identifier..USUBJID.) %>%
  summarise(Cmax.mean = mean(CONC,na.rm = TRUE )) %>%
  left_join(Apr_2023_cohort_std.SDTM, by = "SUBJID") %>% 
  arrange(PEMBRO,DOSE)








#_____________________----- 
#2. Construct safety------
#513 reportable terms
#353 LOWEST LEVEL TERMS - seems ok to work with




##2.1 Construct binary safety------
Apr_2023_AE_Binary <- Apr_2023_AE %>%
  mutate(SUBJID = str_sub(Unique.Subject.Identifier..USUBJID., start = 14L),
         seen = 1)%>%
  filter(`Lowest.Level.Term..AELLT.` != "") %>%
  select(
    SUBJID, 
    `Lowest.Level.Term..AELLT.`,
    `Standard.Toxicity.Grade..AETOXGR.`,
    seen,
    `Start.Date.Time.of.Adverse.Event..AESTDTC.`)%>%
  distinct()%>% #duplicated records (may contain useful info) were truncated
  pivot_wider(
    names_from = `Lowest.Level.Term..AELLT.`,
    values_from = seen,
    values_fill = 0)




##2.2 Construct specific safety marker dataset------



#CRS as an example

###2.2.1 Safety with DATE-------
#CRS only
Apr_2023_AE_CRS <- Apr_2023_AE_Binary %>%
  select(SUBJID,
         `Cytokine release syndrome`,
         `Standard.Toxicity.Grade..AETOXGR.`,
         `Start.Date.Time.of.Adverse.Event..AESTDTC.`) %>%
  filter(`Cytokine release syndrome` != 0)

#CRS or fever
Apr_2023_AE_CRSorFever <- Apr_2023_AE_Binary %>%
  select(SUBJID,
         `Cytokine release syndrome`,
         `Fever`,
         `Standard.Toxicity.Grade..AETOXGR.`,
         `Start.Date.Time.of.Adverse.Event..AESTDTC.`) %>%
  filter(`Cytokine release syndrome`|Fever != 0)

###2.2.2 unique safety per subject-------

#CRS alone
Apr_2023_AE_CRS.Unique <- Apr_2023_AE_CRS %>%
  group_by(SUBJID) %>%
  summarize(CRS.marker = 1,
            count = n(),
            max.grade = max(`Standard.Toxicity.Grade..AETOXGR.`))

#CRS or fever
Apr_2023_AE_CRSorFever.Unique <- Apr_2023_AE_CRSorFever %>%
  group_by(SUBJID) %>%
  summarize(CRSorFever.marker = 1,
            count = n(),
            max.grade = max(`Standard.Toxicity.Grade..AETOXGR.`))




#_____________________----- 
#3. PK-Safety------




##3.1 CRS vs Cmax-----

###3.1.1 dataset -----
#for dot
Apr_2023_AE_CRS.Unique.PK <- cmax.Apr2023.SDTM %>%
  left_join(Apr_2023_AE_CRS.Unique, by = "SUBJID") %>%
  mutate(CRS.marker = if_else(is.na(CRS.marker),0, CRS.marker)) %>%
  ungroup()

#for quartile
Apr_2023_AE_CRS.Unique.PK.qrt <- Apr_2023_AE_CRS.Unique.PK %>%
  select(Cmax.mean, CRS.marker) %>%
  drop_na()%>%
  arrange(Cmax.mean)%>%
  ungroup()%>%
  mutate(qrt = 0) %>%
  mutate(qrt = case_when(
    row_number(qrt) %in% 1:(floor(nrow(Apr_2023_AE_CRS.Unique.PK)/4)) ~ 1,
    row_number(qrt) %in% (floor((nrow(Apr_2023_AE_CRS.Unique.PK)/4))+1):(2*floor((nrow(Apr_2023_AE_CRS.Unique.PK)/4))) ~2,
    row_number(qrt) %in% (2*floor((nrow(Apr_2023_AE_CRS.Unique.PK)/4))+1):(3*floor((nrow(Apr_2023_AE_CRS.Unique.PK)/4)))~3,
    TRUE~4
  )) %>%
  group_by(qrt)%>%
  summarize(Cmax.qrt = mean(Cmax.mean),
            freq.qrt = mean(CRS.marker),
            count = n(),
            freq.se = mean(CRS.marker)/sqrt(n())) 

###3.1.2 Plot----
ggplot()+
  geom_point(data = Apr_2023_AE_CRS.Unique.PK, aes(x = Cmax.mean, y = CRS.marker, color = PEMBRO), position = position_jitter(width=0, height=0.005))+
  geom_smooth(
    data = Apr_2023_AE_CRS.Unique.PK, aes(x = Cmax.mean, y = CRS.marker),
    method="glm", method.args = list(family= "binomial"))+
  geom_point(
    data = Apr_2023_AE_CRS.Unique.PK.qrt, 
    aes(x = Cmax.qrt, y = freq.qrt),
    size = 3, shape = 0)+
  geom_linerange(
    data = Apr_2023_AE_CRS.Unique.PK.qrt, 
    aes(x = Cmax.qrt, ymin = freq.qrt-freq.se, ymax = freq.qrt+freq.se)
  )+
  scale_x_log10()+
  
  labs(y = "Any Grade CRS")+
  
  theme_bw()+
  geom_vline(xintercept = 
               c(
    arrange(Apr_2023_AE_CRS.Unique.PK,Cmax.mean)$Cmax.mean[floor(nrow(Apr_2023_AE_CRS.Unique.PK)/4)],
    arrange(Apr_2023_AE_CRS.Unique.PK,Cmax.mean)$Cmax.mean[2*floor(nrow(Apr_2023_AE_CRS.Unique.PK)/4)],
    arrange(Apr_2023_AE_CRS.Unique.PK,Cmax.mean)$Cmax.mean[3*floor(nrow(Apr_2023_AE_CRS.Unique.PK)/4)]
    ),
    linetype = "dashed"
    )

###3.1.3 STAT ----
summary(glm( CRS.marker ~ Cmax.mean, data = Apr_2023_AE_CRS.Unique.PK, family = binomial))
# Coefficients:
#                Estimate Std. Error z value Pr(>|z|)    
#   (Intercept) -3.317058   0.652318  -5.085 3.68e-07 ***
#   Cmax.mean    0.024953   0.008409   2.968    0.003 ** 


ggsave("./output/20230424_ERsafety/any Grade CRS vs Cmax.png", width = 8, height = 6, dpi = "print", bg = "transparent")
  










##3.2 CRS|Fever vs Cmax-----

###3.2.1 dataset -----
#for dot
Apr_2023_AE_CRSorFever.Unique.PK <- cmax.Apr2023.SDTM %>%
  left_join(Apr_2023_AE_CRSorFever.Unique, by = "SUBJID") %>%
  mutate(CRSorFever.marker = if_else(is.na(CRSorFever.marker),0, CRSorFever.marker)) %>%
  ungroup()

#for quartile
Apr_2023_AE_CRSorFever.Unique.PK.qrt <- Apr_2023_AE_CRSorFever.Unique.PK %>%
  select(Cmax.mean, CRSorFever.marker) %>%
  drop_na()%>%
  arrange(Cmax.mean)%>%
  ungroup()%>%
  mutate(qrt = 0) %>%
  mutate(qrt = case_when(
    row_number(qrt) %in% 1:(floor(nrow(Apr_2023_AE_CRSorFever.Unique.PK)/4)) ~ 1,
    row_number(qrt) %in% (floor((nrow(Apr_2023_AE_CRSorFever.Unique.PK)/4))+1):(2*floor((nrow(Apr_2023_AE_CRSorFever.Unique.PK)/4))) ~2,
    row_number(qrt) %in% (2*floor((nrow(Apr_2023_AE_CRSorFever.Unique.PK)/4))+1):(3*floor((nrow(Apr_2023_AE_CRSorFever.Unique.PK)/4)))~3,
    TRUE~4
  )) %>%
  group_by(qrt)%>%
  summarize(Cmax.qrt = mean(Cmax.mean),
            freq.qrt = mean(CRSorFever.marker),
            count = n(),
            freq.se = mean(CRSorFever.marker)/sqrt(n())) 

###3.2.2 Plot----
ggplot()+
  geom_point(data = Apr_2023_AE_CRSorFever.Unique.PK, aes(x = Cmax.mean, y = CRSorFever.marker, color = PEMBRO), position = position_jitter(width=0, height=0.005))+
  geom_smooth(
    data = Apr_2023_AE_CRSorFever.Unique.PK, aes(x = Cmax.mean, y = CRSorFever.marker),
    method="glm", method.args = list(family= "binomial"))+
  geom_point(
    data = Apr_2023_AE_CRSorFever.Unique.PK.qrt, 
    aes(x = Cmax.qrt, y = freq.qrt),
    size = 3, shape = 0)+
  geom_linerange(
    data = Apr_2023_AE_CRSorFever.Unique.PK.qrt, 
    aes(x = Cmax.qrt, ymin = freq.qrt-freq.se, ymax = freq.qrt+freq.se)
  )+
  scale_x_log10()+
  
  labs(y = "Any Grade CRS or Fever")+
  
  theme_bw()+
  geom_vline(xintercept = 
               c(
                 arrange(Apr_2023_AE_CRS.Unique.PK,Cmax.mean)$Cmax.mean[floor(nrow(Apr_2023_AE_CRS.Unique.PK)/4)],
                 arrange(Apr_2023_AE_CRS.Unique.PK,Cmax.mean)$Cmax.mean[2*floor(nrow(Apr_2023_AE_CRS.Unique.PK)/4)],
                 arrange(Apr_2023_AE_CRS.Unique.PK,Cmax.mean)$Cmax.mean[3*floor(nrow(Apr_2023_AE_CRS.Unique.PK)/4)]
               ),
             linetype = "dashed"
  )

###3.2.3 STAT ----
summary(glm( CRSorFever.marker ~ Cmax.mean, data = Apr_2023_AE_CRSorFever.Unique.PK, family = binomial))
# Coefficients:
#                Estimate Std. Error z value Pr(>|z|)    
#   (Intercept) -1.673011   0.409445  -4.086 4.39e-05 ***
#   Cmax.mean    0.021543   0.007666   2.810  0.00495 ** 


ggsave("./output/20230424_ERsafety/any Grade CRS or Fever vs Cmax.png", width = 8, height = 6, dpi = "print", bg = "transparent")







###Additional Plots----

#CRS or Fever ~ PEMBRO
ggplot()+
  geom_point(data = Apr_2023_AE_CRSorFever.Unique.PK, aes(x = Cmax.mean, y = CRSorFever.marker, color = PEMBRO), position = position_jitter(width=0, height=0.005))+
  geom_smooth(
    data = filter(Apr_2023_AE_CRSorFever.Unique.PK,PEMBRO=="Y"), aes(x = Cmax.mean, y = CRSorFever.marker),
    method="glm", method.args = list(family= "binomial"),
    color ="#00BFC4" )+
  geom_smooth(
    data = filter(Apr_2023_AE_CRSorFever.Unique.PK,PEMBRO=="N"), aes(x = Cmax.mean, y = CRSorFever.marker),
    method="glm", method.args = list(family= "binomial"),
    color = "#F8766D")+
  scale_x_log10()+
  
  theme_bw()+
  labs(y = "Any Grade CRS or Fever")

ggsave("./output/20230424_ERsafety/any Grade CRS or Fever vs Cmax_PEMBRO.png", width = 8, height = 6, dpi = "print", bg = "transparent")

summary(glm( CRSorFever.marker ~ Cmax.mean, data = filter(Apr_2023_AE_CRSorFever.Unique.PK, PEMBRO == "Y"), family = binomial))
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)   
#   (Intercept) -2.16735    0.73865  -2.934  0.00334 **
#   Cmax.mean    0.04214    0.01839   2.291  0.02195 * 
#3.5mg
expit(-2.16735+0.04214*(52+46)/2) 
#47%

#5mg
expit(-2.16735+0.04214*(52+46)/2*5/3.5)
#69%

#10mg
expit(-2.16735+0.04214*(52+46)/2*10/3.5)
#97%


summary(glm( CRSorFever.marker ~ Cmax.mean, data = filter(Apr_2023_AE_CRSorFever.Unique.PK, PEMBRO == "N"), family = binomial))
# Coefficients:
#                Estimate Std. Error z value Pr(>|z|)   
#   (Intercept) -1.748699   0.594358  -2.942  0.00326 **
#   Cmax.mean    0.017938   0.008349   2.149  0.03167 * 



#3.5mg
expit(-1.748699+0.017938*(51+40)/2) 
#28%

#5mg
expit(-1.748699+0.017938*(51+40)/2*5/3.5)
#36%

#10mg
expit(-1.748699+0.017938*(51+40)/2*10/3.5)
#64%




#CRS ~ PEMBRO






#others----
## 95% CI-----
logit.predict.footlength <- seq(min(drop_na(Apr_2023_AE_CRS.Unique.PK,Cmax.mean)$Cmax.mean),max(drop_na(Apr_2023_AE_CRS.Unique.PK,Cmax.mean)$Cmax.mean),by=0.5)

logit.pred.test <-tibble(
  x = logit.predict.footlength,
  y = expit(-3.317058+0.024953*logit.predict.footlength),
  y95 = expit(-3.317058+1.96*0.652318 +(0.024953+1.96*0.008409)*logit.predict.footlength),
  y5 = expit(-3.317058-1.96*0.652318+(0.024953-1.96*0.008409)*logit.predict.footlength)
)

ggplot()+
  geom_ribbon(data =logit.pred.test, aes(x = x, ymin = y5, ymax = y95), alpha = 0.2)+
  geom_line(data =logit.pred.test, aes(x = x, y = y) )+
  geom_line(data =logit.pred.test, aes(x = x, y = y95), linetype="dashed" )+
  geom_line(data =logit.pred.test, aes(x = x, y = y5), linetype="dashed" )+
  
  scale_x_log10()
