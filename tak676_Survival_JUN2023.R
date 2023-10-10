library(tidyverse)
library(ggplot2)
library(magrittr)
library(stringr)
library(scales)
library(lubridate)
library(survival)
library(survminer)


#AIM: Survival analysis stratified by Cmax, with time to discontinuation as end point




#_____________________----- 
#0. Import data-----


##Disposition----

May_2023_SS <- read.csv('./Data/20230607_Survival/DS_TAK-676-1002_Dynamic Listings_01May2023.csv',na.string = "n/a")
#`Cytokine.Release.Syndrome.Flag..CRSFL.` is a good variable to use

#Apr_2023_PC <- read.csv('./Data/20230404_PK/TAK-676-1002_Dynamic Listings_03-Apr-2023.csv',na.string = "n/a")

##Demographics for deriving start date

May_2023_DM <- read.csv('./Data/20230607_Survival/DM_TAK-676-1002_Dynamic Listings_01May2023.csv',na.string = "n/a")


##PK----
cmax.05042023.SDTM<- 
  May_2023_PC %>%
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
  left_join(May_2023_cohort_std.SDTM, by = "SUBJID") %>% 
  arrange(PEMBRO,DOSE)



#_____________________----- 
#1. Construct Survival------




##Derive Start Date----

May_2023_DM.startDate <- May_2023_DM%>%
  mutate(SUBJID = str_sub(Unique.Subject.Identifier..USUBJID., start = 14L))%>%
  filter(Date.Time.of.First.Study.Treatment..RFXSTDTC.!="")%>%
  mutate(startDate = ymd_hm(Date.Time.of.First.Study.Treatment..RFXSTDTC.)) %>%
  mutate(timeTillMay1 = int_length(interval(startDate,ymd_hm("2023-05-01T00:00"))) %/% ddays(1)) %>%
  select(SUBJID,timeTillMay1)
  
  


May_2023_SS.transformed <- May_2023_SS %>%
  filter(Subcategory.for.Disposition.Event..DSSCAT. %in% c("END OF STUDY","END OF TREATMENT")) %>%
  mutate(SUBJID = str_sub(Unique.Subject.Identifier..USUBJID., start = 14L)) %>%
  filter(!is.na(Study.Day.of.Start.of.Disposition.Event..DSSTDY.)) %>%
  group_by(SUBJID) %>% #end of treatment vs end of study
  summarise(time = min(Study.Day.of.Start.of.Disposition.Event..DSSTDY., na.rm = TRUE)) #Picking the smaller one between EOS and EOT (usually treatment ends earlier), if there was a crossover, the time crossed was selected
  


##Survival dataset----
May_2023_SS.transformed.survival <- cmax.05042023.SDTM %>%
  left_join(May_2023_SS.transformed, by = "SUBJID") %>%
  
  left_join(May_2023_DM.startDate, by = "SUBJID")%>%
  
  mutate(censor = if_else(is.na(time), 0, 1)) %>% #need to consider this
  
  mutate(TTD = if_else(is.na(time),timeTillMay1, time)) %>%
  
  filter(!is.na(Cmax.mean))



May_2023_SS.transformed.survival.qrt <-May_2023_SS.transformed.survival%>%
  ungroup()%>%
  arrange(Cmax.mean)%>%
  mutate(Cmax.qrt = case_when(
    row_number(Cmax.mean) %in% 1:(floor(nrow(May_2023_SS.transformed.survival)/4)) ~ "1st qrt" ,
    row_number(Cmax.mean) %in% (1+(floor(nrow(May_2023_SS.transformed.survival)/4))):(2*floor(nrow(May_2023_SS.transformed.survival)/4)) ~"2nd qrt",
    row_number(Cmax.mean) %in% (1+(2*floor(nrow(May_2023_SS.transformed.survival)/4))):(3*floor(nrow(May_2023_SS.transformed.survival)/4)) ~"3rd qrt",
    row_number(Cmax.mean) %in% (1+(3*floor(nrow(May_2023_SS.transformed.survival)/4))):(4*floor(nrow(May_2023_SS.transformed.survival)/4)) ~"4th qrt"
  ))
  
  
  






#2. Survival analysis------


##2.0 LUMPED------
surv_object0 <- Surv(time = May_2023_SS.transformed.survival$TTD, event = May_2023_SS.transformed.survival$censor)

fit0 <- survfit(surv_object ~ 1, data = May_2023_SS.transformed.survival)
surv_median(fit0)
#   strata median lower upper
# 1    All     61    58    65

ggsurvplot(fit0, data = May_2023_SS.transformed.survival, pval = TRUE, conf.int = TRUE, pval.method=TRUE)







##2.1 By PEMBRO------
surv_object <- Surv(time = May_2023_SS.transformed.survival$TTD, event = May_2023_SS.transformed.survival$censor)

fit1 <- survfit(surv_object ~ PEMBRO, data = May_2023_SS.transformed.survival)
surv_median(fit1)
# strata median lower upper
# 1 PEMBRO=N     58    51    70
# 2 PEMBRO=Y     63    59    99



ggsurvplot(fit1, data = May_2023_SS.transformed.survival, pval = TRUE, conf.int = TRUE, pval.method=TRUE)


ggsave( "./output/20230607_Survival/Survival_treatment.png", width = 8, height = 6, dpi = "print", bg = "transparent")









###2.1.1 By PEMBRO - COXPH------

cox1 <- coxph(surv_object ~ PEMBRO, data = May_2023_SS.transformed.survival)
summary(cox1)
# n= 80, number of events= 68 
# 
#            coef exp(coef) se(coef)     z Pr(>|z|)
# PEMBROY -0.1836    0.8323   0.2448 -0.75    0.453
# 
# exp(coef) exp(-coef) lower .95 upper .95
# PEMBROY    0.8323      1.201    0.5151     1.345
# 
# Concordance= 0.537  (se = 0.034 )
# Likelihood ratio test= 0.56  on 1 df,   p=0.5
# Wald test            = 0.56  on 1 df,   p=0.5
# Score (logrank) test = 0.56  on 1 df,   p=0.5

#Slightly negative, PEMBRO decrease the risk - but not significantly

PEMBRO_df <- with(May_2023_SS.transformed.survival,
               data.frame(PEMBRO = c("N", "Y")
               )
)


fit.cox.1 <-  survfit(cox1, newdata = PEMBRO_df)

ggsurvplot(fit.cox.1, data = May_2023_SS.transformed.survival, pval = TRUE, conf.int = TRUE, legend.labs=c("PEMBRO = N", "PEMBRO = Y"))




##2.2 COMBINED by Cmax------

surv_object2 <- Surv(time = May_2023_SS.transformed.survival.qrt$TTD, event = May_2023_SS.transformed.survival.qrt$censor)

fit2 <- survfit(surv_object2 ~ Cmax.qrt, data = May_2023_SS.transformed.survival.qrt)
surv_median(fit2)
#             strata median lower upper
# 1 Cmax.qrt=1st qrt   47.5    42    66
# 2 Cmax.qrt=2nd qrt   77.5    61   186
# 3 Cmax.qrt=3rd qrt   60.0    53   121
# 4 Cmax.qrt=4th qrt   60.0    52    NA

ggsurvplot(fit2, data = May_2023_SS.transformed.survival.qrt, pval = TRUE, pval.method=TRUE, test.for.trend = TRUE)

ggsave("./output/20230607_Survival/Survival_Cmax.png", width = 8, height = 6, dpi = "print", bg = "transparent")

write.csv(May_2023_SS.transformed.survival.qrt, "./output/20230607_Survival/May_2023_SS.transformed.survival.qrt.csv", row.names = FALSE)
          




###2.2.1 By Cmax - COXPH------

cox2 <- coxph(surv_object ~ Cmax.qrt, data = May_2023_SS.transformed.survival.qrt)
summary(cox2)


# n= 80, number of events= 68 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)  
# Cmax.qrt2nd qrt -0.01230   0.98777  0.33964 -0.036   0.9711  
# Cmax.qrt3rd qrt  0.07235   1.07503  0.32153  0.225   0.8220  
# Cmax.qrt4th qrt -0.64138   0.52657  0.35720 -1.796   0.0726 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# Cmax.qrt2nd qrt    0.9878     1.0124    0.5076     1.922
# Cmax.qrt3rd qrt    1.0750     0.9302    0.5724     2.019
# Cmax.qrt4th qrt    0.5266     1.8991    0.2615     1.060

PEMBRO_df <- with(May_2023_SS.transformed.survival,
                  data.frame(PEMBRO = c("N", "Y")
                  )
)


fit.cox.1 <-  survfit(cox1, newdata = PEMBRO_df)

ggsurvplot(fit.cox.1, data = May_2023_SS.transformed.survival, pval = TRUE, conf.int = TRUE, legend.labs=c("PEMBRO = N", "PEMBRO = Y"))