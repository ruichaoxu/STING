library(tidyverse)
library(ggplot2)
library(magrittr)
library(stringr)

#0. File-----
#PK file:JAN_PK_2023.R
cmax.01062023.SDTM<- 
  Jan_2023_PC %>%
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
  left_join(Jan_2023_cohort_std.SDTM, by = "SUBJID") %>% 
  arrange(PEMBRO,DOSE)







#PD file:tak-676_biomarkerFishing.R
##IFNG=====
IFNG.data.02012023.raw <- filter(Feb_2023_cytokine.filtered.IFNG,
                                 VISIT %in% c("Cycle 1 Day 1", "Cycle 1 Day 4", "Cycle 1 Day 8","Cycle 1 Day 9" ,"Cycle 1 Day 15","Cycle 1 Day 16")) %>%
  
  mutate(RESULT.C1D1 = if_else(LBTPT == "PREDOSE", 	
                               LBORRES.num,
                               NA_real_))%>%
  group_by(SUBJID)%>%
  fill(RESULT.C1D1, .direction ="down")%>%
  mutate(FC.C1D1 = LBORRES.num/RESULT.C1D1)%>%
  filter(LBTPT == "24 HOURS POST-EOI")%>% #same method with previous analysis, could do max as well
 
  summarise(IFNG.mean = mean(FC.C1D1, na.rm = TRUE))

#using max FC rather than 24hr fold change
IFNG.data.02012023.raw.2 <- filter(Feb_2023_cytokine.filtered.IFNG,
                                 VISIT %in% c("Cycle 1 Day 1", "Cycle 1 Day 4", "Cycle 1 Day 8","Cycle 1 Day 9" ,"Cycle 1 Day 15","Cycle 1 Day 16")) %>%
  
  mutate(RESULT.C1D1 = if_else(LBTPT == "PREDOSE", 	
                               LBORRES.num,
                               NA_real_))%>%
  group_by(SUBJID)%>%
  fill(RESULT.C1D1, .direction ="down")%>%
  mutate(FC.C1D1 = LBORRES.num/RESULT.C1D1)%>%
 
  
  summarise(IFNG.FC.MAX.C1D1 = max(FC.C1D1))





##3+8+4-Ki67+ (%CD8)======
CD8T.data.01062023.raw<- Jan_2023_flow.filtered.CD8 %>%
  filter(Visit.Names %in% c('C1 D2 24HR POST EOI','C1 D4','C1 D8 PREDOSE','C1 D8 6HR POST EOI','C1 D9 24HR POST EOI','C1 D15 PREDOSE','C1 D16 24HR POST EOI','C1 D1 PREDOSE','C1 D1 6HR POST EOI','C1 D1 PREDOSE AMD4','C1 D1 6HR POST EOI AMD4','C1 D8 PREDOSE AMD4','C1 D15 PREDOSE AMD4','C1 D15 3HR POST EOI','C1 D8 6HR POST EOI AMD4')) %>% #dropping data beyond cycle 1
  mutate(RESULT.C1D1 = if_else(Visit.Names %in% c("C1 D1 PREDOSE AMD4","C1 D1 PREDOSE"),
                                                  Numeric.Result_Finding.in.Standard.Units,
                                                  NA_real_)
    
  )%>%
  group_by(Subject.Identifier.for.the.Study)%>%
  fill(RESULT.C1D1, .direction ="down") %>%
  mutate(FC.C1D1 = Numeric.Result_Finding.in.Standard.Units/RESULT.C1D1) %>%
  summarise(CD8T.FC.MAX.C1D1 = max(FC.C1D1)) %>%
  mutate(SUBJID = str_sub( Subject.Identifier.for.the.Study,start = 2L),.keep = "unused")
  


# CD8Ki67T.data.10252022.SDTM  <-filter(lab.SDTM.source.xaxis,
#                                       Lab.Test.or.Examination.Short.Name..LBTESTCD. == "FCT54899"
#                                       &Visit.Number..VISITNUM. <3 
#                                       &!is.na(xaxis.num)) %>%
#   filter(Planned.Time.Point.Number..LBTPTNUM. == 24)%>%
#   group_by(Unique.Subject.Identifier..USUBJID.) %>%
#   summarise(CD8KI67T.mean = mean(RESULT.FC.C1D1, na.rm = TRUE))





##IP10------

IP10.data.02012023.raw <- filter(Feb_2023_cytokine.filtered.IP10,
                                   VISIT %in% c("Cycle 1 Day 1", "Cycle 1 Day 4", "Cycle 1 Day 8","Cycle 1 Day 9" ,"Cycle 1 Day 15","Cycle 1 Day 16")) %>%
  
  mutate(RESULT.C1D1 = if_else(LBTPT == "PREDOSE", 	
                               LBORRES.num,
                               NA_real_))%>%
  group_by(SUBJID)%>%
  fill(RESULT.C1D1, .direction ="down")%>%
  mutate(FC.C1D1 = LBORRES.num/RESULT.C1D1)%>%
  
  
  summarise(IP10.FC.MAX.C1D1 = max(FC.C1D1))






#1. Plot----

pkpd.plot.data.Feb2023 <- cmax.01062023.SDTM %>%
  left_join(IFNG.data.02012023.raw.2, by = "SUBJID") %>%
  left_join(CD8T.data.01062023.raw, by = "SUBJID" ) %>%
  left_join(IP10.data.02012023.raw, by = "SUBJID")%>%
  #filter(DOSE<5.0) %>%
  mutate(PEMBRO = if_else(PEMBRO == 'N', "1. Single Agent", "2. Combination" ))

write.csv(pkpd.plot.data.Feb2023, "./output/20230220_PKPD/Cmax_PK_PD.csv", row.names = F)



##1.1 Cmax vs INFG----
ggplot(data = filter(pkpd.plot.data.Feb2023,DOSE<7))+
  geom_point(aes(x = Cmax.mean, y = IFNG.FC.MAX.C1D1, shape = PEMBRO, color = factor(DOSE)), size = 4)+
  geom_hline(yintercept = 1.75, linetype = "dashed")+
  geom_smooth(aes(x = Cmax.mean, y = IFNG.FC.MAX.C1D1), method = "lm", se = FALSE,  color = "Black")+
  scale_color_brewer(palette = "Spectral", direction = -1)+
  scale_x_log10()+
  #scale_y_continuous(breaks = c(0,2,seq(5,40,by = 5)), labels = c("0","2",as.character(seq(5,40,by = 5))))+
  labs(color = "Dose (mg)",
       shape = "Treatment Arm",
       
       x = "Mean Cmax (ng/mL)",
       y = "Fold change of INF-G (Maximum VS baseline)")+
  facet_wrap(ncol = 1,vars(PEMBRO))+
  theme(
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.line = element_line()
  ) +
  theme(
    panel.border = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )

ggsave("./output/20230220_PKPD/IFG.max_vs_meanCmax_faceted.png", width = 8, height = 6, dpi = "print", bg = "transparent")


##Coeff----
lm(data = filter(pkpd.plot.data.Feb2023, PEMBRO == "1. Single Agent"), IFNG.mean ~ Cmax.mean)$coefficients[2]
#SLOPE = 0.004935081
lm(data = filter(pkpd.plot.data.Feb2023, PEMBRO == "2. Combination"), IFNG.mean ~ Cmax.mean)$coefficients[2]
#SLOPE = 0.0353907
summary(lm(data = filter(pkpd.plot.data.Feb2023, PEMBRO == "1. Single Agent"),IFNG.mean ~ Cmax.mean))$coefficients[2,4]
#p = 0.695, NS
summary(lm(data = filter(pkpd.plot.data.Feb2023, PEMBRO == "2. Combination"),IFNG.mean ~ Cmax.mean))$coefficients[2,4]
#p = 0.237, NS




##Coeff.2----
lm(data = filter(pkpd.plot.data.Feb2023, PEMBRO == "1. Single Agent"), IFNG.FC.MAX.C1D1 ~ Cmax.mean)$coefficients[2]
#SLOPE = 0.027673 
lm(data = filter(pkpd.plot.data.Feb2023, PEMBRO == "2. Combination"), IFNG.FC.MAX.C1D1 ~ Cmax.mean)$coefficients[2]
#SLOPE = 0.1489111
summary(lm(data = filter(pkpd.plot.data.Feb2023, PEMBRO == "1. Single Agent"),IFNG.FC.MAX.C1D1 ~ Cmax.mean))$coefficients[2,4]
#p = 0.149, NS
#Multiple R-squared:  0.1011, r = 0.32
summary(lm(data = filter(pkpd.plot.data.Feb2023, PEMBRO == "2. Combination"),IFNG.FC.MAX.C1D1 ~ Cmax.mean))$coefficients[2,4]
#p = 0.051, NS
#Multiple R-squared:  0.1438, r = 0.38




##1.2 Cmax vs Ki67CD8----
ggplot(data = filter(pkpd.plot.data.Feb2023,DOSE<7))+
  geom_point(aes(x = Cmax.mean, y = CD8T.FC.MAX.C1D1, shape = PEMBRO, color = factor(DOSE)), size = 4)+
  geom_hline(yintercept = 1.75, linetype = "dashed")+
  geom_smooth(aes(x = Cmax.mean, y = CD8T.FC.MAX.C1D1), method = "lm", se = FALSE,  color = "Black")+
  scale_color_brewer(palette = "Spectral", direction = -1)+
  scale_x_log10()+
  scale_y_continuous(limits= c(0,6), breaks = c(0,1,1.75,seq(2,6,by = 1)), labels = c("0","1","1.75",as.character(seq(2,6,by = 1))))+
  labs(color = "Dose (mg)",
       shape = "Treatment Arm",
       
       x = "Mean Cmax (ng/mL)",
       y = "Fold change of CD8+ Ki67+ T cell (Maximum VS baseline)")+
  facet_wrap(ncol = 1,vars(PEMBRO))+
  theme(
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.line = element_line()
  ) +
  theme(
    panel.border = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )


ggsave("./output/20230220_PKPD/CD8T.max_vs_meanCmax_faceted.png", width = 8, height = 6, dpi = "print", bg = "transparent")



##Coeff----
# lm(data = filter(pkpd.plot.data.Feb2023, PEMBRO == "1. Single Agent"), CD8T.FC.MAX.C1D1 ~ Cmax.mean)$coefficients[2]
# #SLOPE = 0.00714
# lm(data = filter(pkpd.plot.data, PEMBRO == "2. Combination"), CD8KI67T.mean ~ Cmax.mean)$coefficients[2]
# #SLOPE = -0.005292058 

summary(lm(data = filter(pkpd.plot.data.Feb2023, PEMBRO == "2. Combination"), CD8T.FC.MAX.C1D1 ~ Cmax.mean))
#-0.01510, p = 0.415, NS
summary(lm(data = filter(pkpd.plot.data.Feb2023, PEMBRO == "1. Single Agent"), CD8T.FC.MAX.C1D1 ~ Cmax.mean))
#-0.002550, p = 0.636, NS

ggsave("./output/20221111_PKPD/meanKi67CD8_vs_meanCmax_faceted.png", width = 8, height = 6, dpi = "print", bg = "transparent")










##1.3 Cmax vs IP10----

ggplot(data = filter(pkpd.plot.data.Feb2023,DOSE<7))+
  geom_point(aes(x = Cmax.mean, y = IP10.FC.MAX.C1D1, shape = PEMBRO, color = factor(DOSE)), size = 4)+
  #geom_hline(yintercept = 1.75, linetype = "dashed")+ 
  geom_smooth(aes(x = Cmax.mean, y = IP10.FC.MAX.C1D1), method = "lm", se = FALSE,  color = "Black")+
  scale_color_brewer(palette = "Spectral", direction = -1)+
  scale_x_log10()+
  #scale_y_continuous(breaks = c(0,2,seq(5,40,by = 5)), labels = c("0","2",as.character(seq(5,40,by = 5))))+
  labs(color = "Dose (mg)",
       shape = "Treatment Arm",
       
       x = "Mean Cmax (ng/mL)",
       y = "Fold change of IP-10 (Maximum VS baseline)")+
  facet_wrap(ncol = 1,vars(PEMBRO))+
  theme(
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.line = element_line()
  ) +
  theme(
    panel.border = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )

ggsave("./output/20230220_PKPD/IP10.max_vs_meanCmax_faceted.png", width = 8, height = 6, dpi = "print", bg = "transparent")


##Coeff----
#lm(data = filter(pkpd.plot.data.Feb2023, PEMBRO == "1. Single Agent"), IFNG.mean ~ Cmax.mean)$coefficients[2]
#SLOPE = 0.004935081
#lm(data = filter(pkpd.plot.data.Feb2023, PEMBRO == "2. Combination"), IFNG.mean ~ Cmax.mean)$coefficients[2]
#SLOPE = 0.0353907
summary(lm(data = filter(pkpd.plot.data.Feb2023, PEMBRO == "1. Single Agent"),IP10.FC.MAX.C1D1 ~ Cmax.mean))
#0.017, p = 0.056, NS
#Multiple R-squared:  0.133, r = 0.36
summary(lm(data = filter(pkpd.plot.data.Feb2023, PEMBRO == "2. Combination"),IP10.FC.MAX.C1D1 ~ Cmax.mean))
#0.024, p = 0.051, NS
#Multiple R-squared:  0.01593, r = 0.13