library(tidyverse)
library(ggplot2)
library(ggpattern)
library(magrittr)
library(stringr)
library(datapasta) #get string expression from vector

#_____________________----- 
#0. Import data-----

geneSig0802_2022 <- read.csv("./Data/20220802_GeneSig/TAK-676-1002_STING24_GMScore07292022.csv")
glimpse(geneSig0802_2022)

geneSig0802_2022%>%
  group_by(Timepoint) %>%
  summarise(Counts = n())

#Change from last file 0610 (shared on July28)
# <chr>      <int>
# 1 24Hr         128 (+15)
# 2 3Hr          132 (+13)
# 3 6Hr          170 (+20)
# 4 Pre          219 (+19)


# 1 24Hr         113 (+12)
# 2 3Hr          119
# 3 6Hr          150
# 4 Pre          200 (+24)

unique(geneSig0802_2022$Dose)

unique(geneSig0324$Dose)


#1. Data Wrangling-----

##1.1 PEMBRO-------

#Extract PEMBRO list
#method 1 - summarise, groupby, slice
pembroSummary <- geneSig0727_2022 %>%
  summarise(SUBJID, PEMBRO) %>%
  group_by(SUBJID) %>%
  slice(1) #utilizing previous gs data PEMBRO info to construct a PEMBRO worksheet

#method2 - select, distinct (ungroup in the beginning may be necessary to avoid attachment of other grouping variable)
pembroSummary2 <- forNCA_July_2022 %>%
  ungroup() %>% #prevet grouping variable be added automatically
  select(SUBJID, PEMBRO) %>% #extract all PEMBRO status
  distinct(SUBJID, PEMBRO)  #Take unique SUBJID number

pembroSummary_comb <- pembroSummary %>%
  bind_rows(pembroSummary2) %>%
  distinct(SUBJID,PEMBRO)
  
  
geneSig0802_2022 <- geneSig0802_2022 %>%
  left_join (pembroSummary_comb, by = "SUBJID") #Fill in PEMBRO info

ggplot(data = geneSig0802_2022)
  
  
  

##1.2.Get MaxFC -----

#get subject max FC compared to C1D1
summarize_geneSig0802_2022 <- geneSig0802_2022%>%
  group_by(PEMBRO,Dose,SUBJID) %>%
  summarise(MaxFC_C1D1 = max(FC.C1D1))

#Another way to get subject max FC compared to C1D1, sorted for more natural order layout
reorder_geneSig0802_2022 <- geneSig0802_2022 %>%
  arrange(PEMBRO,Dose,SUBJID,Cycle,Day)


summarize_reorder_geneSig0802_2022 <- reorder_geneSig0802_2022%>%
  ungroup%>%
  group_by(PEMBRO,Dose,SUBJID) %>%
  summarise(MaxFC_C1D1 = max(FC.C1D1))

#Some subjects with NA values are excluded. NA maybe they do not have pre dose level, so no fold change can be calculated

dropNA_summarize_reorder_geneSig0802_2022 <-summarize_reorder_geneSig0802_2022%>%
  drop_na()

#join 2 data sets - using the filtered data set to filter out the full data set
filtered_summarize_geneSig0802_2022 <- geneSig0802_2022 %>%
  semi_join(dropNA_summarize_reorder_geneSig0802_2022,by = c("SUBJID","FC.C1D1" = "MaxFC_C1D1"))


write.csv(filtered_summarize_geneSig0802_2022,"./output/20220802_Genesig/dropNA_summarize_reorder_geneSig0802_2022.csv", row.names = FALSE)
# colname is still FC.C1D1 but it was actually recording MaxFC_C1D1







#2. PLOT- DOSE - Max FC_C1D1-----

length(unique(geneSig0802_2022$SUBJID)) #total 50 subj
setdiff( unique(geneSig0802_2022$SUBJID),unique(filtered_summarize_geneSig0802_2022$SUBJID))

# 6 subjects in genesig dataset but not in this FCmax dataset
# Mainly due to no baseline level to calculate FC.CxDx for (later on time point could present)
# [1] "58006-004" "58005-004" "58007-009", (previously not in)
#"58007-012" "07001-007" "58004-008" 

##2.1 Correct-24hr only-jitter -Max FC_C1D1-----


ggplot(data =filter(filtered_summarize_geneSig0802_2022,Timepoint=="24Hr"),
       aes(x=factor(Dose),
           y=FC.C1D1)) +
  geom_jitter(alpha = 0.8, width = 0.05,size = 4,aes(color = factor(Dose),shape = factor(Timepoint)))+
  scale_color_brewer(palette = "Reds") +
  labs(color = "DOSE(mg)", shape = "Time Point")+
  scale_y_log10(breaks = c(1,3,10,30,1.56), labels = c("1","3","10","30","1.56")) +
  xlab("Dose (mg)") +
  ylab("Fold Change of gene signature score")+
  ggtitle("24Hr fold change comparing to C1D1 baseline\n(July.28 cutoff)")+
  scale_shape_manual(values = 15)+
  
  
  theme(
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    
    axis.line = element_line()
  ) +
  
  theme(
    panel.border = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    # panel.grid.major = element_line(color = "grey"), # get rid of major grid
    # panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.key = element_rect(fill = "transparent")
  )+
  facet_wrap(~ PEMBRO, 
             labeller = label_both,
             scales = "free_x") +
  geom_hline(yintercept=1.56, linetype="dashed", color = "black")

ggsave("./output/20220802_Genesig/24hr_MaxFC_vs_Dose_July.28.2022.png", width = 8, height = 6, dpi = "print", bg = "transparent")






#_____________________----- 
#3. Appending Cmax-----

##3.1. Plotting TIME vs GeneSig ALL-----




plot_geneSig0802_2022 <- geneSig0802_2022 %>%  
  
  unite(col = "occasion",Visit, Timepoint,sep="_",remove=FALSE)%>% #create x axis
  unite(col = "lineGraphGroup",SUBJID, Visit,remove=FALSE) %>%
  filter(Visit != "C9D1") #irrelevant




###Y-2.0-----

pl.y.20<-
  ggplot(data = filter(plot_geneSig0802_2022,
                       PEMBRO == "Y",
                       Dose == 2.0),
         aes(x=factor(occasion, 
                      levels = c("C1D1_Pre",'C1D1_3Hr',"C1D1_6Hr","C1D1_24Hr",
                                 "C1D8_Pre","C1D8_3Hr","C1D8_6Hr", "C1D8_24Hr",
                                 "C1D15_Pre","C1D15_3Hr","C1D15_24Hr",
                                 "C2D1_Pre","C2D1_6Hr","C3D1_Pre","C3D1_6Hr",
                                 "C4D1_Pre","C4D1_6Hr")),
             y=STING24GMScore))+
  geom_point(aes(color = SUBJID ))+
  
  scale_y_log10(limits = c(NA,500))+
  scale_shape_manual(values = c(16,1))+
  geom_line(aes(group =lineGraphGroup,color = SUBJID)) +
  xlab("occasion")+
  ylab("STING 24 GM Score")+
  labs(color = "Subject")+
  #facet_wrap(~Visit,nrow=1, strip.position = "bottom") +
  
  theme(
    
    panel.spacing.x =unit(0, "lines")
    ,strip.background = element_blank()
    ,panel.border = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    #plot.background = element_rect(fill = "transparent", color = "black"), # bg of the plot
    #panel.grid.major = element_line(color = "grey"), # get rid of major grid
    #panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    
    axis.text.x = element_text(angle = 45,vjust = 1,hjust=1)
  )+
  geom_vline(xintercept=c(4.5,8.5,11.5,13.5), linetype="dashed", color = "black")

###export-----
ggsave("./output/20220802_GeneSig/timeVSgeneSig/Y_2.0.png", width = 8, height = 2, dpi = "print", bg = "transparent")







###N-2.5-----

pl.n.25<-
  ggplot(data = filter(plot_geneSig0802_2022,
                       PEMBRO == "N",
                       Dose == 2.5),
         aes(x=factor(occasion, 
                      levels = c("C1D1_Pre",'C1D1_3Hr',"C1D1_6Hr","C1D1_24Hr",
                                 "C1D8_Pre","C1D8_3Hr","C1D8_6Hr", "C1D8_24Hr",
                                 "C1D15_Pre","C1D15_3Hr","C1D15_24Hr",
                                 "C2D1_Pre","C2D1_6Hr","C3D1_Pre","C3D1_6Hr",
                                 "C4D1_Pre","C4D1_6Hr")),
             y=STING24GMScore))+
  geom_point(aes(color = SUBJID ))+
  
  scale_y_log10(limits = c(NA,500))+
  scale_shape_manual(values = c(16,1))+
  geom_line(aes(group =lineGraphGroup,color = SUBJID)) +
  xlab("occasion")+
  ylab("STING 24 GM Score")+
  labs(color = "Subject")+
  #facet_wrap(~Visit,nrow=1, strip.position = "bottom") +
  
  theme(
    
    panel.spacing.x =unit(0, "lines")
    ,strip.background = element_blank()
    ,panel.border = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    #plot.background = element_rect(fill = "transparent", color = "black"), # bg of the plot
    #panel.grid.major = element_line(color = "grey"), # get rid of major grid
    #panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    
    axis.text.x = element_text(angle = 45,vjust = 1,hjust=1)
  )+
  geom_vline(xintercept=c(4.5,8.5,11.5,13.5), linetype="dashed", color = "black")

###export-----
ggsave("./output/20220802_GeneSig/timeVSgeneSig/n_2.5.png", width = 8, height = 2, dpi = "print", bg = "transparent")









##3.2 Appending Cmax_All-----

#Get Cmax form July data
July2022_Cmax_allVisit <- forNCA_July_2022_all %>%
  filter(TIME == 1&CYCLE == 1&(DAY == 1|DAY == 8|DAY == 15)) %>% #filer out C1D 1 & 8 & 15 where Cmax lies
  mutate(Dose = as.numeric(DOSE), 
         Cycle = CYCLE, 
         Day = DAY,
         .keep = "unused")%>% #align with gene signature data column name for easier combining
  
  bind_rows(Mar2022_Cmax_allVisit) #Making it cumulative 


AddedCmax_allVisit_reorder_geneSig0802 <- reorder_geneSig0802_2022%>%
  left_join(July2022_Cmax_allVisit,by = c("SUBJID", "Dose","Cycle","Day")) %>%
  # drop_na() this is not suitable since all new PK data introduces NA value to column PCORRES
  filter(!is.na(CONC))
#Using reorder_geneSig0802_2022 to retain full information
#Due to the PK data format change, July PK data was not cumulative, hence the bind_row operation in previous step


##3.3 prepare data for plot-----
plot_AddedCmax_allVisit_reorder_geneSig0802  <-AddedCmax_allVisit_reorder_geneSig0802%>%
  select(-PEMBRO.y)%>%
  mutate(PEMBRO = PEMBRO.x)%>%
  mutate(Cmax = CONC)%>%
  select(-TIME,-CONC)%>%
  group_by(PEMBRO,Dose,SUBJID,Cmax,Visit) %>%
  filter(FC.C1D1 == max(FC.C1D1)) %>% 
  mutate(MaxFC_C1D1 = FC.C1D1 )%>%  #clarify
  select(PEMBRO,Dose,SUBJID,Cmax,Visit,MaxFC_C1D1,Timepoint)









#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#_____________________----- 




#4. Plot-----







##4.2 Cmax_allVisit VS Max(FC_C1D1)_24hr-----
ggplot(data =filter(plot_AddedCmax_allVisit_reorder_geneSig0802, Timepoint == "24Hr"),
       aes(x=Cmax,
           y=MaxFC_C1D1)) +
  
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black")+
  geom_point(alpha = 0.8, size= 4, aes(color = factor(Dose),shape = PEMBRO))+
  
  scale_color_brewer(palette = "Reds") +
  labs(color = "Dose(mg)", shape = "PEMBRO", title ="Fold change of GS vs Cmax" )+
  #scale_y_log10(breaks = c(1,3,10,30,1.56), labels = c("1","3","10","30","1.56")) +
  scale_y_log10()+
  scale_x_log10()+
  xlab("Cmax (ng/mL)") +
  ylab("Fold Change of gene signature score")+
  #ggtitle("Max fold change compare to baseline")+
  #scale_size(range = c(4,7))+
  #scale_shape_manual(values = c(15,16,17))+
  
  theme(
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.line = element_line()
  ) +
  
  theme(
    panel.border = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_line(color = "grey"), # get rid of major grid
    # panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )+
  facet_grid(PEMBRO~factor(Visit, levels = c("C1D1","C1D8","C1D15" ,"C2D1" ,"C3D1","C4D1" )),labeller = labeller(Visit = label_value, PEMBRO = label_both))
#geom_hline(yintercept=1.56, linetype="dashed", color = "black")

ggsave("./output/2022_0802_GeneSig/24hr_MaxFC_vs_Cmax_all_wrapPEMBRO.png", width = 10, height = 4, dpi = "print", bg = "transparent")









###4.2.1 get slope and coeff-----------

plot_AddedCmax_allVisit_reorder_geneSig0802_coeff<- plot_AddedCmax_allVisit_reorder_geneSig0802 %>%
  filter(Timepoint == "24Hr")%>% 
  ungroup() %>%
  group_by(PEMBRO,Visit) %>%
  mutate(
    slope = round(lm(MaxFC_C1D1 ~Cmax)$coefficients[2],3),
    significance = summary(lm(MaxFC_C1D1 ~Cmax))$coefficients[2,4],
    n = n()
  ) %>% 
  slice(1) %>%
  ungroup()%>%
  arrange(PEMBRO,factor(Visit,levels = c("C1D1","C1D8","C1D15"))) %>%
  select(-Dose,-SUBJID,-Cmax,-MaxFC_C1D1)


write.csv(plot_AddedCmax_allVisit_reorder_geneSig0802_coeff,"output/20220802_GeneSig/coeff_summ.csv", row.names = FALSE)









##4.3. Mean FC vs Cmax-----

###4.3.1. Average GS-----

meanGS_0802 <- AddedCmax_allVisit_reorder_geneSig0802 %>%
  select(-PEMBRO.y)%>%
  mutate(PEMBRO = PEMBRO.x)%>%
  mutate(Cmax = CONC)%>%
  select(-TIME,-CONC) %>% #Trimming down unwanted columns
  filter(Timepoint %in% c("Pre","24Hr"))   #510 -> 274 records

#used for plot 4.3.2
summ_meanGS_0802 <- meanGS_0802  %>%
  group_by(SUBJID, Timepoint)%>%
  summarise( mean.GS = mean(STING24GMScore),
             mean.Conc = mean(Cmax),
             counts = n()) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("Pre","24Hr"))) 
summ_meanConc_0802 <- meanGS_0802 %>%
  group_by(SUBJID)%>%
  summarise( mean.conc = mean(Cmax)) #Create mean Cmax for each subject

#used for plot 4.3.3
wide_summ_meanGS_0802 <-summ_meanGS_0802 %>% #pivot to wider table for calculation of FC
  select(SUBJID,Timepoint,mean.GS) %>%
  pivot_wider(names_from = Timepoint, values_from = mean.GS) %>% #58007-008 does not have 24hr level
  mutate(mean.FC = `24Hr`/Pre) %>%
  left_join(select(summ_meanConc_0802, SUBJID, mean.conc), 
            by = "SUBJID") %>%
  left_join(select(meanGS_0802, SUBJID, PEMBRO, Dose),
            by = "SUBJID") %>% #Join cohort info
  distinct() %>% #removed duplicated row
  drop_na()



###4.3.2 GS VS Time-----
ggplot(data = summ_meanGS_0506) +
  geom_point(aes(x = Timepoint, y = mean.GS)) +
  geom_line(aes(x = Timepoint, y = mean.GS, group = SUBJID))






###4.3.3.FC VS Cmax-----
ggplot(data = wide_summ_meanGS_0802) +
  #N = 42, removing 1 with missing values
  
  
  geom_point(aes(x = mean.conc, y = mean.FC, shape = PEMBRO, color = factor(Dose)), size = 4)+
  geom_hline(yintercept = 1.56, linetype = "dashed")+
  geom_smooth(aes(x = mean.conc, y = mean.FC), method = "lm", se = FALSE,  color = "Black") +
  
  scale_color_brewer(palette = "Reds")+
  
  scale_x_log10()+
  scale_y_continuous(breaks = c(0,1.56,seq(4,12,by = 4)), labels = c("0","1.56","4","8","12"))+
  labs(color = "Dose (mg)",
       x = "Mean Cmax (ng/mL)",
       y = "Mean fold change (24hr VS baseline)")+
  
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

ggsave("./output/20220802_Genesig/meanFC_vs_meanCmax_faceted.png", width = 8, height = 6, dpi = "print", bg = "transparent")



lm.wide_summ_meanGS_0802 = lm(mean.FC~mean.conc, 
                              data = wide_summ_meanGS_0802) #linear regression
summary(lm.wide_summ_meanGS_0802)

# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.09040    0.59213   1.841 0.072973 .  
# mean.conc    0.11247    0.02683   4.192 0.000148 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




# Coefficients 0506_2022:
#                Estimate Std. Error t value Pr(>|t|)   
#   (Intercept)  1.20546    0.50056   2.408  0.02197 * 
#   mean.conc    0.07146    0.02449   2.917  0.00641 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


lm.wide_summ_meanGS_0802_N = lm(mean.FC~mean.conc, 
                              data = filter(wide_summ_meanGS_0802,PEMBRO == "N")) #linear regression
summary(lm.wide_summ_meanGS_0802_N)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)  1.57882    1.15234   1.370    0.187
# mean.conc    0.07902    0.04586   1.723    0.101
# 
# Residual standard error: 3.361 on 19 degrees of freedom
# Multiple R-squared:  0.1351,	Adjusted R-squared:  0.08961 
# F-statistic: 2.969 on 1 and 19 DF,  p-value: 0.1011




lm.wide_summ_meanGS_0802_Y = lm(mean.FC~mean.conc, 
                                data = filter(wide_summ_meanGS_0802,PEMBRO == "Y")) #linear regression
summary(lm.wide_summ_meanGS_0802_Y)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.32468    0.29983   1.083    0.292    
# mean.conc    0.17796    0.01619  10.994 1.12e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.8882 on 19 degrees of freedom
# Multiple R-squared:  0.8642,	Adjusted R-squared:  0.857 
# F-statistic: 120.9 on 1 and 19 DF,  p-value: 1.12e-09



###4.3.4.New data from MARCH-----

a <- setdiff(wide_summ_meanGS_0802$SUBJID,wide_summ_meanGS_0506$SUBJID)
# "07001-005" "07001-006" "07001-007" "58002-008" "58004-008" "58005-004" "58007-008" "58007-009"
#New subjects

wide_summ_meanGS_0802 <- wide_summ_meanGS_0802 %>%
  mutate(newDat = if_else(SUBJID %in% a, 1,0))



ggplot(data = wide_summ_meanGS_0802) +
  #N = 42, removing 1 with missing values
  
  geom_point(aes(x = mean.conc, y = mean.FC, shape = PEMBRO, color = factor(newDat, labels = c("March", "July"))),size = 4)+
  geom_hline(yintercept = 1.56, linetype = "dashed")+
  geom_smooth(aes(x = mean.conc, y = mean.FC), method = "lm", se = FALSE,  color = "Black") +
  
  scale_color_manual(values = c("black", "red"))+
  
  scale_x_log10()+
  scale_y_continuous(breaks = c(0,1.56,seq(4,12,by = 4)), labels = c("0","1.56","4","8","12"))+
  labs(color = "Date",
       x = "Mean Cmax (ng/mL)",
       y = "Mean fold change (24hr VS baseline)")+
  
  
  
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


ggsave("./output/20220802_Genesig/meanFC_vs_meanCmax_NewDat.png", width = 8, height = 6, dpi = "print", bg = "transparent")


newGS_subject_July_2022 <- 
  wide_summ_meanGS_0802%>%
  filter(SUBJID %in% a) %>%
  arrange(PEMBRO, desc(Dose), SUBJID)

write.csv(newGS_subject_July_2022, "output/20220802_GeneSig/new_subject.csv", row.names = FALSE)













##4.4.Add in STING1-----

###4.4.1 Combined-----
#Shared STING genotype from Zhenqiang,09162022, focusing on c.695A>G

STING_genotype <- read.csv("././Data/20220916_STING_WES/STING1.csv")


wide_summ_meanGS_0802_addedSTING695AG <- wide_summ_meanGS_0802 %>%
  left_join(STING_genotype,by = "SUBJID")%>%
  drop_na()
#N = 34


#Combined plot
ggplot(data = wide_summ_meanGS_0802_addedSTING695AG) +
  
  geom_point(aes(x = mean.conc, y = mean.FC, shape = PEMBRO, color = factor(c.695A.G_GeneCode, labels = c("High", "Low"))),size = 4)+
  geom_hline(yintercept = 1.56, linetype = "dashed")+
  geom_smooth(aes(x = mean.conc, y = mean.FC), method = "lm", se = FALSE,  color = "Black") +
  
  scale_color_manual(values = c("blue","black"))+
  
  scale_x_log10()+
  scale_y_continuous(breaks = c(0,1.56,seq(4,12,by = 4)), labels = c("0","1.56","4","8","12"))+
  labs(color = "Sting expression",
       x = "Mean Cmax (ng/mL)",
       y = "Mean fold change (24hr VS baseline)")+
  
  
  
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


ggsave("./output/20220916_STING_WES/meanFC_vs_meanCmax_STING695A.G.png", width = 8, height = 6, dpi = "print", bg = "transparent")


###4.4.2 By STING -----
#Seperate plot
ggplot(data = wide_summ_meanGS_0802_addedSTING695AG) +
  #N = 42, removing 1 with missing values
  
  
  geom_point(aes(x = mean.conc, y = mean.FC, shape = PEMBRO, color = factor(c.695A.G_GeneCode, labels = c("High", "Low"))), size = 4)+
  geom_hline(yintercept = 1.56, linetype = "dashed")+
  geom_smooth(aes(x = mean.conc, y = mean.FC), method = "lm", se = FALSE,  color = "Black") +
  
  scale_color_manual(values  = c("blue","black"))+
  
  scale_x_log10()+
  scale_y_continuous(breaks = c(0,1.56,seq(4,12,by = 4)), labels = c("0","1.56","4","8","12"))+
  labs(color = "STING expression",
       x = "Mean Cmax (ng/mL)",
       y = "Mean fold change (24hr VS baseline)")+
  
  facet_wrap(ncol = 1,vars(factor(c.695A.G_GeneCode, labels = c("High", "Low"))))+
  
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

ggsave("./output/20220916_STING_WES/meanFC_vs_meanCmax_STING695A.G_bySTING1.png", width = 8, height = 6, dpi = "print", bg = "transparent")

#Stats

#1.High expresser
lm.wide_summ_meanGS_0802_addedSTING695AG_High = lm(mean.FC~mean.conc, 
                                                   data = filter(wide_summ_meanGS_0802_addedSTING695AG,c.695A.G_GeneCode == 1)) #linear regression
summary(lm.wide_summ_meanGS_0802_addedSTING695AG_High)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)   1.2890     0.9418   1.369    0.201
# mean.conc     0.0689     0.0538   1.281    0.229
# 
# Residual standard error: 1.58 on 10 degrees of freedom
# Multiple R-squared:  0.1409,	Adjusted R-squared:  0.05499 
# F-statistic:  1.64 on 1 and 10 DF,  p-value: 0.2292




#2.Low expresser
lm.wide_summ_meanGS_0802_addedSTING695AG_Low = lm(mean.FC~mean.conc, 
                                                  data = filter(wide_summ_meanGS_0802_addedSTING695AG,c.695A.G_GeneCode == 2)) #linear regression
summary(lm.wide_summ_meanGS_0802_addedSTING695AG_Low)

# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)   
#   (Intercept)  1.74103    0.54329   3.205  0.00445 **
#   mean.conc    0.02492    0.02744   0.908  0.37471   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.889 on 20 degrees of freedom
# Multiple R-squared:  0.03959,	Adjusted R-squared:  -0.008433 
# F-statistic: 0.8244 on 1 and 20 DF,  p-value: 0.3747















###4.4.3 By other STING -----
#All STING polymorphism
STING_full <- read.csv("./Data/20220916_STING_WES/STING_full.csv")

STING_full_GeneCode <- STING_full %>%
  mutate(across(
    starts_with("rs"), ~ case_when
                        (.x =="0/0"~"0",
                         .x =="0/1"~"1",
                         .x =="1/1"~"2")),
    .keep = "all")


wide_summ_meanGS_0802_addedSTING_full <-wide_summ_meanGS_0802 %>%
  left_join(STING_full_GeneCode,by = "SUBJID")%>%
  drop_na()
#N = 34 
  
####4.4.3.1 rs7380824------
ggplot(data = wide_summ_meanGS_0802_addedSTING_full) +
  
  
  geom_point(aes(x = mean.conc, y = mean.FC, shape = PEMBRO, color = factor(rs7380824)), size = 4)+
  geom_hline(yintercept = 1.56, linetype = "dashed")+
  geom_smooth(aes(x = mean.conc, y = mean.FC), method = "lm", se = FALSE,  color = "Black") +
  
  scale_color_manual(values  = c("black","cyan","blue3"))+
  
  scale_x_log10()+
  scale_y_continuous(breaks = c(0,1.56,seq(4,12,by = 4)), labels = c("0","1.56","4","8","12"))+
  labs(color = "STING rs7380824 genotype",
       x = "Mean Cmax (ng/mL)",
       y = "Mean fold change (24hr VS baseline)")+
  
  facet_wrap(ncol = 1,vars(factor(PEMBRO, labels = c("Single agent", "Combination"))))+
  
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


ggsave("./output/20220916_STING_WES/otherSTING/meanFC_vs_meanCmax_rs7380824.png", width = 8, height = 6, dpi = "print", bg = "transparent")


####4.4.3.2 rs78233829------
ggplot(data = wide_summ_meanGS_0802_addedSTING_full) +
  
  
  geom_point(aes(x = mean.conc, y = mean.FC, shape = PEMBRO, color = factor(rs78233829)), size = 4)+
  geom_hline(yintercept = 1.56, linetype = "dashed")+
  geom_smooth(aes(x = mean.conc, y = mean.FC), method = "lm", se = FALSE,  color = "Black") +
  
  scale_color_manual(values  = c("black","cyan","blue3"))+
  
  scale_x_log10()+
  scale_y_continuous(breaks = c(0,1.56,seq(4,12,by = 4)), labels = c("0","1.56","4","8","12"))+
  labs(color = "STING rs78233829 genotype",
       x = "Mean Cmax (ng/mL)",
       y = "Mean fold change (24hr VS baseline)")+
  
  facet_wrap(ncol = 1,vars(factor(PEMBRO, labels = c("Single agent", "Combination"))))+
  
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

ggsave("./output/20220916_STING_WES/otherSTING/meanFC_vs_meanCmax_rs78233829.png", width = 8, height = 6, dpi = "print", bg = "transparent")









####4.4.3.3 rs7380272------
ggplot(data = wide_summ_meanGS_0802_addedSTING_full) +
  
  
  geom_point(aes(x = mean.conc, y = mean.FC, shape = PEMBRO, color = factor(rs7380272)), size = 4)+
  geom_hline(yintercept = 1.56, linetype = "dashed")+
  geom_smooth(aes(x = mean.conc, y = mean.FC), method = "lm", se = FALSE,  color = "Black") +
  
  scale_color_manual(values  = c("black","cyan","blue3"))+
  
  scale_x_log10()+
  scale_y_continuous(breaks = c(0,1.56,seq(4,12,by = 4)), labels = c("0","1.56","4","8","12"))+
  labs(color = "STING rs7380272 genotype",
       x = "Mean Cmax (ng/mL)",
       y = "Mean fold change (24hr VS baseline)")+
  
  facet_wrap(ncol = 1,vars(factor(PEMBRO, labels = c("Single agent", "Combination"))))+
  
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

ggsave("./output/20220916_STING_WES/otherSTING/meanFC_vs_meanCmax_rs7380272.png", width = 8, height = 6, dpi = "print", bg = "transparent")
















####4.4.3.3 rs151074578------
ggplot(data = wide_summ_meanGS_0802_addedSTING_full) +
  
  
  geom_point(aes(x = mean.conc, y = mean.FC, shape = PEMBRO, color = factor(	
    rs151074578)), size = 4)+
  geom_hline(yintercept = 1.56, linetype = "dashed")+
  geom_smooth(aes(x = mean.conc, y = mean.FC), method = "lm", se = FALSE,  color = "Black") +
  
  scale_color_manual(values  = c("black","cyan","blue3"))+
  
  scale_x_log10()+
  scale_y_continuous(breaks = c(0,1.56,seq(4,12,by = 4)), labels = c("0","1.56","4","8","12"))+
  labs(color = "STING rs151074578 genotype",
       x = "Mean Cmax (ng/mL)",
       y = "Mean fold change (24hr VS baseline)")+
  
  facet_wrap(ncol = 1,vars(factor(PEMBRO, labels = c("Single agent", "Combination"))))+
  
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

ggsave("./output/20220916_STING_WES/otherSTING/meanFC_vs_meanCmax_rs151074578.png", width = 8, height = 6, dpi = "print", bg = "transparent")


















####4.4.3.4 rs11554776------
ggplot(data = wide_summ_meanGS_0802_addedSTING_full) +
  
  
  geom_point(aes(x = mean.conc, y = mean.FC, shape = PEMBRO, color = factor(	
    rs11554776)), size = 4)+
  geom_hline(yintercept = 1.56, linetype = "dashed")+
  geom_smooth(aes(x = mean.conc, y = mean.FC), method = "lm", se = FALSE,  color = "Black") +
  
  scale_color_manual(values  = c("black","cyan","blue3"))+
  
  scale_x_log10()+
  scale_y_continuous(breaks = c(0,1.56,seq(4,12,by = 4)), labels = c("0","1.56","4","8","12"))+
  labs(color = "STING rs11554776 genotype",
       x = "Mean Cmax (ng/mL)",
       y = "Mean fold change (24hr VS baseline)")+
  
  facet_wrap(ncol = 1,vars(factor(PEMBRO, labels = c("Single agent", "Combination"))))+
  
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

ggsave("./output/20220916_STING_WES/otherSTING/meanFC_vs_meanCmax_rs11554776.png", width = 8, height = 6, dpi = "print", bg = "transparent")



















####4.4.3.4 rs7447927------
ggplot(data = wide_summ_meanGS_0802_addedSTING_full) +
  
  
  geom_point(aes(x = mean.conc, y = mean.FC, shape = PEMBRO, color = factor(	
    rs7447927)), size = 4)+
  geom_hline(yintercept = 1.56, linetype = "dashed")+
  geom_smooth(aes(x = mean.conc, y = mean.FC), method = "lm", se = FALSE,  color = "Black") +
  
  scale_color_manual(values  = c("black","cyan","blue3"))+
  
  scale_x_log10()+
  scale_y_continuous(breaks = c(0,1.56,seq(4,12,by = 4)), labels = c("0","1.56","4","8","12"))+
  labs(color = "STING rs7447927 genotype",
       x = "Mean Cmax (ng/mL)",
       y = "Mean fold change (24hr VS baseline)")+
  
  facet_wrap(ncol = 1,vars(factor(PEMBRO, labels = c("Single agent", "Combination"))))+
  
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

ggsave("./output/20220916_STING_WES/otherSTING/meanFC_vs_meanCmax_rs7447927.png", width = 8, height = 6, dpi = "print", bg = "transparent")







###4.4.2 By PEMBRO -----
#Seperate plot
ggplot(data = wide_summ_meanGS_0802_addedSTING695AG) +
  #N = 42, removing 1 with missing values
  
  
  geom_point(aes(x = mean.conc, y = mean.FC, shape = PEMBRO, color = factor(c.695A.G_GeneCode, labels = c("High", "Low"))), size = 4)+
  geom_hline(yintercept = 1.56, linetype = "dashed")+
  geom_smooth(aes(x = mean.conc, y = mean.FC), method = "lm", se = FALSE,  color = "Black") +
  
  scale_color_manual(values  = c("blue","black"))+
  
  scale_x_log10()+
  scale_y_continuous(breaks = c(0,1.56,seq(4,12,by = 4)), labels = c("0","1.56","4","8","12"))+
  labs(color = "STING expression",
       x = "Mean Cmax (ng/mL)",
       y = "Mean fold change (24hr VS baseline)")+
  
  facet_wrap(ncol = 1,vars(factor(PEMBRO, labels = c("Single agent", "Combination"))))+
  
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

ggsave("./output/20220916_STING_WES/meanFC_vs_meanCmax_STING695A.G_byPEMBRO.png", width = 8, height = 6, dpi = "print", bg = "transparent")




#CONCLUSION ------- 

#THIS gene signature score is shared on AUG 2, data curoff is 0729_2022
#3 subjects added, 2 in 2.0 COMBO arm and 1 in 2.5 SA.