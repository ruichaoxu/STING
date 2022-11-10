library(tidyverse)
library(ggplot2)
library(magrittr)
library(stringr)

#1.0 Data Processing-----

lab.SDTM.source <- read.csv("./Data/20221108_labVisualizer/LB_TAK-676-1002_Dynamic Listings_25-Oct-2022.csv")

lab.SDTM.source.xaxis <- lab.SDTM.source %>%
  mutate(
         xaxis.num = if_else(is.na(Planned.Time.Point.Number..LBTPTNUM.),
                             Visit.Number..VISITNUM.,
                             Visit.Number..VISITNUM.+0.0001*Planned.Time.Point.Number..LBTPTNUM.),
         xaxis.description = str_c(Visit.Name..VISIT.,Planned.Time.Point.Name..LBTPT., sep="_"),
         PEMBRO = case_when(
           str_detect(Description.of.Planned.Arm..ARM.,"\\-(.+)\\+") == 1 ~ "1. Single Agent", #assign PEMBRO
           str_detect(Description.of.Planned.Arm..ARM.,"\\+") == 1 ~ "2. COMBO",
           TRUE ~ "1. Single Agent"),
         DOSE = str_extract(Description.of.Planned.Arm..ARM., #extract dose
                            "([:digit:]\\.[:digit:])"),
         RESULT.C1D1 = if_else(xaxis.num == 1.0098, 	
                          Numeric.Result.Finding.in.Standard.Units..LBSTRESN.,
                          NA_real_)
         )%>%
  group_by(Unique.Subject.Identifier..USUBJID.,Lab.Test.or.Examination.Short.Name..LBTESTCD.)%>%
  fill(RESULT.C1D1, .direction ="down")%>%
  mutate(RESULT.FC.C1D1 = Numeric.Result.Finding.in.Standard.Units..LBSTRESN./RESULT.C1D1)
#add unique x axis coordinate






#2.0 Plotting-----

#data filter
#for line
plot.data<- filter(lab.SDTM.source.xaxis,
                    Lab.Test.or.Examination.Short.Name..LBTESTCD. == "FCT54899"
                    &Visit.Number..VISITNUM. <3 
                    &!is.na(xaxis.num)
                   )
#For dot indicating Predose level
plot.data.2<- filter(lab.SDTM.source.xaxis,
                     Lab.Test.or.Examination.Short.Name..LBTESTCD. == "FCT54899"
                     &xaxis.num %in% c(1.0098,1.0798,1.1498,2.0098) 
                     &!is.na(xaxis.num)
)
#FCT54899:3+8+4-Ki67+(%8);  FCT54898:3+4+8-Ki67+(%4)

#Auto generate x-axis label based on visit and time point
xaxis.vector.1<-plot.data%>%
  group_by(xaxis.num) %>%
  slice(1) %>%
  pull(xaxis.description)

xaxis.vector.2 <- plot.data.2 %>%
  group_by (xaxis.num)%>%
  slice(1)%>%
  pull(xaxis.description)

#Plot
ggplot(data = plot.data) +
  
  geom_line(aes(x = factor(xaxis.num, labels=xaxis.vector.1), 
                y = Numeric.Result.Finding.in.Standard.Units..LBSTRESN., 
                group = Unique.Subject.Identifier..USUBJID., 
                color = Unique.Subject.Identifier..USUBJID.)) +
  geom_point(
    data = plot.data.2,
    aes(x = factor(xaxis.num, labels=xaxis.vector.2), 
                 y = Numeric.Result.Finding.in.Standard.Units..LBSTRESN., 
                 group = Unique.Subject.Identifier..USUBJID., 
                 color = Unique.Subject.Identifier..USUBJID.))+
 
  facet_grid(PEMBRO~DOSE)+
  labs(y = paste0(
    unique(	plot.data$Lab.Test.or.Examination.Name..LBTEST.),
    " (",
    unique(	plot.data$Original.Units..LBORRESU.),
    ")"),
    x = "Time")+
  theme(axis.text.x = element_text(angle = 330, hjust = 0))

