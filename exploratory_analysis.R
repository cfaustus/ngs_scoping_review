library(readr)
ngsscop_dbv4 <- read_csv("data/NGSscoping_dbv4.csv")
spec(ngsscop_dbv4)
samples = grepl('SAMPLE_SIZE',names(ngsscop_dbv4))
ngsscop_dbv4 %>% mutate_if(samples, as.double)  
View(ngsscop_dbv4)
head(ngsscop_dbv4)

hist(log10(ngsscop_dbv4$ANI_SAMPLE_SIZE+1))
hist(log10(as.numeric(ngsscop_dbv4$HUM_SAMPLE_SIZE)))

# need to fix column name importations .. not showing up as correct values
sort(unique(ngsscop_dbv4$PUBLICATION_YEAR))
library(ggplot2)
ggplot(ngsscop_dbv4, aes(x=PUBLICATION_YEAR, fill=AGENT_TYPE)) +
  geom_histogram(col = 'white', position="stack", bins = 12)+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="year published", y = "publications")+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,23,2))+
  scale_x_continuous(breaks=seq(2011,2022,1))
  
sort(unique())
names(ngsscop_dbv4)
ggplot(ngsscop_dbv4, aes(x=PUBLICATION_YEAR, fill=NGS_PLATFORM_SHORT_PRIMARY)) +
  geom_histogram(col = 'white', position="stack", bins = 12)+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="year published", y = "publications")+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,23,2))+
  scale_x_continuous(breaks=seq(2011,2022,1))

ggplot(ngsscop_dbv4, aes(x=PUBLICATION_YEAR, fill=STUDY_AIM_1)) +
  geom_histogram(col = 'white', position="stack", bins = 12)+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="year published", y = "publications")+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,23,2))+
  scale_x_continuous(breaks=seq(2011,2022,1))

ggplot(ngsscop_dbv4, aes(x=PUBLICATION_YEAR, fill=STUDY_AIM_1)) +
  geom_histogram(col = 'white', position="stack", bins = 12)+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="year published", y = "publications")+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,23,2))+
  scale_x_continuous(breaks=seq(2011,2022,1))

sample_sizes =  sum(ngsscop_dbv4$HUM_SAMPLE_SIZE,
                ngsscop_dbv4$ANI_SAMPLE_SIZE, 
                ngsscop_dbv4$LIVESTOCK_SAMPLE_SIZE,
                ngsscop_dbv4$POULTRY_SAMPLE_SIZE,
                ngsscop_dbv4$WILDLIFE_SAMPLE_SIZE, 
                ngsscop_dbv4$ARACHNIDA_SAMPLE_SIZE,
                ngsscop_dbv4$INSECTA_SAMPLE_SIZE,
                ngsscop_dbv4$`BIOTIC-ENV_SAMPLE_SIZE`, 
                ngsscop_dbv4$ABIOTIC_SAMPLE_SIZE,
                na.rm = T)
  
[13] "ANI_MULTISPECIES"           "ANI_SAMPLE_SIZE"           
[15] "LIVESTOCK_SAMPLE_SIZE"      "PET_SAMPLE_SIZE"           
[17] "POULTRY_SAMPLE_SIZE"        "WILDLIFE_SAMPLE_SIZE"      
[19] "ENV_DOMAIN"                 "ARACHNIDA_SAMPLE_SIZE"     
[21] "INSECTA_SAMPLE_SIZE"        "BIOTIC-ENV_SAMPLE_SIZE"    
[23] "ABIOTIC_SAMPLE_SIZE"  



