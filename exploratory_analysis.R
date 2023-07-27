# code to explore and visualize dataset

# libraries
library(ggplot2)

# importing and correctly formatting data
db4 <- read.csv("data/NGSscoping_dbv4.csv", header = TRUE)
names(db4)<-tolower(names(db4))
db4$biotic.env_sample_size
colnames(db4)[22] <- 'biotic_env_sample_size'

str(db4) # lots of columns coming up as characters
col.chr.ids = c('hum_sample_size','livestock_sample_size', 'pet_sample_size','poultry_sample_size',
                'wildlife_sample_size','env_domain',
db4[col.chr.ids] <- sapply(db4[col.chr.ids],as.numeric)


db4$hum_sample_size = as.numeric(db4$hum_sample_size)
sample_cols = grep('sample_size',names(db4)) # identifying sample size columns
rowSums(db4[,sample_cols], na.rm=TRUE)



db4 %>% mutate_if(samples, as.double)  
View(db4)
head(db4)

hist(log10(db4$ANI_SAMPLE_SIZE+1))
hist(log10(as.numeric(db4$HUM_SAMPLE_SIZE)))

# need to fix column name importations .. not showing up as correct values
sort(unique(db4$PUBLICATION_YEAR))
library(ggplot2)
ggplot(db4, aes(x=PUBLICATION_YEAR, fill=AGENT_TYPE)) +
  geom_histogram(col = 'white', position="stack", bins = 12)+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="year published", y = "publications")+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,23,2))+
  scale_x_continuous(breaks=seq(2011,2022,1))
  
sort(unique())
names(db4)
ggplot(db4, aes(x=PUBLICATION_YEAR, fill=NGS_PLATFORM_SHORT_PRIMARY)) +
  geom_histogram(col = 'white', position="stack", bins = 12)+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="year published", y = "publications")+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,23,2))+
  scale_x_continuous(breaks=seq(2011,2022,1))

ggplot(db4, aes(x=PUBLICATION_YEAR, fill=STUDY_AIM_1)) +
  geom_histogram(col = 'white', position="stack", bins = 12)+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="year published", y = "publications")+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,23,2))+
  scale_x_continuous(breaks=seq(2011,2022,1))

ggplot(db4, aes(x=PUBLICATION_YEAR, fill=STUDY_AIM_1)) +
  geom_histogram(col = 'white', position="stack", bins = 12)+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="year published", y = "publications")+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,23,2))+
  scale_x_continuous(breaks=seq(2011,2022,1))

sample_sizes =  sum(db4$HUM_SAMPLE_SIZE,
                db4$ANI_SAMPLE_SIZE, 
                db4$LIVESTOCK_SAMPLE_SIZE,
                db4$POULTRY_SAMPLE_SIZE,
                db4$WILDLIFE_SAMPLE_SIZE, 
                db4$ARACHNIDA_SAMPLE_SIZE,
                db4$INSECTA_SAMPLE_SIZE,
                db4$`BIOTIC-ENV_SAMPLE_SIZE`, 
                db4$ABIOTIC_SAMPLE_SIZE,
                na.rm = T)
  
[13] "ANI_MULTISPECIES"           "ANI_SAMPLE_SIZE"           
[15] "LIVESTOCK_SAMPLE_SIZE"      "PET_SAMPLE_SIZE"           
[17] "POULTRY_SAMPLE_SIZE"        "WILDLIFE_SAMPLE_SIZE"      
[19] "ENV_DOMAIN"                 "ARACHNIDA_SAMPLE_SIZE"     
[21] "INSECTA_SAMPLE_SIZE"        "BIOTIC-ENV_SAMPLE_SIZE"    
[23] "ABIOTIC_SAMPLE_SIZE"  



