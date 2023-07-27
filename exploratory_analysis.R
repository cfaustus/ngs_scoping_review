# code to explore and visualize dataset

# libraries
library(ggplot2)

# importing and correctly formatting data
db4 <- read.csv("data/NGSscoping_dbv4.csv", header = TRUE)
names(db4)<-tolower(names(db4))
db4$biotic.env_sample_size
colnames(db4)[22] <- 'biotic_env_sample_size'

#'env_domain' looks like a sample size.... 

str(db4) # lots of columns coming up as characters
col.chr.ids = c('hum_sample_size','livestock_sample_size', 'pet_sample_size','poultry_sample_size',
                'wildlife_sample_size','arachnida_sample_size', 'insecta_sample_size', 'env_domain',
                'biotic_env_sample_size','abiotic_sample_size')
db4[col.chr.ids] <- sapply(db4[col.chr.ids],as.numeric)
str(db4)

sample_cols = grep('sample_size',names(db4)) # identifying sample size columns
db4$total_samples = rowSums(db4[,sample_cols], na.rm=TRUE)
db4$env_bin = ifelse(db4$env_domain>0,1,0)
db4$ani_bin = ifelse(db4$ani_sample_size>0,1,0)
db4$human_bin = ifelse(db4$hum_sample_size>0,1,0)
db4$total_domains = rowSums(db4[,c('env_bin','ani_bin','human_bin')], na.rm=T)
# some are 0?!? and some are 1?

# by publication year
sort(unique(db4$publication_year))
library(ggplot2)
ggplot(db4, aes(x=publication_year, fill=agent_type)) +
  geom_histogram(col = 'white', position="stack", bins = 12)+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="year published", y = "publications")+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,23,2))+
  scale_x_continuous(breaks=seq(2011,2022,1))

ggplot(db4, aes(x=publication_year, fill=ngs_platform_short_primary)) +
  geom_histogram(col = 'white', position="stack", bins = 12)+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="year published", y = "publications")+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,23,2))+
  scale_x_continuous(breaks=seq(2011,2022,1))

ggplot(db4, aes(x=publication_year, fill=study_aim_1)) +
  geom_histogram(col = 'white', position="stack", bins = 12)+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="year published", y = "publications")+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,23,2))+
  scale_x_continuous(breaks=seq(2011,2022,1))

# by sample size year
names(db4)
ggplot(db4, aes(x=total_samples, y = total_domains, col=study_aim_1)) +
  geom_point()+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="total samples", y = "domains")+
  scale_x_continuous(trans='log10')

ggplot(db4, aes(x=total_samples, y = total_domains, col=agent_type)) +
  geom_point()+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="total samples", y = "domains")+
  scale_x_continuous(trans='log10')
