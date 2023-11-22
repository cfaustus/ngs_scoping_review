# code to explore predictors of sample size

# LIBRARIES REQUIRED
library(ggplot2)
library(stats)
library(MASS)
library(dplyr, warn.conflicts= FALSE)
library(GlmSimulatoR)


# IMPORTING DATA & FORMATTING
db <- read.csv("data/NGSscoping_dbv6.csv", header = TRUE)
names(db)<-tolower(names(db))

colnames(db)[29] <- 'biotic_env_sample_size' # changing names to 

str(db) # a lot of numeric columns are read in as characters
col.chr.ids = c('hum_sample_size','livestock_sample_size', 'pet_sample_size','poultry_sample_size',
                'wildlife_sample_size','arachnida_sample_size', 'insecta_sample_size', 'env_domain',
                'biotic_env_sample_size','abiotic_sample_size')
db[col.chr.ids] <- sapply(db[col.chr.ids],as.numeric)
str(db)


# TOTAL SAMPLE SIZES
db$total_samples = rowSums(db[,c('hum_sample_size','ani_sample_size', 'env_domain')], na.rm=TRUE)
db$env_bin = ifelse(db$env_domain>0,1,0)
db$ani_bin = ifelse(db$ani_sample_size>0,1,0)
db$human_bin = ifelse(db$hum_sample_size>0,1,0)
db$total_domains = rowSums(db[,c('env_bin','ani_bin','human_bin')], na.rm=T)

# one row is 0
ind = which(db$total_samples == 0) # row that neeeds manual correction
# annotation in excel clarified there were 276 samples from human and environment but specfic split not given
db$total_domains[ind]<-2
db$total_samples[ind]<-276
db$env_bin[ind] = 1
db$human_bin[ind] = 1

# visualizing database
names(db)
ggplot(db, aes(x=total_samples, fill=agent_type)) +
  geom_histogram()+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="total samples", y = "total")+
  scale_x_continuous(trans='log10')

ggplot(db, aes(x=total_samples, fill=as.factor(total_domains))) +
  geom_density(alpha = 0.5)+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="total samples", y = "density")+
  scale_x_continuous(trans='log10')


# model construction

glmP1 <- glm(total_samples ~ 1, data = db5, family = poisson(link = "log"))
summary(glmP1)

glmNB <- glm.nb(total_samples ~ 1, data = db5, link = "log")
summary(glmNB)

names(db5)
glmNB <- glm.nb(total_samples ~ env_bin + ani_bin + agent_hse_cat + study_aim_1, data = db5, link = "log")
summary(glmNB)

glmNB2 <- glm.nb(total_samples ~ study_aim_1, data = db5, link = "log")
summary(glmNB2)

glmNB3 <- glm.nb(total_samples ~ ani_bin + env_bin, data = db5, link = "log")
summary(glmNB3)

glmNB4 <- glm.nb(total_samples ~ agent_hse_cat, data = db5, link = "log")
summary(glmNB4)

glmNB4 <- glm.nb(total_samples ~ phylo_results, data = db5, link = "log")
summary(glmNB4)
