# code to explore predictors of sample size

# LIBRARIES REQUIRED
library(ggplot2)
library(stats)

# IMPORTING DATA & FORMATTING
db <- read.csv("data/NGSscoping-db-final.csv", header = TRUE)
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
hist(log10(db$total_samples))

names(db)
glmP1 <- glm(total_samples ~ study_aim_1+ total_domains + agent_type, data = db, family = poisson(link = "log"))
summary(glmP1)

glmNB <- glm.nb(total_samples ~ study_aim_1+ total_domains + agent_type, data = db, link = "log")
summary(glmNB)

glmQP1 <- glm(adj_samples ~ study_aim_1+ total_domains + agent_type, data = db, family = quasipoisson(link = "log"))
summary(glmQP1)

with(glmP1, cbind(res.deviance = deviance, df = df.residual,
               p = pchisq(deviance, df.residual, lower.tail=FALSE)))


# negative bionomial models
glmNB1 <- glm.nb(total_samples ~ env_bin + ani_bin + study_aim_1, data = db, link = "log")
summary(glmNB1)
plot(glmNB1)
# residuals look terrible

glmNB2 <- glm.nb(total_samples ~ study_aim_1, data = db, link = "log")
summary(glmNB2)

glmNB3 <- glm.nb(total_samples ~ total_domains, data = db, link = "log")
summary(glmNB3)

glmNB4 <- glm.nb(total_samples ~ agent_hse_cat, data = db, link = "log")
summary(glmNB4)

glmNB5 <- glm.nb(total_samples ~ phylo_results, data = db, link = "log")
summary(glmNB5)
unique(db$phylo_results)
# interesting that temporally resolved isn't significantly higher than genetic distace, but relatedness is

glmNB6 <- glm.nb(total_samples ~ agent_genus, data = db, link = "log")
summary(glmNB6)
# some significant, but low df

glmNB7 <- glm.nb(total_samples ~ ngs_platform_primary, data = db, link = "log")
summary(glmNB7)
# no significance

glmNB8 <- glm.nb(total_samples ~ ngs_geo, data = db, link = "log")
summary(glmNB8)

glmNB9 <- glm.nb(total_samples ~ publication_year, data = db, link = "log")
summary(glmNB9)
plot(glmNB9)
# 74 is outlier

glmNB10 <- glm.nb(total_samples ~ publication_year + umic_hic, data = db, link = "log")
summary(glmNB10)
plot(glmNB10)

glmNB11 <- glm.nb(total_samples ~ publication_year + lic_lmic, data = db, link = "log")
summary(glmNB11)

names(db)
glmNB12 <- glm.nb(total_samples ~ publication_year + agent_type, data = db, link = "log")
summary(glmNB12)

glmNB13 <- glm.nb(total_samples ~ agent_sapo, data = db, link = "log")
summary(glmNB13)

glmNB14 <- glm.nb(total_samples ~ agent_type, data = db, link = "log")
summary(glmNB14)

glmNB15 <- glm.nb(total_samples ~ log10(hum_sample_size), data = db, link = "log")
summary(glmNB15)

glmNB16 <- glm.nb(total_samples ~ publication_year + agent_type, data = db, link = "log")
summary(glmNB16)


# removing humans

db$adj_samples = db$total_samples-db$hum_sample_size


# negative bionomial models
glmNB1 <- glm.nb(adj_samples ~ env_bin + ani_bin + study_aim_1, data = db, link = "log")
summary(glmNB1)
plot(glmNB1)
# residuals look terrible

glmNB2 <- glm.nb(adj_samples ~ study_aim_1, data = db, link = "log")
summary(glmNB2)

glmNB3 <- glm.nb(adj_samples ~ total_domains, data = db, link = "log")
summary(glmNB3)

glmNB4 <- glm.nb(adj_samples ~ agent_hse_cat, data = db, link = "log")
summary(glmNB4)

glmNB5 <- glm.nb(adj_samples ~ phylo_results, data = db, link = "log")
summary(glmNB5)
unique(db$phylo_results)
# interesting that temporally resolved isn't significantly higher than genetic distace, but relatedness is

glmNB6 <- glm.nb(adj_samples ~ agent_genus, data = db, link = "log")
summary(glmNB6)
# some significant, but low df

glmNB7 <- glm.nb(adj_samples ~ ngs_platform_primary, data = db, link = "log")
summary(glmNB7)
# no significance

glmNB8 <- glm.nb(adj_samples ~ ngs_geo, data = db, link = "log")
summary(glmNB8)

glmNB9 <- glm.nb(adj_samples ~ publication_year, data = db, link = "log")
summary(glmNB9)
plot(glmNB9)
# more outliers

glmNB10 <- glm.nb(adj_samples ~ publication_year + umic_hic, data = db, link = "log")
summary(glmNB10)

glmNB11 <- glm.nb(adj_samples ~ publication_year + lic_lmic, data = db, link = "log")
summary(glmNB11)


glmNB12 <- glm.nb(adj_samples ~ publication_year + agent_type, data = db, link = "log")
summary(glmNB12)
plot(glmNB12)

glmNB13 <- glm.nb(adj_samples ~ agent_sapo, data = db, link = "log")
summary(glmNB13)

glmNB14 <- glm.nb(adj_samples ~ agent_type, data = db, link = "log")
summary(glmNB14)

glmNB15 <- glm.nb(adj_samples ~ log10(hum_sample_size), data = db, link = "log")
summary(glmNB15)

glmNB16 <- glm.nb(adj_samples ~ publication_year + agent_type + log10(hum_sample_size), data = db, link = "log")
summary(glmNB16)
