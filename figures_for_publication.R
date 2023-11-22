
# code to explore and visualize dataset

# libraries
library(ggplot2)
library(ggalluvial)
library(dplyr)

# importing and correctly formatting data
db <- read.csv("data/NGSscoping_dbv6.csv", header = TRUE)
names(db)<-tolower(names(db))
db$biotic.env_sample_size
colnames(db)[29] <- 'biotic_env_sample_size'

#'env_domain' looks like a sample size.... 

str(db) # lots of columns coming up as characters
col.chr.ids = c('hum_sample_size','livestock_sample_size', 'pet_sample_size','poultry_sample_size',
                'wildlife_sample_size','arachnida_sample_size', 'insecta_sample_size', 'env_domain',
                'biotic_env_sample_size','abiotic_sample_size')
db[col.chr.ids] <- sapply(db[col.chr.ids],as.numeric)
str(db)

#sample_cols = grep('sample_size',names(db)) # identifying sample size columns - WRONG double counting
db$total_samples = rowSums(db[,c('hum_sample_size','ani_sample_size', 'env_domain')], na.rm=TRUE)
db$env_bin = ifelse(db$env_domain>0,1,0)
db$ani_bin = ifelse(db$ani_sample_size>0,1,0)
db$human_bin = ifelse(db$hum_sample_size>0,1,0)
db$total_domains = rowSums(db[,c('env_bin','ani_bin','human_bin')], na.rm=T)

# some are 0?!? and some are 1?
# annotation in excel clarified there were 276 samples from human and environment but split not given

ind = which(db$total_samples == 0)
db$total_domains[ind]<-2
db$total_samples[ind]<-276
db$env_bin[ind] = 1
db$human_bin[ind] = 1




# by publication year
sort(unique(db$publication_year))
library(ggplot2)
ggplot(db, aes(x=publication_year, fill=agent_type)) +
  geom_histogram(col = 'white', position="stack", bins = 12)+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="year published", y = "publications")+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,23,2))+
  scale_x_continuous(breaks=seq(2011,2022,1))
db$ngs_platform_primary
ggplot(db, aes(x=publication_year, fill=ngs_platform_primary)) +
  geom_histogram(col = 'white', position="stack", bins = 12)+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="year published", y = "publications")+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,23,2))+
  scale_x_continuous(breaks=seq(2011,2022,1))

ggplot(db, aes(x=publication_year, fill=study_aim_1)) +
  geom_histogram(col = 'white', position="stack", bins = 12)+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="year published", y = "publications")+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,23,2))+
  scale_x_continuous(breaks=seq(2011,2022,1))

# by sample size year
names(db)
ggplot(db, aes(x=total_samples, y = total_domains, col=study_aim_1)) +
  geom_point()+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="total samples", y = "domains")+
  scale_x_continuous(trans='log10')

ggplot(db, aes(x=total_samples, y = total_domains, col=agent_type)) +
  geom_point()+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="total samples", y = "domains")+
  scale_x_continuous(trans='log10')

# https://cran.r-project.org/web/packages/ggalluvial/vignettes/ggalluvial.html
head(db)
db_allu = is_alluvia_form(as.data.frame(db), axes = 1:3, silent = TRUE)
ggplot(data = db,
       aes(x = publication_year, y = total_samples, alluvium = study_aim_1)) +
  geom_alluvium(aes(fill = study_aim_1, colour = study_aim_1),
                alpha = .75, decreasing = FALSE) 
scale_x_continuous(breaks = seq(2003, 2013, 2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -30, hjust = 0)) +
  scale_fill_brewer(type = "qual", palette = "Set3") +
  scale_color_brewer(type = "qual", palette = "Set3") +
  facet_wrap(~ region, scales = "fixed") +
  ggtitle("refugee volume by country and region of origin")


db$ngs_platform_short_primary
head(db)
db_summary = as.data.frame(db %>% group_by(ngs_platform_primary, publication_year) %>%
                             summarise(papers = n()))
db_summary$ngs_platform_primary=as.factor(db_summary$ngs_platform_primary)

db_summary = db_summary[1:45,]
db_summary$ngs_platform_short_primary
unique(db$ngs_platform_primary)

db_summary$ngs_platform_primary = as.factor(db_summary$ngs_platform_primary)

levels(db_summary$ngs_platform_primary)
db_summary$ngs_platform_primary <- factor(db_summary$ngs_platform_primary, 
                                          levels = c("BGI_MGISEQ2000", "Illumina_iSeq100","Illumina_NovaSeq",
                                                     "BGI_BGISEQ50", 
                                                     "BGI_BGISEQ500", "ONT_MinION","Illumina_NextSeq", 
                                                     "Illumina_MiSeq", "IonTorrent_IonPGM", 
                                                     "Illumina_HiSeq", "Illumina_GAIIx", 
                                                     "Roche_454"))

ggplot(db_summary, aes(x = publication_year, y = ngs_platform_short_primary )) +
  geom_point(aes(size = papers, color = ngs_platform_company))+ # aes(dotsize = papers),
  scale_x_continuous(breaks=seq(2011,2022,1))+
  theme_classic()+
  labs(x="publication year", y = "sequencer")+
  scale_size_continuous(breaks = c(min(db_summary$papers),
                                   3,
                                   max(db_summary$papers)),
                        labels = c("1",
                                   "3",
                                   "9"))+
  geom_rect(aes(xmin = 2011, xmax = 2022, ymin = "Illumina_HiSeq", ymax = "Illumina_HiSeq"))+
  annotate("rect", x = 2011, y = 2022, ymin = 12, ymax = 28)
# need to combine short read and long read sequencing

# is within host diversity studies predicted by mutation rate?
# maybe not relevant b/c no host diversity

# alluvial plot to make
# open data 
# open software
# open code
# colored by country?
db$public_data_binary = !db$public_data == 'N/A'
db$public_code_binary = db$public_code == 'yes'
db$phylo_software_binary = !db$phylo_software == 'Geneious'
db$sampling_geo

db_summary2 = as.data.frame(db %>% group_by(public_data_binary, public_code_binary, phylo_software_binary,sampling_geo) %>%
                              summarise(Freq = n()))

ggplot(data = db_summary2,
       aes(axis1 = as.factor(phylo_software_binary), axis2 = public_data_binary,axis3 = public_code_binary, 
           y = Freq)) +
  scale_x_discrete(limits = c("software", "data", "code"), expand = c(.2, .05)) +
  xlab("open science") +
  geom_alluvium(aes(fill = sampling_geo)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() 
# need to update with stefano's linkages
# figure - sample size by HG status?
# software updates? license 

db5 = read.csv('data/NGSscoping_dbv5.csv')
names(db5)<-tolower(names(db5))
db5$geo_scope = ifelse(db5$lic_lmic == 'yes' & db5$umic_hic == 'yes', 'both',
                       ifelse(db5$lic_lmic == 'yes', 'lmic', 'hic'))

db5$public_code
db5$public_data_binary = !db5$public_data == 'N/A'
db5$public_code_binary = db5$public_code == 'yes'
db5$phylo_software_binary = db5$assembly.mapping_license == 'free'
db5_summary2 = as.data.frame(db5 %>% group_by(public_data_binary, public_code_binary, phylo_software_binary,geo_scope) %>%
                               summarise(Freq = n()))

db5_summary2$geo_scope = as.factor(db5_summary2$geo_scope)
db5_summary2$geo_scope = factor(db5_summary2$geo_scope, levels = c('lmic','both','hic'))
ggplot(data = db5_summary2,
       aes(axis1 = phylo_software_binary, axis2 = public_data_binary,axis3 = public_code_binary, 
           y = Freq)) +
  scale_x_discrete(limits = c("software", "data", "code"), expand = c(.2, .05)) +
  xlab("open science") +
  geom_alluvium(aes(fill = geo_scope)) +
  scale_fill_manual(values = c('#0c2c84','#7fcdbb', '#edf8b1'))+
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_classic() 




names(db5)
db5$time_to_pub = db5$publication_year-as.numeric(db5$last_sample_date)
ggplot(db5, aes(x=time_to_pub, fill=study_aim_1)) +
  geom_density(alpha = 0.4)+
  theme(legend.position="top")+
  # facet_wrap(~study_aim_1, nrows = 5) +
  theme_classic()+
  labs(x="years from last collection to publication", y = "papers")


