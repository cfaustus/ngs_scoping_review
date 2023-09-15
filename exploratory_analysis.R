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
db$agent_hse_cat
ggplot(db, aes(x=total_samples, y = total_domains, col= agent_hse_cat)) +
  geom_point()+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="total samples", y = "domains")+
  scale_x_continuous(trans='log10')

# https://cran.r-project.org/web/packages/ggalluvial/vignettes/ggalluvial.html
head(db)

  


db_summary = as.data.frame(db %>% group_by(ngs_platform_primary, publication_year) %>%
                                summarise(papers = n()))
db_summary = db_summary[!db_summary$ngs_platform_primary == 'Illumina_N/S',]
db_summary = db_summary[!db_summary$ngs_platform_primary == 'N/S',]

db_summary_v2 = as.data.frame(db %>% group_by(ngs_platform_complement, publication_year) %>%
                             summarise(papers_sec = n()))
db_summary_v2 = db_summary_v2[!db_summary_v2$ngs_platform_complement == 'N/S',]
db_summary_v2 = db_summary_v2[!db_summary_v2$ngs_platform_complement == 'no',]

db_platform <- merge(x=db_summary,y=db_summary_v2, by.x=c("ngs_platform_primary","publication_year"), 
             by.y=c("ngs_platform_complement","publication_year"), all = TRUE)
db_platform$tot_papers = rowSums(db_platform[,c('papers_sec','papers')], na.rm = TRUE)

db_platform$ngs_platform_primary=as.factor(db_platform$ngs_platform_primary)

levels(db_platform$ngs_platform_primary)
db_platform$ngs_platform_primary <- factor(db_platform$ngs_platform_primary, 
                levels = c("BGI_MGISEQ2000", "Illumina_iSeq100","Illumina_NovaSeq",
                           "BGI_BGISEQ50", 
                           "BGI_BGISEQ500", "ONT_MinION","Illumina_NextSeq", "PacBio_RSII",
                           "Illumina_MiSeq", "IonTorrent_IonPGM", 
                           "Illumina_HiSeq", "Illumina_GAIIx", 
                            "Roche_454"))

seq = read.csv('data/sequencing_tech_v2.csv', header = TRUE)
names(seq)
seq_trun = seq[,c('company', 'sequence')]
db_platform <- merge(x=db_platform,y=seq_trun, by.x=c("ngs_platform_primary"), 
                     by.y=c("sequence"), all = TRUE)
seq3 = read.csv('data/sequencing_tech_v3.csv', header = TRUE)

ggplot(db_platform, aes(x = publication_year, y = ngs_platform_primary )) +
  geom_point(aes(size = tot_papers, color = company))+ # aes(dotsize = papers),
  scale_x_continuous(limits = c(2006,2024), breaks=seq(2006,2024,1))+
  theme_classic()+
  labs(x="year", y = "sequencer")+
  scale_size_continuous(breaks = c(min(db_platform$tot_papers),
                                    3,
                                    max(db_platform$tot_papers)),
                         labels = c("1",
                                    "3",
                                    "9"))+
  geom_line(data = seq3, aes(y = ngs_platform_primary,x = publication_year, color = company))
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


db$geo_scope = ifelse(db$lic_lmic == 'yes' & db$umic_hic == 'yes', 'both',
                       ifelse(db$lic_lmic == 'yes', 'lmic', 'hic'))

db$public_code
db$public_data_binary = !db$public_data == 'N/A'
db$public_code_binary = db$public_code == 'yes'
db$phylo_software_binary = db$assembly.mapping_license == 'free'
db_summary2 = as.data.frame(db %>% group_by(public_data_binary, public_code_binary, phylo_software_binary,geo_scope) %>%
                               summarise(Freq = n()))

db_summary2$geo_scope = as.factor(db_summary2$geo_scope)
db_summary2$geo_scope = factor(db_summary2$geo_scope, levels = c('lmic','both','hic'))
ggplot(data = db_summary2,
       aes(axis1 = phylo_software_binary, axis2 = public_data_binary,axis3 = public_code_binary, 
           y = Freq)) +
  scale_x_discrete(limits = c("software", "data", "code"), expand = c(.2, .05)) +
  xlab("open science") +
  geom_alluvium(aes(fill = geo_scope)) +
  scale_fill_manual(values = c('#0c2c84','#7fcdbb', '#edf8b1'))+
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_classic() 

# 
# 
# col.chr.ids = c('hum_sample_size','livestock_sample_size', 'pet_sample_size','poultry_sample_size',
#                 'wildlife_sample_size','arachnida_sample_size', 'insecta_sample_size', 'env_domain',
#                 'biotic.env_sample_size','abiotic_sample_size')
# names(db)
# db[col.chr.ids] <- sapply(db[col.chr.ids],as.numeric)
# str(db)

db$total_samples = rowSums(db[,c('hum_sample_size','ani_sample_size', 'env_domain')], na.rm=TRUE)
db$env_bin = ifelse(db$env_domain>0,1,0)
db$ani_bin = ifelse(db$ani_sample_size>0,1,0)
db$human_bin = ifelse(db$hum_sample_size>0,1,0)
db$total_domains = rowSums(db[,c('env_bin','ani_bin','human_bin')], na.rm=T)

# some are 0?!? and some are 1?
# annotation in excel clarified there were 276 samples from human and environment but split not given
db$total_domains[19]<-2
db$total_samples[19]<-276
db$env_bin[19] = 1
db$human_bin[19] = 1

names(db)
ggplot(db, aes(x=total_samples, fill=study_aim_1)) +
  geom_histogram()+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="total samples", y = "publications")+
  scale_x_continuous(trans='log10')

# The 4 I’m about to paste should be good, papers by sequencer needs a tiny amount of re-working and inclusion of start-end dates for the various platforms. Christina will focus on this.
# Figure 5: histogram of total study’s sample size (irrespective of the domain) coloured by agent_type. Jayna would mind taking care of this?
#   Figure 6: time gap between last_sample_date and publication_year coloured by study_aim_1. Christina is that OK for you to do it?
#   Figure 7: histogram/test of phylo_results connected to total study’s sample size (irrespective of the domain). Jayna I’d leave this with you if OK.
# M&M + Results sections: for you both to write down elements connected with these analyses/figures.
# I’ll check that the EMBO journal’s guidelines suit this paper, we agreed on giving it a try. I’ll focus on formatting the article, the rest of the writing, suppl. mat., Table 3 (a summary of zoonotic agents and study aims).

names(db)
db$time_to_pub = db$publication_year-as.numeric(db$last_sample_date)
ggplot(db, aes(x=time_to_pub, fill=study_aim_1)) +
  geom_histogram()+
  theme(legend.position="top")+
  # facet_wrap(~study_aim_1, nrows = 5) +
  theme_classic()+
  labs(x="years from last collection to publication", y = "papers")
