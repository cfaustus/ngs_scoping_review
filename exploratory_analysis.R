# code to explore and visualize dataset

# libraries
library(ggplot2)
library(ggalluvial)
library(dplyr)

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

#sample_cols = grep('sample_size',names(db4)) # identifying sample size columns - WRONG double counting
db4$total_samples = rowSums(db4[,c('hum_sample_size','ani_sample_size', 'env_domain')], na.rm=TRUE)
db4$env_bin = ifelse(db4$env_domain>0,1,0)
db4$ani_bin = ifelse(db4$ani_sample_size>0,1,0)
db4$human_bin = ifelse(db4$hum_sample_size>0,1,0)
db4$total_domains = rowSums(db4[,c('env_bin','ani_bin','human_bin')], na.rm=T)

# some are 0?!? and some are 1?
# annotation in excel clarified there were 276 samples from human and environment but split not given
db4$total_domains[19]<-2
db4$total_samples[19]<-276
db4$env_bin[19] = 1
db4$human_bin[19] = 1



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

# https://cran.r-project.org/web/packages/ggalluvial/vignettes/ggalluvial.html
head(db4)
db4_allu = is_alluvia_form(as.data.frame(db4), axes = 1:3, silent = TRUE)
ggplot(data = db4,
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
  

db4$ngs_platform_short_primary
head(db4)
db4_summary = as.data.frame(db4 %>% group_by(ngs_platform_short_primary, publication_year, ngs_platform_company) %>%
                                summarise(papers = n()))
db4_summary$ngs_platform_short_primary=as.factor(db4_summary$ngs_platform_short_primary)

db4_summary = db4_summary[1:45,]
db4_summary$ngs_platform_short_primary
unique(db4$ngs_platform_long)

db4_summary$ngs_platform_short_primary = as.factor(db4_summary$ngs_platform_short_primary)

levels(db4_summary$ngs_platform_short_primary)
db4_summary$ngs_platform_short_primary <- factor(db4_summary$ngs_platform_short_primary, 
                levels = c("Illumina_NovaSeq","BGI_MGISEQ2000", 
                           "Illumina_iSeq100","BGI_BGISEQ50", 
                           "BGI_BGISEQ500", "Illumina_NextSeq", 
                           "Illumina_MiSeq", "IonTorrent_IonPGM", 
                           "Illumina_HiSeq", "Illumina_GAIIx", 
                            "Roche_454"))

ggplot(db4_summary, aes(x = publication_year, y = ngs_platform_short_primary )) +
  geom_point(aes(size = papers, color = ngs_platform_company))+ # aes(dotsize = papers),
  scale_x_continuous(breaks=seq(2011,2022,1))+
  theme_classic()+
  labs(x="publication year", y = "sequencer")+
  scale_size_continuous(breaks = c(min(db4_summary$papers),
                                    3,
                                    max(db4_summary$papers)),
                         labels = c("1",
                                    "3",
                                    "9"))
# need to combine short read and long read sequencing

# is within host diversity studies predicted by mutation rate?
# maybe not relevant b/c no host diversity

# alluvial plot to make
# open data 
# open software
# open code
# colored by country?
db4$public_data_binary = !db4$public_data == 'N/A'
db4$public_code_binary = db4$public_code == 'yes'
db4$phylo_software_binary = !db4$phylo_software == 'Geneious'
db4$sampling_geo

db4_summary2 = as.data.frame(db4 %>% group_by(public_data_binary, public_code_binary, phylo_software_binary,sampling_geo) %>%
                              summarise(Freq = n()))

ggplot(data = db4_summary2,
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

db5 = read.csv('data/NGSscoping_dbv5.csv')
names(db5)<-tolower(names(db5))

col.chr.ids = c('hum_sample_size','livestock_sample_size', 'pet_sample_size','poultry_sample_size',
                'wildlife_sample_size','arachnida_sample_size', 'insecta_sample_size', 'env_domain',
                'biotic.env_sample_size','abiotic_sample_size')
names(db5)
db5[col.chr.ids] <- sapply(db5[col.chr.ids],as.numeric)
str(db5)

db5$total_samples = rowSums(db5[,c('hum_sample_size','ani_sample_size', 'env_domain')], na.rm=TRUE)
db5$env_bin = ifelse(db5$env_domain>0,1,0)
db5$ani_bin = ifelse(db5$ani_sample_size>0,1,0)
db5$human_bin = ifelse(db5$hum_sample_size>0,1,0)
db5$total_domains = rowSums(db5[,c('env_bin','ani_bin','human_bin')], na.rm=T)

# some are 0?!? and some are 1?
# annotation in excel clarified there were 276 samples from human and environment but split not given
db5$total_domains[19]<-2
db5$total_samples[19]<-276
db5$env_bin[19] = 1
db5$human_bin[19] = 1

names(db5)
ggplot(db5, aes(x=total_samples, fill=agent_hse_cat)) +
  geom_histogram()+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="total samples", y = "domains")+
  scale_x_continuous(trans='log10')

# The 4 I’m about to paste should be good, papers by sequencer needs a tiny amount of re-working and inclusion of start-end dates for the various platforms. Christina will focus on this.
# Figure 5: histogram of total study’s sample size (irrespective of the domain) coloured by agent_type. Jayna would mind taking care of this?
#   Figure 6: time gap between last_sample_date and publication_year coloured by study_aim_1. Christina is that OK for you to do it?
#   Figure 7: histogram/test of phylo_results connected to total study’s sample size (irrespective of the domain). Jayna I’d leave this with you if OK.
# M&M + Results sections: for you both to write down elements connected with these analyses/figures.
# I’ll check that the EMBO journal’s guidelines suit this paper, we agreed on giving it a try. I’ll focus on formatting the article, the rest of the writing, suppl. mat., Table 3 (a summary of zoonotic agents and study aims).

names(db5)
db5$time_to_pub = db5$publication_year-as.numeric(db5$last_sample_date)
ggplot(db5, aes(x=time_to_pub, fill=study_aim_1)) +
  geom_histogram()+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="years from last collection to publication", y = "papers")
  scale_x_continuous(trans='log2')
