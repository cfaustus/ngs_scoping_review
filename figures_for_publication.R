# code to explore and visualize dataset

# libraries
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(plyr)
library(wesanderson)
library(patchwork)

# importing and correctly formatting data
db <- read.csv("data/NGSscoping-db-final.csv", header = TRUE)
names(db)<-tolower(names(db))
db$biotic.env_sample_size
colnames(db)[29] <- 'biotic_env_sample_size'

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

# annotation in excel clarified there were 276 samples from human and environment but split not given
ind = which(db$total_samples == 0)
db$total_domains[ind]<-2
db$total_samples[ind]<-276
db$env_bin[ind] = 1
db$human_bin[ind] = 1

# temporal trends
sort(unique(db$publication_year))
# summarizing agent type into brader categories
db = db |> mutate(agent_type = ifelse(str_detect(agent_type, "bacterium") == TRUE, "bacterium", agent_type)) |>
  mutate(agent_type = ifelse(str_detect(agent_type, "virus") == TRUE, "virus", agent_type))

# figure 1 - annual trend by pathogen type
ggplot(db, aes(x=publication_year, fill=agent_type)) +
  geom_histogram(col = 'white', position="stack", bins = 12)+
  theme(legend.position="top")+
  scale_fill_ordinal()+
  theme_classic()+
  labs(x="year published", y = "publications")+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,23,2))+
  scale_x_continuous(breaks=seq(2011,2022,1))

# extra figure - annual trend by NGS sequencer
ggplot(db, aes(x=publication_year, fill=ngs_platform_primary)) +
  geom_histogram(col = 'white', position="stack", bins = 12)+
  theme(legend.position="top")+
  theme_classic()+
  labs(x="year published", y = "publications")+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,23,2))+
  scale_x_continuous(breaks=seq(2011,2022,1))

# extra figure - annual trend by study aim - in SI
ggplot(db, aes(x=publication_year, fill=study_aim_1)) +
  geom_histogram(col = 'white', position="stack", bins = 12)+
  theme(legend.position="top")+
  facet_wrap(~study_aim_1, ncol = 1)+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x="year published", y = "publications")+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,23,2))+
  scale_x_continuous(breaks=seq(2011,2022,1))


# NGS Sequencers over time
db_platform = as.data.frame(db %>% group_by(ngs_platform_primary, publication_year) %>%
                              dplyr::summarise(papers = n()))
db_platform$ngs_platform_primary=as.factor(db_platform$ngs_platform_primary)
missing = grep('N/S', db_platform$ngs_platform_primary)
db_platform= db_platform %>%  filter(!row_number() %in% missing)

seq = read.csv('data/sequencing_tech_long.csv', header = TRUE)
names(seq)
seq_trun = seq[,c('company', 'ngs_platform_primary','type')]
db_platform <- merge(x=db_platform,y=seq_trun, by = "ngs_platform_primary", all = TRUE)

db_platform$type = as.factor(db_platform$type)
db_platform$type <- factor(db_platform$type, 
                           levels = c("MGISEQ2000", "iSeq100","NovaSeq",
                                      "BGISEQ50", 
                                      "BGISEQ500", "MinION","NextSeq",'RSII',
                                      "MiSeq", "IonPGM", 
                                      "HiSeq", "GAIIx", 
                                      "454"))
seq_col = wes_palette('Darjeeling1', 6, type = c("continuous"))

ggplot(db_platform, aes(x = publication_year, y = type)) +
  geom_point(aes(size = papers, color = company))+ # aes(dotsize = papers),
  scale_x_continuous(breaks=seq(2006,2022,1))+
  theme_classic()+
  labs(x="publication year", y = "sequencer")+
  scale_color_manual(values = seq_col)+
  geom_line(data = seq, aes(y = type,x = publication_year, color = company))+
  scale_size_continuous(breaks = c(min(db_platform$papers),
                                   3,
                                   max(db_platform$papers)),
                        labels = c("1",
                                   "3",
                                   "9"))

# alluvial plot to make
db$public_data_binary = !db$public_data == 'N/A'
db$public_code_binary = db$public_code == 'yes'
db$phylo_software_binary = !db$phylo_software == 'Geneious'
db$sampling_geo

db_summary2 = as.data.frame(db %>% group_by(public_data_binary, public_code_binary, phylo_software_binary,sampling_geo) %>%
                              dplyr::summarise(Freq = n()))

ggplot(data = db_summary2,
       aes(axis1 = as.factor(phylo_software_binary), axis2 = public_data_binary,axis3 = public_code_binary, 
           y = Freq)) +
  scale_x_discrete(limits = c("software", "data", "code"), expand = c(.2, .05)) +
  xlab("open science") +
  geom_alluvium(aes(fill = sampling_geo)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() 

db$geo_scope = ifelse(db$lic_lmic == 'yes' & db$umic_hic == 'yes', 'both',
                      ifelse(db$lic_lmic == 'yes', 'lmic', 'hic'))

db$public_code
db$public_data_binary = !db$public_data == 'N/A'
db$public_code_binary = db$public_code == 'yes'
db$phylo_software_binary = db$assembly.mapping_license == 'free'
db_summary2 = as.data.frame(db %>% group_by(public_data_binary, public_code_binary, phylo_software_binary,geo_scope) %>%
                              dplyr::summarise(Freq = n()))

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


db$time_to_pub = db$publication_year-as.numeric(db$last_sample_date)
aim_col = wes_palette('AsteroidCity3', 5, type = c("continuous"))
db_me <- ddply(db, .(study_aim_1), numcolwise(median, na.rm = TRUE))
db_mean <- ddply(db, .(study_aim_1), numcolwise(mean, na.rm = TRUE))

ggplot(db, aes(x=time_to_pub, fill=study_aim_1)) +
  geom_density(alpha = 0.80)+
  theme(legend.position="top")+
  facet_wrap(~study_aim_1, nrow= 5) +
  scale_fill_manual(values = aim_col)+
  theme_classic()+
  geom_vline(data=db_me, aes(xintercept=time_to_pub, colour=study_aim_1),
             linetype="solid", size=0.5) +
  geom_vline(data=db_mean, aes(xintercept=time_to_pub, colour=study_aim_1),
             linetype="dashed", size=0.5) +
  theme(legend.position = "none")+
  scale_color_manual(values = aim_col)+
  labs(x="years from last collection to publication", y = "papers")


study_aim_dens = ggplot(db, aes(x=time_to_pub, fill=study_aim_1)) +
  geom_density(alpha = 0.80)+
  facet_wrap(~study_aim_1, nrow= 5) +
  scale_fill_manual(values = aim_col)+
  theme_classic()+
  geom_vline(data=db_me, aes(xintercept=time_to_pub, colour=study_aim_1),
             linetype="solid", size=0.5) +
  geom_vline(data=db_mean, aes(xintercept=time_to_pub, colour=study_aim_1),
             linetype="dashed", size=0.5) +
  theme(legend.position = "none")+
  scale_color_manual(values = aim_col)+
  labs(x="years from last sample collection to publication", y = "density")

study_aim_hist = ggplot(db, aes(x=publication_year, fill=study_aim_1)) +
  geom_histogram(col = 'white', position="stack", bins = 12)+
  facet_wrap(~study_aim_1, nrow= 5) +
  scale_fill_manual(values = aim_col)+
  theme_classic()+
  labs(x="year published", y = "publications")+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,10,2))+
  theme(legend.position = "none")+
  scale_x_continuous(breaks=seq(2011,2022,2))

# Labels for each plot
study_aim_hist + study_aim_dens + plot_annotation(tag_levels = "A")
