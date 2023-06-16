library(readr)
ngsscop_dbv4 <- read_csv("data/NGSscoping_dbv4.csv")
View(ngsscop_dbv4)
head(ngsscop_dbv4)

hist(log10(ngsscop_dbv4$ANI_SAMPLE_SIZE+1))
hist(log10(as.numeric(ngsscop_dbv4$HUM_SAMPLE_SIZE)))

# need to fix column name importations .. not showing up as correct values
