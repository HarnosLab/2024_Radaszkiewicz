# REPRINT analysis

# Libraries
library(dplyr)
library(here)

# Data import
data <- read.delim(here('data', '20240502_report.txt'))

subset <- data %>%
  select(UniProt_Protein.ID,
         UniProt_Gene.Names..primary.,
         UniProt_Length,
         starts_with("Precursor.Id"))

colnames(subset) <- gsub("Precursor.Id.", "", colnames(subset))

subset <- subset %>%
  rename('PROTID' = 'UniProt_Protein.ID',
         'GENEID' = 'UniProt_Gene.Names..primary.',
         'PROTEIN' = 'UniProt_Length',
         'CTRL_1_NUMSPECSTOT' = '1_REP2_CTRL' ,
         'PRK1_1_NUMSPECSTOT' = '2_REP2_PRK1' ,
         'PRK2_1_NUMSPECSTOT' = '3_REP2_PRK2' ,
         'PRK3_1_NUMSPECSTOT' = '4_REP2_PRK3' ,
         'CTRL_2_NUMSPECSTOT' = '5_REP3_CTRL' ,
         'PRK1_2_NUMSPECSTOT' = '6_REP3_PRK1' ,
         'PRK2_2_NUMSPECSTOT' = '7_REP3_PRK2' ,
         'PRK3_2_NUMSPECSTOT' = '8_REP3_PRK3' ,
         'CTRL_3_NUMSPECSTOT' = '9_REP4_CTRL' ,
         'PRK1_3_NUMSPECSTOT' = '10_REP4_PRK1' ,
         'PRK2_3_NUMSPECSTOT' = '11_REP4_PRK2' ,
         'PRK3_3_NUMSPECSTOT' = '12_REP4_PRK3' ,
         'CTRL_4_NUMSPECSTOT' = '13_REP5_CTRL' ,
         'PRK1_4_NUMSPECSTOT' = '14_REP5_PRK1' ,
         'PRK2_4_NUMSPECSTOT' = '15_REP5_PRK2' ,
         'PRK3_4_NUMSPECSTOT' = '16_REP5_PRK3' ,
         'CTRL_5_NUMSPECSTOT' = '17_REP6_CTRL' ,
         'PRK1_5_NUMSPECSTOT' = '18_REP6_PRK1' ,
         'PRK2_5_NUMSPECSTOT' = '19_REP6_PRK2' ,
         'PRK3_5_NUMSPECSTOT' = '20_REP6_PRK3')

load(here('data', "genenames_20240502.RData"))
update_subset <- checkGeneSymbols(subset$GENEID, species = "human", map = genenames_newest, unmapped.as.na = FALSE)
subset <- left_join(subset, update_subset %>% select(x, Suggested.Symbol), by = c("GENEID" = "x"))
subset$Suggested.Symbol[subset$Suggested.Symbol == "MASP2 /// BRD8 /// C11orf58 /// KIFAP3"] <- "KIFAP3"
subset$Suggested.Symbol[subset$Suggested.Symbol == "C10orf88 /// SLC38A4"] <- "C10orf88"
subset$GENEID <- subset$Suggested.Symbol
subset$Suggested.Symbol <- NULL

vec <- colnames(subset)
vec[1:3] <- "--"
vec <- gsub("(\\_\\d+\\_[A-Z]+)", "", vec)
vec <- gsub("CTRL", "C", vec)

# join the row with dataframe
joined <- rbind (vec, subset)

# write txt file
write.table(joined, here('data', 'REPRINT-data-input.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
