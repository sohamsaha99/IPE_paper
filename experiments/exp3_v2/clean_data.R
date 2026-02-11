rm(list = ls())
# Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(survival)
# library(SurvRegCensCov)

# Read GBM data --- AVAGLIO
df_GBM_trt <- read.csv("./experiments/exp3_v2/GBM_files/AVAGLIO_CleanRR.csv", header = TRUE, stringsAsFactors = FALSE)
# Only keep necessary columns
df_GBM_trt <- df_GBM_trt %>%
  select(Patient, age, sex, kps90_100, eor, mgmt, os_status, os, pfs) %>%
  rename(kps = kps90_100) %>% # Binary: 0: 0-80, 1: 90-100
  filter(eor %in% c("GTR", "STR")) %>%
  mutate(eor = factor(eor)) %>% # Levels: BX, GTR, STR
  mutate(os = os / 30, pfs = pfs / 30) # Measure `os` in months instead of days
# Additionally, we note that os_status = 0 implies alive at last follow-up,
# sex = 0 implies female
# Delete rows with NA (missing) entry
df_GBM_trt <- df_GBM_trt %>% drop_na(age, sex, kps, eor, mgmt, os_status, os)
# After above operations, there should be 351 rows and 9 columns
print(dim(df_GBM_trt))


df_GBM_trt %>% write.csv(file = "./experiments/exp3_v2/GBM_files/GBM_trt.csv", row.names = FALSE)

# Read GBM data --- DFCI
df_GBM_con <- read.csv("./experiments/exp3_v2/GBM_files/DFCI_CleanRR_3-7-21.csv", header = TRUE, stringsAsFactors = FALSE)
# Only keep necessary columns
df_GBM_con <- df_GBM_con %>%
  select(Patient_ID, age, sex, kps90_100, eor, mgmt, os_status, os_rtstart, pfs_rtstart) %>%
  rename(Patient = Patient_ID, kps = kps90_100, os = os_rtstart, pfs = pfs_rtstart) %>%
  filter(eor %in% c("GTR", "STR")) %>%
  mutate(eor = factor(eor)) %>% # Levels: BX, GTR, STR
  mutate(os = os / 30, pfs = pfs / 30) # Measure `os` in months instead of days
# Additionally, we note that os_status = 0 implies alive at last follow-up,
# sex = 0 implies female
# Delete rows with NA (missing) entry
df_GBM_con <- df_GBM_con %>% drop_na(age, sex, kps, eor, mgmt, os_status, os)
# After above operations, there should be 351 rows and 9 columns
print(dim(df_GBM_con))


df_GBM_con %>% write.csv(file = "./experiments/exp3_v2/GBM_files/GBM_con.csv", row.names = FALSE)

