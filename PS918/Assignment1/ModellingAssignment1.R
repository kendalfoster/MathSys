##############################
##### PS918 Assingment 1 #####
####### Kendal Foster ########
##############################


# Preamble ----------------------------------------------------------------
## Need to install packages "rtdists" and "tidyverse" ##
library("rtdists")
library("tidyverse")
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom", 
                  panel.grid.major.x = element_blank(), 
                  panel.grid.minor.x = element_blank()))

## Set Working Directory ##
#setwd("G:/My Drive/Scholastic/Grad School/Warwick/Term 2/PS918/Workshops/Week 3-4")
setwd("/home/kfoster/Documents/PS918/Workshops/Week 3-4")

## Import Data ##
med <- read_csv("medical_dm.csv")



# Simple Model; Groups -----------------------------------------------------
## Data Processing
# Split individuals into Novice and Expert categories
med_novice <- filter(med, group == "novice")
med_expert <- filter(med, group == "experienced" | group == "inexperienced")





