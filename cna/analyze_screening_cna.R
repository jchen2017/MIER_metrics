library(cna)
library(readxl)

datadir<-"./"

# read data
exceldata = read_excel(paste(datadir,"clinic_characteristis.xlsx", sep=""), sheet = "factors_cna")                                                                            
dfdata = data.frame(exceldata)

dfdata

# rename variables
URBAN <- dfdata$Service_area == "Urban"
URBAN
Medical <- dfdata$Specialty == "Medical"
Medical
LARGE <- dfdata$Encounters_per_year == 1
LARGE
HIGH_ALT_EXP <- dfdata$Team_alert_exposure == 1
HIGH_ALT_EXP
HIGH_ALT_ENC_REL <- dfdata$Alert_Encounter_relevance == 1 
HIGH_ALT_ENC_REL
HI_SCREEN_RATE <- dfdata$Screening_alert_completion_rate == 1
HI_SCREEN_RATE

df <- data.frame(dfdata$Clinic, URBAN, Medical, LARGE, HIGH_ALT_EXP, HIGH_ALT_ENC_REL, HI_SCREEN_RATE)
df


# Data resulting from theoretical factor selection
dat <- df[,c("HI_SCREEN_RATE","URBAN", "LARGE", "Medical", "HIGH_ALT_EXP", "HIGH_ALT_ENC_REL")]
dat

# CNA
result <-subset(asf(mvcna(dat, ordering=list("HI_SCREEN_RATE"), strict = T, con.msc=0.6, con=0.8,cov=0.8, maxstep=c(3,3,3), rm.dup.factors = FALSE)), outcome=="HI_SCREEN_RATE=1")
result
sol.list <- list(result,dat)

# List of solutions with corresponding datasets
sink("cna_results.txt")
sol.list
sink()

