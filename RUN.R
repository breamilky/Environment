################################################################################
# RUN_ALL.R  
# Usage:  source("RUN_ALL.R")
################################################################################

t_start <- Sys.time()
cat("Pipeline started:", format(t_start), "\n\n")

for (p in c("Output", "Output/RDS_files", "Output/Figures", "Output/Maps",
            "Output/Tables", "Output/tables", "Output/Robustness/Tables",
            "Output/PM10_Robustness/Tables", "Output/DoseResponse",
            "Output/EventStudy", "Output/Counterfactual",
            "Data/Mid_process_data")) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

#-------------------------------------------------------------------------------
# STAGE 1 - Data
#-------------------------------------------------------------------------------
source("01_Clean_Merge.R")              # Table 1                       

#-------------------------------------------------------------------------------
# STAGE 2 - Core results  
#-------------------------------------------------------------------------------
source("0202_first.R")                  # Tables 2-6, 11, 12, 13, 14,
# A1a, A1b, A2; Figure 2       
source("00_Bridge_RDS.R")               # RDS filename bridge           
source("0202_second.R")                 # Tables 7-10, A3-A6 (PM2.5)     
source("01_Feb2_4_PM10.R")              # PM10 columns                   
source("01_Feb2_3.R")                   # Figure 1, Figure A1            
source("96_Pairs_Bootstrap.R")          # Table A2 row 1 (PRIMARY SE)

#-------------------------------------------------------------------------------
# STAGE 3 - Gap scripts   (read the merged CSV only; safe to run in a
#                          second R session in parallel with Stage 2)
#-------------------------------------------------------------------------------
source("91_Exclude2020.R")              # Table 11, COVID row         
source("92_NegativeBinomial.R")         # Table 11, NegBin row        
source("93_Regional_Heterogeneity.R")   # Section 5.3                    
source("94_FireDays_Above5.R")          # Table 12, final column         
source("95_Welfare_Sensitivity.R")      # Section 6.7                    

#-------------------------------------------------------------------------------
# STAGE 4 - Collect
#-------------------------------------------------------------------------------
source("99_Build_Publication_Folder.R") # Publication/             

cat("\nFinished in",
    round(as.numeric(difftime(Sys.time(), t_start, units = "hours")), 2),
    "hours\n")

