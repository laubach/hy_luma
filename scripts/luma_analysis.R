###########################################################
########   Spotted Hyena Global DNA Methylation    ########
########               LUMA Analysis               ########
########             By: Zach Laubach              ########
########       last updated: 21 January 2017       ########
###########################################################

### PURPOSE: This code is desingned to analyze LUMA data.

### Special thanks to all contributors and commentors on stackoverflow 
### for creative code suggestions and solutions.


    # Code Blocks
    # 1: Load Packages
    # 2: Import Data
    # 3: Set Global Values
    # 4: Query and Select Sessions 
    # 5: Save and Write Out Data File to .csv



###########################################################
####                 1  Load Packages                  ####
###########################################################

  # clear global environment
    rm(list = ls())

  # Require - load a packages  into a current R session


  ### 1.1  Data Manipulation and Descriptive Stats Packages
    require(plyr)
    require(dplyr)
    require(reshape)
    require(reshape2)
    require(tidyr)
    options(gsubfn.engine = "R") #fixes tcltk bug; run before require sqldf
    require(sqldf)

    
  ### 1.2 Get Version and Session Info
    R.Version()
    sessionInfo()    
    # Ran with...
    #R version 3.3.1 (2016-06-21)
    #Platform: x86_64-apple-darwin13.4.0 (64-bit)
    #Running under: OS X 10.12 (Sierra)

###########################################################
####                  2  Import Data                   ####
###########################################################

  ### 2.1 # Set a working directory
    #setwd("~/R/R_wd/basic_tools/maternal_care") # when wd is in home folder
    # set workding directory on L: Drive
    setwd("~/Git/fisi_lab/luma_hyena")
###test again
    
    
    ### let's try this again 
    #### hmm
