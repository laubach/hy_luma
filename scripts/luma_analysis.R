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
    
    require(dplyr)
    require(plyr)
    require(purrr)
    require(readr)
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
####            2  Set Working Directory               ####
###########################################################

  ### 2.1 Set a working directory
    setwd("~/Git/luma")
  

  ### 2.2 Create path to a LUMA data folder
    luma_data_path <- "~/R/R_wd/fisi/project/hy_GR_global_DNA_meth/LUMA/data/"
    
    
  ### 2.3 Create path to Access Fise data folder
    access_data_path <- "~/R/R_wd/fisi/access_fisi_data/"
    
    
  ### 2.4 Create path to sample selection data folder (from repository query)
    select_data_path <- paste("~/R/R_wd/fisi/project/hy_GR_global_DNA_meth/",
                              "sample_selection/output/", sep = '')
    
###########################################################
####                  3  Import Data                   ####
###########################################################      
  
  ### 3.1 Import LUMA data files
    
    ## a) Create a list of file names (note the files are not s)
      files <- list.files(paste(luma_data_path), pattern = "*.csv")
  
    ## b) Import LUMA files and bind them row by row into a single data frame
      luma_raw <- files %>%
        map (function(x) read_csv(file.path(luma_data_path, x))) %>% 
        reduce(rbind)
    
        
  ### 3.2 Import Sample Selection table
    # This file is a record of all samples selected for a project. Samples
    # are generated from querying the bio_repository and tblDarting. This
    # file in the sample_selection output folder.
      samp_select <- read_csv(paste(select_data_path,
                                    "sample_request_03May2016.csv", sep = ''))
      
      
  ### 3.3 Import Access Fisi Tables
      
    ## a.1) Import table hyenas, which contains basic deomgraphic data
      # read in both Talek and Serena tblHyena and save as a data frame
        tblHyenas_T <- read_csv(paste(access_data_path,
                                      "tblHyenas_Talek_11Jan16.csv", sep = ''))
        tblHyenas_S <- read_csv(paste(access_data_path,
                                      "tblHyenas_serena_13Jan16.csv", sep = ''))
      
    ## a.2) combine tblHyenas data frames into one by appending as new rows
        tbl_hyenas <- bind_rows(tblHyenas_T, tblHyenas_S)
      
    ## a.3) clean up workspace by removing Talek and Serena tblHyenas
        rm(tblHyenas_T, tblHyenas_S)
        
    ## b) Import table Darting, which contains darting record data
      # read in tblHyenas access backend file and save as a data frame
        tbl_darting <- read_csv(paste(access_data_path,
                                      "tblDarting_3jan17.csv", sep = '')) 
        
    ## c) Import table Darting, which contains darting record data
      # read in tblHyenas access backend file and save as a data frame
        tbl_weather <- read_csv(paste(access_data_path,
                                      "tblweather_3jan17.csv", sep = '')) 
      
  
           
###########################################################
####                4  Massage Data                    ####
###########################################################   
