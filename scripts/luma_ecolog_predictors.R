###########################################################
########   Spotted Hyena Global DNA Methylation    ########
########         LUMA Ecological Predictors        ########
########             By: Zach Laubach              ########
########       last updated: 12 February 2017      ########
###########################################################

### PURPOSE: This code is desingned to analyze LUMA data.

### Special thanks to all contributors and commentors on stackoverflow 
### and other open source communities for creative and informative coding
### solutions.


  # Code Blocks
    # 1: Configure Workspace
    # 2: Set Working Directory
    # 3: Import Data
    # 4: Prepare Raw Data
    # 5: 
    # 6: 
    # 7: 
    # 8: 


###########################################################
####            1  Configure Workspace                 ####
###########################################################

  ### 1.1 clear global environment
    rm(list = ls())


  ### 1.2 Require - load a packages  into a current R session
    ## a) Data Manipulation and Descriptive Stats Packages
      require(dplyr)
      require(plyr)
      require(purrr)
      require(readr)
      require(reshape)
      require(reshape2)
      require(tidyr)
      options(gsubfn.engine = "R") #fixes tcltk bug; run before require sqldf
      require(sqldf)
    
    ## b) Graph Plotting and Visualization Packages
      require(ggplot2)
      require(sjPlot)
      require(sciplot)
      require(effects)
      require(plotrix)
      require(scatterplot3d)
      require(texi2dvi)
      require(stargazer)
      require(xtable)
      require(gridExtra)
    
      
    ## c) Modeling Packages
      # Modeling Packages
      require(lme4)
      require(arm)
      require(car)
      require(lsmeans)
      require(boot)
      require(bbmle)

    
  ### 1.3 Get Version and Session Info
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
    setwd("~/Git/fisi_lab/hy_luma")
  

  ### 2.2 Create path to a LUMA data folder
    luma_data_path <- paste("~/R/R_wd/fisi/project/hy_GR_global_DNA_meth/",
                            "LUMA/output/", sep = '')
    
    
  ### 2.3 Create path to sample selection data folder (from repository query)
    access_data_path <- paste("~/R/R_wd/fisi/project/access_fisi_data/",
                              sep = '')
    
    
       
###########################################################
####                  3  Import Data                   ####
###########################################################      
  
  ### 3.1 Import LUMA data files
    
    ## a) Import LUMA data, which has undergone QAQC (see R script 
      # luma_prep_analysis.R)
      luma_data <- read_csv(paste(luma_data_path,
                                    "luma_data.csv", sep = ''))
      
        
  ### 3.2 Import Sample Selection table
    # This file is a record of all samples selected for a project. Samples
    # are generated from querying the bio_repository and tblDarting. This
    # file in the sample_selection output folder.
      samp_select <- read_csv(paste(select_data_path,
                                    "sample_request_03May2016.csv", sep = ''))
      
    # clean up/standardize categorical age names
      samp_select$Age[grepl("s", samp_select$Age)]<-"subadult"    
      
  
           
###########################################################
####              4  Prepare Raw Data                  ####
###########################################################   


###########################################################
####                5  Assess Controls                 ####
###########################################################  
        
 
###########################################################
####                6  Build Data Set                  ####
###########################################################          
        
          
###########################################################
####      7  Descriptive Statistics (Univariate)       ####
###########################################################         

    
###########################################################
####            8  Identify Sample Re-Runs             ####
###########################################################         
 
###########################################################
####             9. Save Intermediate Tables           ####
####                  as Spreadsheets                  ####
###########################################################
          
  # Save intermediate tables as spreadsheets with a .cvs extension and today's
  # date. Files are saved in the 'data' folder or the 'output' folder
  # in the working directory.
  
                  
  ### 9.1 Set up date parameters
    # print today's date
      today <- Sys.Date()
      date <- format(today, format="%d%b%Y")
          
          
  ### 9.2 Generate File Names
    # For each table that will be saved as a .csv file, first generate a file 
    # name to save each table
          
    ## a) File name for sample_request table
      csv.file.name.luma <- paste (luma_data_out_path, "luma_data",
                                   ".csv", sep= "")   
    
    ## b) File name for sample_request table
      csv.file.name.re_runs <- paste (luma_data_out_path, "sample_re_runs",
                                      date, ".csv", sep= "")  
      
      
  ### 9.3 Save Tables 
    # Save each data frame as a .csv file (a spreadsheet/table) into the 
    # data folder in the working directory.
      
    ## a) Save luma_data_no_out table
      write.csv (luma_data_no_out, file = csv.file.name.luma)
      
    ## b) Save re_runs table
      write.csv (re_runs, file = csv.file.name.re_runs)
        
        
        
        
        
        
        
        
            
