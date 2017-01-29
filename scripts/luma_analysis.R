###########################################################
########   Spotted Hyena Global DNA Methylation    ########
########               LUMA Analysis               ########
########             By: Zach Laubach              ########
########       last updated: 21 January 2017       ########
###########################################################

### PURPOSE: This code is desingned to analyze LUMA data.

### Special thanks to all contributors and commentors on stackoverflow 
### and other open source communities for creative and informative coding
### solutions.


    # Code Blocks
    # 1: Configure Workspace
    # 2: Set Working Directory
    # 3: Import Data
    # 4: Massage Data
    # 5: Assess Controls 
    # 6: 



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
      #require(texi2dvi)
      #require(stargazer)
      #require(xtable)
      
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

  ### 4.1 Remove Failed Wells
      # update luma_raw by removing any well that failed due to high CV, 
      # non-specific peaks or low pyrosequencing signal (see meta data file
      # for specific details)
      luma_raw <- filter (luma_raw, analysis_status == "y")
        
        
  ### 4.2 Controls
    ## a) Subset Controls 
      # subset luma_raw into a data frame that contains all LUMA controls
        luma_cntrl <- filter (luma_raw, grepl("hy", sample_ID))
      
    ## b) Gather luma_cntrl
      # gather the luma_cntrl data so that methylation duplicates are in 
      # long format
        luma_cntrl <- luma_cntrl %>%
          gather (meth_dup, dup, dup1:dup2)
        
    
  ### 4.3 Linearization Standards
    ## a) Subset Linearization Standards
      #subset luma_raw into a data frame that contains all LUMA linearization 
      #standards
        luma_linearz <- filter (luma_raw, grepl("lam", sample_ID)) 
      
      
  ### 4.4 Data
    ## a)  Subset Data
      # subset luma_raw into a data frame that contains all LUMA data. 
      # sample_IDs are numbers
        luma_data <- filter (luma_raw, grepl("^[[:digit:]]", sample_ID)) 
        


###########################################################
####                5  Assess Controls                 ####
###########################################################  
        
  ### 5.1 Intra-Class CV
    ## a) hy_pool and hy_100% Intra-Class CV
      # use ddply to calulate the intra-class CV (within plate variation)
      # for pooled samples.
        intra_CV <- ddply (luma_cntrl, .(plate_rxn_ID, sample_ID),
                           summarise,
                           N = sum(!is.na(methylation)),
                           mean = round (mean(methylation, na.rm = T), 2),
                           median =  round (quantile(methylation, c(.5),
                                                     na.rm = T), 2),
                           sd = round (sd(methylation, na.rm = T), 2),
                           cv = 100* (sd/mean))
    
       
                               
        
  ### 5.2 Inter-Class CV
    ## a) hy_pool and hy_100% Inter-Class CV
      # use ddply to calulate the inter-class CV (between plate variation) 
      # for pooled samples.  
        inter_CV <- ddply (luma_cntrl, .( sample_ID), summarise,
                           N = sum(!is.na(methylation)),
                           mean = round (mean(methylation, na.rm = T), 2),
                           median =  round (quantile(methylation, c(.5),
                                                     na.rm = T), 2),
                           sd = round (sd(methylation, na.rm = T), 2),
                           cv = 100* (sd/mean))
        
        
  ### 5.3 Check Linearization
    ## a) Predicted Values
      # arrange the luma_linearz descending order by plate_pos_seq
        luma_linearz <- arrange(luma_linearz, plate_pos_seq)
      
      # make a vector of predicted and append those to the linearization 
      # data frame
        luma_linearz$pred_meth <- as.numeric(c(100, 75, 65, 50, 25, 0))
        
    ## b) Graph Linearization Data
      # Graph the actual and predicted methylation values from Lambda phage
      # DNA methylation standard curves
        ggplot(luma_linearz, aes (x = sample_ID, group=1)) +
                 geom_line(aes(y = pred_meth, color = "pred_meth"), 
                           size = 1) +
                 geom_line(aes(y = methylation, color = "methylation"), 
                           size = 1) +
          scale_color_manual(values = c("red", "dark grey")) +
          # scale_color_hue(l = 50) + # Use a slightly darker palette than normal
          scale_x_discrete(limits = c(luma_linearz$sample_ID)) +
          labs (title = "LUMA Global 
                DNA Mehtylation Standard Curve") +
          ylab ("% Global DNA Methylation") +
          xlab ("Standard Curve")
        
    ## c) Save Plot
      # use ggsave to save the linearization plot
        ggsave("luma_linearization.pdf", plot = last_plot(), device = NULL, 
               path = "./output", scale = 1, width = 6, height = 4, 
               units = c("in"), dpi = 300, limitsize = TRUE)

        
  ### 5.4 Check Plates for Drift    
    ## a) Graph Controls Across Plates
      # Graph the control on each reaction plate to assess for any plate drift
        ggplot(luma_cntrl, aes (x = as.factor(plate_pos_factor), 
                                y = methylation, color = sample_ID, 
                                group = sample_ID)) +
          geom_line(size = 1) +
          facet_grid(. ~ plate_rxn_ID) +
          scale_color_manual(values = c("red", "dark grey")) +
          scale_x_discrete() +
          labs (title = "LUMA Global DNA Mehtylation 
                Plate by Plate Controls") +
          ylab ("% Global DNA Methylation") +
          xlab ("Plate Positions")
    
    ## b) Save Plot    
      # use ggsave to save the control drift plot
        ggsave("control_drift.pdf", plot = last_plot(), device = NULL, 
               path = "./output", scale = 1, width = 9, height = 3, 
               units = c("in"), dpi = 300, limitsize = TRUE)   
        
    ## c) Graph Controls Across Plates 
      # Graph the control on each reaction plate to assess for any plate drift 
        ggplot(luma_cntrl, aes (x = well, y = methylation,
                                color = sample_ID, group = sample_ID)) +
          geom_line(size = 1) +
          facet_grid(plate_rxn_ID ~ sample_ID) +
          scale_color_manual(values = c("red", "dark grey")) +
          scale_x_discrete() +
          labs (title = "LUMA Global DNA Mehtylation 
                Plate by Plate Controls") +
          ylab ("% Global DNA Methylation") +
          xlab ("Plate Positions")
        
    ## d) Save Plot
      # use ggsave to save the linearization plot
        ggsave("double_panel_control_drift.pdf", plot = last_plot(), device = NULL, 
               path = "./output", scale = 1, width = 7, height = 7, 
               units = c("in"), dpi = 300, limitsize = TRUE)   
        
    ## c) Graph Controls Across Plates 
      # Graph the controls on each reaction plate to assess for any plate drift
      # least squares regression is used for the fit function 
        ggplot(luma_cntrl, aes (x = well, y = methylation,
                                color = sample_ID, group = sample_ID)) +
          geom_point(size = 1) +
          geom_smooth(method = lm, se = F) +
          facet_grid(plate_rxn_ID ~ sample_ID) +
          scale_color_manual(values = c("red", "dark grey")) +
          scale_x_discrete() +
          labs (title = "LUMA Global DNA Mehtylation 
                Plate by Plate Controls") +
          ylab ("% Global DNA Methylation") +
          xlab ("Plate Positions")
  
    ## d) Save Plot
      # use ggsave to save the linearization plot
        ggsave("double_panel_control_drift.pdf", plot = last_plot(), device = NULL, 
               path = "./output", scale = 1, width = 7, height = 7, 
               units = c("in"), dpi = 300, limitsize = TRUE) 
   
        
    ## e) Subset Controls 
      # subset luma_raw into a data frame that contains one of the LUMA control
      # dupliates
        luma_pool <- filter (luma_cntrl, grepl("pool", sample_ID))     
        luma_pool <- filter (luma_pool, grepl("dup1", meth_dup))
    
    ## f) Plate by Plate Linear Regression
      # For each plate model the possible drift between controls; output is
      # a list of lm objects
        plate_drift = dlply(luma_pool, .(plate_rxn_ID), lm,
                            formula = methylation ~ plate_pos_factor)

      # For use in plyr's ldply() function, the utility function should
      # return a data frame. We save some effort in simple linear regression
      # by noting that the two-sided p-value of the t-test of zero slope is the
      # same as that of the overall F test:
        extractfun <- function(m) {
          cf <- coef(m)
          tinfo <- summary(m)$coefficients[2, c(2, 4)]
          r2 <- summary(m)$r.squared
          data.frame(intercept = cf[1], slope = cf[2], n = length(resid(m)),
                     slope.se = tinfo[1], pval = tinfo[2], Rsq = r2)
        }
        
      # Take a list (of models) as input and output a data frame:
        ldply(plate_drift, extractfun)
      