###########################################################
########   Spotted Hyena Global DNA Methylation    ########
########               LUMA Analysis               ########
########             By: Zach Laubach              ########
########       last updated: 12 February 2017       ########
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
    # 5: Assess Controls 
    # 6: Build Data Set
    # 7: Descriptive Statistics (Univariate)
    # 8: Bi-Variate Analysis


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
    luma_data_path <- "~/R/R_wd/fisi/project/hy_GR_global_DNA_meth/LUMA/data/"
    
    
  ### 2.3 Create path to sample selection data folder (from repository query)
    select_data_path <- paste("~/R/R_wd/fisi/project/hy_GR_global_DNA_meth/",
                              "sample_selection/output/", sep = '')
 
    
  ### 2.4 Create path to a LUMA data folder
    luma_data_out_path <- paste("~/R/R_wd/fisi/project/hy_GR_global_DNA_meth/",
                                "LUMA/output/", sep = '')  
    
       
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
      
    # clean up/standardize categorical age names
      samp_select$Age[grepl("s", samp_select$Age)]<-"subadult"    
      
  
           
###########################################################
####              4  Prepare Raw Data                  ####
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
        
    ## b) save the data frame of summary stats out as a pdf into output file
        pdf("output/intra_CV.pdf", height = 6, width = 7)
        grid.table(intra_CV)
        dev.off()
        
        
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
        
    ## b) save the data frame of summary stats out as a pdf into output file
        pdf("output/inter_CV.pdf", height = 6, width = 7)
        grid.table(inter_CV)
        dev.off()
        
              
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
  
    ## b) Save Plot
      # use ggsave to save the linearization plot
        ggsave("double_panel_control_drift.pdf", plot = last_plot(), device = NULL, 
               path = "./output", scale = 1, width = 7, height = 7, 
               units = c("in"), dpi = 300, limitsize = TRUE) 
   
    ## e) Subset Controls 
      # subset luma_raw into a data frame that contains one of the LUMA control
      # duplicates
        luma_pool <- filter (luma_cntrl, grepl("pool", sample_ID))     
        luma_pool <- filter (luma_pool, grepl("dup1", meth_dup))
    
    ## f) Plate by Plate Linear Regression
      # For each plate model the possible drift between controls; output is
      # a list of lm objects
        plate_drift = dlply(luma_pool, .(plate_rxn_ID), lm,
                            formula = methylation ~ plate_pos_seq)

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
        drift_coef <- ldply(plate_drift, extractfun)
   
    ## g) save the data frame of summary stats out as a pdf into output file
        pdf("output/drift_coef.pdf", height = 4, width = 8)
        grid.table(drift_coef)
        dev.off()
      

        
###########################################################
####                6  Build Data Set                  ####
###########################################################          
        
  ### 6.1 Calibrate Methylation Values
        
    ## a) Join drift_coef to luma_data
      # A Left join of 'luma_data' with 'drift_coef', making an updated
      # 'luma_data' dataframe which includes the drift slope, 'slope,'.
      # parent table. Parent tables are linked on 'plate_rxn_ID.'
        luma_data <- sqldf("SELECT
                            luma_data.*           
                            , drift_coef.slope   
                            FROM luma_data      
                            LEFT JOIN drift_coef       
                            ON luma_data.plate_rxn_ID = 
                            drift_coef.plate_rxn_ID")
        
    ## b) Weigthed Plate Calibration
      # Use ddply calculate a weighted plate calibration, with the result of 
      # shrinking drift towards the plate center (hy_pool control mean); 
      # a symmetrical shrinkage.
        calibration <- ddply (luma_data, .(plate_rxn_ID, sample_ID),
                           summarise,
                           meth_adjust = ifelse(plate_pos_factor == 1, 
                                      (((1-(plate_pos_seq/24))*slope) + 
                                         methylation), 
                                      (methylation - ((plate_pos_seq/24)-1)
                                       * slope))) 
      
    ## c) Join calibrated methylation to luma_data
      # A Left join of 'luma_data' with 'calibration', making an updated
      # 'luma_data' dataframe which includes the drift slope, 'slope,'.
      # parent table. Parent tables are linked on 'sample_ID.'
        luma_data <- sqldf("SELECT
                           luma_data.*           
                           , calibration .meth_adjust   
                           FROM luma_data      
                           LEFT JOIN calibration       
                           ON luma_data.sample_ID = 
                           calibration.sample_ID")    

                
  ### 6.2 Link Samples to Hyena ID
    # Join the luma_data back to samp_select so that methylation values can
    # be linked to hyena ID and kay code by the sample ID
    
    ## a) Join luma_data to samp_select
      # A Left join of 'luma_data' with 'samp_select', making an updated
      # 'luma_data' dataframe which includes the drift slope, 'slope,'.
      # parent table. Parent tables are linked on 'sample_ID.'
        luma_data <- sqldf("SELECT
                           luma_data. plate_rxn_ID, plate_pos_seq, 
                            plate_pos_factor, well, methylation,
                            analysis_status, assay_notes, dup1, dup2, stdev, cv,
                            slope, meth_adjust, stock_notes, cc_notes, 
                            rxn_notes,           
                            samp_select.*  
                           FROM luma_data      
                           LEFT JOIN samp_select       
                           ON luma_data.sample_ID = 
                           samp_select.sample_ID")   
    
    ## b) Reorder Age variable for graphing
        luma_data$Age <- factor(luma_data$Age, levels = c("cub", "subadult",
                                "adult"))
 
        
               
###########################################################
####      7  Descriptive Statistics (Univariate)       ####
###########################################################         

  ### 7.1 Outcome Univariate 
    ## a) Descriptive Stats Outcome
      # calculate the mean, median and standard deviation of % methylation
      # and the adjusted % methylation values
        univar_meth = ddply(luma_data, .(), summarise,
                            N = sum(!is.na(methylation)),
                            mean = mean(methylation, na.rm = T),
                            median =  quantile(methylation, c(.5),
                                                      na.rm = T),
                            sd = sd(methylation, na.rm = T),
                            N_adjust = sum(!is.na(meth_adjust)),
                            mean_adjust = mean(meth_adjust, 
                                                      na.rm = T),
                            median_adjust =  quantile(meth_adjust, c(.5),
                                                      na.rm = T),
                            sd_adjust = sd(meth_adjust, na.rm = T))
        
    ## b) save the data frame of summary stats out as a pdf into output file
        pdf("output/univar_meth.pdf", height = 4, width = 8)
        grid.table(univar_meth)
        dev.off()
        
    ## c) Histogram Outcome (adjust_meth)
        ggplot(data=luma_data, aes(x=meth_adjust, y = ..density..)) + 
          geom_histogram(breaks=seq(60, 100, by = 0.5), 
                         col="black",
                         aes(fill = ..count..)) +
          scale_fill_gradient("Count", low = "light green", high = 
                                "dark blue") +
          geom_density() +
          xlim(c(20,90)) +
          labs(title="Histogram for % Methylation") +
          labs(x="% Methylation", y="Frequency")
    
    ## d) Save Plot
      # use ggsave to save the linearization plot
        ggsave("meth_histogram.pdf", plot = last_plot(), device = NULL, 
               path = "./output", scale = 1, width = 7, height = 5, 
               units = c("in"), dpi = 300, limitsize = TRUE)
        
    ## e) Remove Outliers
      RemoveOutlier <- function (data, nos_sd, sd, mean) {
        low_cut <- mean - (nos_sd * sd) 
        hi_cut <- mean + (nos_sd * sd)
        data <- filter(data, methylation > low_cut & methylation < hi_cut)
      }  
    
    ## f) Run RemoverOutlier function to generate another data set      
      luma_data_no_out <- RemoveOutlier(data = luma_data, nos_sd = 2, 
                                        sd = univar_meth$sd,
                                        mean = univar_meth$mean)   

  
    ## g) Histogram Outcome (adjust_meth and outliers removed)
      ggplot(data=luma_data_no_out, aes(x=meth_adjust, y = ..density..)) + 
        geom_histogram(breaks=seq(60, 100, by = 0.5), 
                       col="black",
                       aes(fill = ..count..)) +
        scale_fill_gradient("Count", low = "light green", high = 
                              "dark blue") +
        geom_density() +
        xlim(c(60,90)) +
        labs(title="Histogram for % Methylation
             (No Outliers") +
        labs(x="% Methylation", y="Frequency")
      
    ## h) Save Plot
      # use ggsave to save the linearization plot
      ggsave("meth_histogram_no_out.pdf", plot = last_plot(), device = NULL, 
             path = "./output", scale = 1, width = 7, height = 5, 
             units = c("in"), dpi = 300, limitsize = TRUE)  
      
    
    ## a) Descriptive Stats Outcome
      # calculate the mean, median and standard deviation of % methylation
      # and the adjusted % methylation values
      univar_meth_no_out = ddply(luma_data_no_out, .(), summarise,
                                 N = sum(!is.na(methylation)),
                                 mean = mean(methylation, na.rm = T),
                                 median =  quantile(methylation, c(.5),
                                             na.rm = T),
                                 sd = sd(methylation, na.rm = T),
                                 N_adjust = sum(!is.na(meth_adjust)),
                                 mean_adjust = mean(meth_adjust, 
                                             na.rm = T),
                                 median_adjust =  quantile(meth_adjust, c(.5),
                                                    na.rm = T),
                                 sd_adjust = sd(meth_adjust, na.rm = T))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf("output/univar_meth_no_out.pdf", height = 4, width = 8)
      grid.table(univar_meth_no_out)
      dev.off()  
    

  ### 7.2 Predictive Variables Univariate      
    ## a) Descriptive Stats Sex
        univar_sex <- ddply(luma_data, .(Sex), summarise,
                            N = sum(!is.na(Sex)))
    
    ## b) save the data frame of summary stats out as a pdf into output file
        pdf("output/univar_sex.pdf", height = 4, width = 8)
        grid.table(univar_sex)
        dev.off() 
        
    ## c) Descriptive Stats Age
      # calculate the N, mean, median and standard deviation of hyena ages 
        univar_age <- ddply(luma_data, .(Age), summarise,
                            N = sum(!is.na(Age)),
                            mean = round(mean(AgeMonths, na.rm = T), 2),
                            median =  round (quantile(AgeMonths, c(.5),
                                               na.rm = T), 2),
                            sd = round(sd(AgeMonths, na.rm = T), 2))
       
    ## d) save the data frame of summary stats out as a pdf into output file
          pdf("output/univar_age.pdf", height = 5, width = 8)
          grid.table(univar_age)
          dev.off()
        
    

          
###########################################################
####            8  Identify Sample Re-Runs             ####
###########################################################         
          
  ### 8.1 Sample Re-Runs                                  
    ## a) Anti-join luma_data_no_out with samp_select
      # Identify Sample's from the orignial sample selection that need to be
      # re-run that failed QAQC or that were outliers
      re_runs <- anti_join(samp_select, luma_data_no_out, by = "sample_ID")
    
    ## b) Add notes columns to Re-Runs
        # A Left join of 'luma_data' with 'samp_select', making an updated
        # 'luma_data' dataframe which includes the drift slope, 'slope,'.
        # parent table. Parent tables are linked on 'sample_ID.'
          re_runs <- sqldf("SELECT
                             re_runs .sample_ID, new_box, cell, sample_type,
                             sample_state_quality, orig_est_volume_uL,
                             current_vol_uL, num_freeze_thaw_cyc, 
                             last_freeze_thaw_date, Hyena, ID, KayCode, 
                             DartingDate, Sex, Age,
                             luma_raw .plate_rxn_ID, plate_pos_seq, 
                             plate_pos_factor, well, methylation,
                             analysis_status, assay_notes, dup1, dup2, stdev, 
                             cv, stock_notes, cc_notes, rxn_notes
                             FROM re_runs      
                             LEFT JOIN luma_raw       
                             ON re_runs.sample_ID = 
                             luma_raw.sample_ID")    
        
        
        
        
#################################################
####       9. Save Intermediate Tables       ####
####             as Spreadsheets             ####
#################################################
          
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
        
        
        
        
        
        
        
        
            
