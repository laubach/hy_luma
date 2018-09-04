###############################################################################
##############       Spotted Hyena Global DNA Methylation        ##############
##############               LUMA Data Preparation               ##############
##############                 By: Zach Laubach                  ##############
##############             last updated: 27 Aug 2018             ##############
###############################################################################

### PURPOSE: This code is desingned to clean and prepare LUMA data for 
### downstream analyses.


  # Code Blocks
    # 1: Configure workspace
    # 2: Import data
    # 3: Data management
    # 4: Assess Controls
    # 5: Build Data Set
    # 6: Descriptive Statistics (Univariate)
    # 7: Save Intermediate Tables as .csv
    


###############################################################################
##############             1.  Configure workspace               ##############
###############################################################################

  ### 1.1 clear global environment
    rm(list = ls())


  ### 1.2 Install and load packages 
    ## a) Data Manipulation and Descriptive Stats Packages

      # Check for tidyverse and install if not already installed
        if (!'tidyverse' %in% installed.packages()[,1]){
          install.packages ('tidyverse')
        }
      # load tidyverse packages
        library ('tidyverse')
      
      # Check for sqldf and install if not already installed
        if (!'sqldf' %in% installed.packages()[,1]){
          install.packages ('sqldf')
        }
        options(gsubfn.engine = "R") #fixes tcltk bug; run before require sqldf
      # load tidyverse packages
        library ('sqldf')

      # Check for here and install if not already installed
        if (!'here' %in% installed.packages()[,1]){
          install.packages ('here')
        }
      # load here packages
        library ('here')
    
    ## b) Graph Plotting and Visualization Packages
        
      # Check for ggplot2 and install if not already installed
        if (!'ggplot2' %in% installed.packages()[,1]){
          install.packages ('ggplot2')
        }
      # load ggplot2 packages
        library ('ggplot2')
        
      # Check for gridExtra and install if not already installed
        if (!'gridExtra' %in% installed.packages()[,1]){
          install.packages ('gridExtra')
        }
      # load gridExtra packages
        library ('gridExtra')  
        
      # Check for sjPlot and install if not already installed
        if (!'sjPlot' %in% installed.packages()[,1]){
          install.packages ('sjPlot')
        }
      # load sjPlot packages
        library ('sjPlot')    
        
        
    ## b) Modeling Packages
        
        # Check for broom and install if not already installed
        if (!'broom' %in% installed.packages()[,1]){
          install.packages ('broom')
        }
        # loadbroom packages
        library ('broom')      

    
  ### 1.3 Get Version and Session Info
      R.Version()
      sessionInfo()
      
      # Developed in:   
      # R version 3.4.3 (2017-11-30)
      # Platform: x86_64-apple-darwin15.6.0 (64-bit)
      # Running under: macOS Sierra 10.12.4
      
      
  ### 1.4 Set working directory 
      setwd(here())
      
  
  ### 1.5 Set file paths for data importing and exporting
    ## a) Create path to LUMA data folder
    luma_data_path <- "~/R/R_wd/fisi/project/hy_GR_global_DNA_meth/LUMA/data/"
      
    ## b) Create path to meta data folder (compiled from repository query 
      # sample selection)
    meta_data_path <- paste("~/R/R_wd/fisi/project/hy_GR_global_DNA_meth/",
                              "LUMA/meta_data/", sep = '')
      
    ## c) Create path to a LUMA output folder
    luma_data_out_path <- paste("~/R/R_wd/fisi/project/hy_GR_global_DNA_meth/",
                                "LUMA/output/", sep = '')  

    
    
###############################################################################
##############                  2. Import data                   ##############
###############################################################################    
  
  ### 2.1 Import LUMA data files
    
    ## a) Create a list of raw LUMA file names 
      files <- list.files(paste(luma_data_path), pattern = "*.csv")
  
    ## b) Import LUMA files and bind them row by row into a single data frame
      luma_raw <- files %>%
        map (function(x) read_csv(file.path(luma_data_path, x))) %>% 
        reduce(rbind)

        
  ### 2.2 Import sample selection table
      # This file is a record of all samples selected for a project. Samples
      # are generated from querying the bio_repository and tblDarting.
    ## a) Read in sample record file
      samp_record <- read_csv(paste(meta_data_path,
                                    "LUMA_sample_records.csv", sep = ''))
      
    ## b) fix dates and times functions
      source (file = paste0("/Volumes/Holekamp/code_repository/R/",
                            "4_scripts_source/fix_dates_and_times.R"))
      
    ## c) Convert darting date to formatted date
      samp_record$darting.date <- fix.dates (samp_record$darting.date)
  
      
  ### 2.3 Import Access fisi backend
    # read in tidy Access fisi backend tables and save as data frames
    #  source(paste0("/Volumes/Holekamp/code_repository/R/1_output_tidy_tbls/",
     #              "load_tidy_tbls.R"))
      
    ## a) manually load tblHyenas 
      tblHyenas <- read.csv(paste0("/Volumes/Holekamp/code_repository/R",
                                   "/1_output_tidy_tbls/tblHyenas.csv"),
                            header = T)
    ## b) manually load tblDarting 
      tblDarting <- read_csv(paste0("/Volumes/Holekamp/code_repository/R",
                                   "/1_output_tidy_tbls/tblDarting.csv"))
      
      

            
  
###############################################################################
##############                 3. Data management                ##############
###############################################################################
    
  ### 3.1 QAQC: Remove failed wells
      # Update luma_raw by removing any well that failed due to high CV, 
      # non-specific peaks or low pyrosequencing signal (see meta data file
      # for specific details)
      luma_raw <- filter (luma_raw, analysis_status == "y")
        
        
  ### 3.2 Controls
    ## a) Subset Controls 
      # subset luma_raw into a data frame that contains all LUMA controls
        luma_cntrl <- filter (luma_raw, grepl("hy", sample_ID))
      
    ## b) Gather luma_cntrl
      # gather the luma_cntrl data so that methylation duplicates are in 
      # long format
      # NOTE: Run this if calculating intra-class CV based on each duplicate
#        luma_cntrl <- luma_cntrl %>%
#          gather (meth_dup, dup, dup1:dup2)
        
    
  ### 3.3 Linearization Standards
    ## a) Subset Linearization Standards
      # subset luma_raw into a data frame that contains all LUMA linearization 
      # standards
        luma_linearz <- filter (luma_raw, grepl("lam", sample_ID)) 
    
    ## b) Subset luma_linearz to include only one set of lambda linearizations
        luma_linearz <- filter (luma_linearz, grepl("%", sample_ID)) 
        
      
  ### 3.4 Data
    ## a)  Subset Data
      # subset luma_raw into a data frame that contains all LUMA data. 
      # sample_IDs are numbers
        luma_data <- filter (luma_raw, grepl("^[[:digit:]]", sample_ID)) 
 
        
        
###############################################################################
##############                 4. Assess controls                ##############
###############################################################################
        
  ### 4.1 Intra-Class CV
    ## a) hy_pool and hy_100% Intra-Class CV
      # use dplyr to calulate the intra-class CV (within plate variation)
      # for control samples.
        intra_CV <- luma_cntrl %>%
          group_by(plate_rxn_ID, sample_ID) %>%
          summarize (n = n(),
                     avg = round (mean (methylation, na.rm = T), 2),
                     median =  round (quantile (methylation, c(.5), na.rm = T),
                                      2),
                     stdev = round (sd (methylation, na.rm = T), 2),
                     cv = round (100 * (stdev/avg), 3))
        
    ## b) save the data frame of summary stats out as a pdf into output file
        pdf(paste0 (here(),"/output/output_luma_prep/intra_CV.pdf"), 
                   height = 9, width = 7)
        grid.table(intra_CV)
        dev.off()
        
        
  ### 4.2 Inter-Class CV
    ## a) hy_pool and hy_100% Inter-Class CV
      # use ddply to calulate the inter-class CV (between plate variation) 
      # for control samples. 
        inter_CV <- luma_cntrl %>%
          group_by(sample_ID) %>%
          summarize (n = n(),
                     avg = round (mean (methylation, na.rm = T), 2),
                     median =  round (quantile (methylation, c(.5), na.rm = T),
                                      2),
                     stdev = round (sd (methylation, na.rm = T), 2),
                     cv = round (100 * (stdev/avg), 3))
        
    ## b) save the data frame of summary stats out as a pdf into output file
        pdf(paste0(here(),"/output/output_luma_prep/inter_CV.pdf"), 
            height = 3, width = 7)
        grid.table(inter_CV)
        dev.off()    

              
  ### 4.3 Check Linearization
    ## a) Organize data and create vector of predicted Values
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
               path = paste0(here(), "/output/output_luma_prep"), 
                             scale = 1, width = 6, height = 4, 
               units = c("in"), dpi = 300, limitsize = TRUE)

        
  ### 4.4 Check Plates for Drift    
        
    ## a) Graph Controls Across Plates 
      # Graph the controls on each reaction plate to assess for any plate drift
      # least squares regression is used for the fit function 
        ggplot(luma_cntrl, aes (x = plate_pos_seq, y = methylation,
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
        ggsave("double_panel_control_drift.pdf", plot = last_plot(), 
               device = NULL, 
               path = paste0(here(), "/output/output_luma_prep"),
               scale = 1, width = 8, height = 11, 
               units = c("in"), dpi = 300, limitsize = TRUE) 
   
    ## c) Subset Controls 
      # subset luma_raw into a data frame that contains one of the LUMA control
      # duplicates
        luma_pool <- filter (luma_cntrl, grepl("pool", sample_ID))     
        #luma_pool <- filter (luma_pool, grepl("dup1", meth_dup))
    
    ## d) Plate by Plate Linear Regression
      # For each plate model the possible drift between controls; output is
      # a list of lm (linear models) objects
        plate_drift <- luma_pool %>%
          group_by(plate_rxn_ID) %>%
          do(lm_mods = lm(methylation ~ plate_pos_seq, data = .))
      
      # use 'broom' to extract model coefficients from lists as a data frame      
        drift_coef <- tidy(plate_drift, lm_mods)
        
      # filter to remove intercepts and retain only slope coefs and test stats
        drift_coef <- filter (drift_coef, grepl("plate", term)) 
       
    ## e) save the data frame of summary stats out as a pdf into output file
        pdf(paste0 (here(),"/output/output_luma_prep/drift_coef.pdf"), 
            height = 11, width = 8)
        grid.table(drift_coef)
        dev.off()

        

###############################################################################
##############                 5. Build Data Set                 ##############
###############################################################################       
         
  ### 5.1 Calibrate Methylation Values
        
    ## a) Join drift_coef to luma_data
      # A Left join of 'luma_data' with 'drift_coef', making an updated
      # 'luma_data' dataframe which includes the drift slope, 'slope,'.
      # parent table. Parent tables are linked on 'plate_rxn_ID.'
        luma_data <- sqldf("SELECT
                            luma_data.*           
                            , drift_coef.estimate   
                            FROM luma_data      
                            LEFT JOIN drift_coef       
                            ON luma_data.plate_rxn_ID = 
                            drift_coef.plate_rxn_ID")
        
    ## b) Weigthed Plate Calibrations
      # Use dplyr to calculate a weighted plate calibration, with the result of 
      # shrinking drift towards the plate center (hy_pool control mean); 
      # a symmetrical shrinkage. 
<<<<<<< HEAD
        # NOTE: samples added here with single channel pipette well by well
        # so plate_pos_seq is count from 1-48
        calibration_sing <- luma_data %>%
          filter(plate_rxn_ID == "p1r1" | plate_rxn_ID == "p1r2" |
                   plate_rxn_ID == "p2r3" |  plate_rxn_ID == "p2r4" |
                   plate_rxn_ID == "p3r5" | plate_rxn_ID == "p3r6" |
                   plate_rxn_ID == "p4r7" | plate_rxn_ID == "p4r8" |
                   plate_rxn_ID == "p7r13_15" | plate_rxn_ID == "p6.2r14") %>%
          group_by (plate_rxn_ID, sample_ID) %>%
          summarize(meth_adjust = ifelse(plate_pos_factor == 1,
                                         (((1-(plate_pos_seq/30))*estimate) +
                                            methylation),
                                         (methylation - ((plate_pos_seq/30)-1)
                                          * estimate))) 
        
        # NOTE: samples added here with single and multi channel pipette 
        # so plate_pos_seq is count from 1-17
        calibration_mult1 <- luma_data %>%
          filter(plate_rxn_ID == "p5r9") %>%
          group_by (plate_rxn_ID, sample_ID) %>%
          summarize(meth_adjust = ifelse(plate_pos_factor == 1,
                                         (((1-(plate_pos_seq/14))*estimate) +
                                            methylation),
                                         (methylation - ((plate_pos_seq/14)-1)
                                          * estimate))) 
        
        # NOTE: samples added here with multi channel pipette 
        # so plate_pos_seq is count from 1-9
        calibration_mult2 <- luma_data %>%
          filter(plate_rxn_ID == "p5r10" | plate_rxn_ID == "p6r11" |
                   plate_rxn_ID == "p6r12") %>%
=======
        # NOTE: samples added here with single channel pippet well by well
        # so plate_pos_seq is count from 1-48
        calibration <- luma_data %>%
          #******************************************filter by here
>>>>>>> bcc99c27b6282d228a325396aaae17b00cac78d0
          group_by (plate_rxn_ID, sample_ID) %>%
          summarize(meth_adjust = ifelse(plate_pos_factor == 1,
                                         (((1-(plate_pos_seq/6))*estimate) +
                                            methylation),
                                         (methylation - ((plate_pos_seq/6)-1)
                                          * estimate))) 
<<<<<<< HEAD
        
    ## c) Combine calibration data frames (row by row)
        calibration <- rbind(calibration_sing, calibration_mult1, 
                             calibration_mult2)
=======
        # NOTE: samples added here with single channel pippet well by well
        # so plate_pos_seq is count from 1-48
        calibration2 <- luma_data %>%
          group_by (plate_rxn_ID, sample_ID) %>%
          summarize(meth_adjust = ifelse(plate_pos_factor == 1,
                                         (((1-(plate_pos_seq/2))*estimate) +
                                            methylation),
                                         (methylation - ((plate_pos_seq/2)-1)
                                          * estimate))) 
>>>>>>> bcc99c27b6282d228a325396aaae17b00cac78d0
      
    ## d) Join calibrated methylation to luma_data
      # A Left join of 'luma_data' with 'calibration', making an updated
      # 'luma_data' dataframe which includes the drift slope, 'slope,'.
      # parent table. Parent tables are linked on 'sample_ID.'
        luma_data <- sqldf("SELECT
                           luma_data.*           
                           , calibration .meth_adjust   
                           FROM luma_data      
                           LEFT JOIN calibration       
                           ON luma_data.sample_ID = 
                           calibration.sample_ID AND
                           luma_data.plate_rxn_ID = 
                           calibration.plate_rxn_ID")    

        
  ### 5.1 Need to remove the duplicates that were run twice on accident
        # see rxn notes
    ## a) make a list of duplicates    
        luma_data %>%
          filter(duplicated(.[["sample_ID"]]))
    ## b) Group rows with same sample_id to reduce duplicate
        luma_data <- luma_data %>% 
          group_by (sample_ID) %>% # set grouping same ID within same cat age
          summarise (plate_rxn_ID = first(plate_rxn_ID),
                     plate_pos_seq = first(plate_pos_seq),
                     plate_pos_factor = first(plate_pos_factor),
                     well = first(well),
                     methylation = mean(methylation),
                     analysis_status = first(analysis_status),
                     assay_notes = first(assay_notes),
                     dup1 = mean(dup1),
                     dup2 = mean(dup2),
                     stdev = mean(stdev),
                     cv = mean(cv),
                     stock_notes = first(stock_notes),
                     cc_notes = first(cc_notes),
                     rxn_notes = first(rxn_notes),
                     drift_est = mean(estimate), # update variable name here
                     meth_adjust = mean(meth_adjust))

        
  ### 5.2 Link Samples to Hyena ID
    # Join the luma_data back to samp_record so that methylation values can
    # be linked to hyena ID and kay code by the sample ID
    
    ## a) Join luma_data to samp_record
      # A Left join of 'luma_data' with 'samp_record', making an updated
      # 'luma_data' dataframe. Parent tables are linked on 'sample_ID.'
        luma_data <- sqldf("SELECT
                           luma_data. plate_rxn_ID, plate_pos_seq, 
                            plate_pos_factor, well, methylation,
                            analysis_status, assay_notes, dup1, dup2, stdev, cv,
                            drift_est, meth_adjust, stock_notes, cc_notes, 
                            rxn_notes,           
                            samp_record.*  
                           FROM luma_data      
                           LEFT JOIN samp_record       
                           ON luma_data.sample_ID = 
                           samp_record.sample_ID
                           GROUP BY plate_rxn_ID, well") # to get rid of dups 
                                                         # created in join
        
        
    ## b) Convert darting date to formatted date in luma_data 
        luma_data$darting.date <- fix.dates (luma_data$darting.date)
    
        
    ## c) Join tblHyenas to tblDarting
        # A Left join of 'tblDarting' with select columns from 'tblHyenas', 
        # making a new dataframe, 'tbl_dart_hy'. 
<<<<<<< HEAD
        # Parent tables are linked on 'id'. 
        tblDarting <- tblDarting %>%
          select(-c(clan))
      
#************************ fix table darting hack *****************************
# replaces full name with abbreviated name
        tblDarting$id <- ifelse((nchar(tblDarting$id) >= 4 & 
                                   nchar(tblDarting$hyena) <= 4), 
                                 tblDarting$hyena,
                                 tblDarting$id)
        
      tbl_dart_hy <- tblDarting %>%
        left_join(select(tblHyenas, c(id, first.seen, den.grad, disappeared, 
                                      mom, birthdate, number.littermates,
                                      litrank, mortality.source, death.date,
                                      weaned, clan, park)), by = "id")
       
    ## d) convert darting.date to Date format
        tbl_dart_hy$darting.date <- as.Date(tbl_dart_hy$darting.date)
        
    ## e) Join luma_data to tbl_dart_hy
=======
        # Parent tables are linked on 'id'.   
        tbl_dart_hy <- tblDarting %>%
          select(-one_of("clan")) %>%
          left_join(select(tblHyenas, c(id, first.seen, den.grad, disappeared, 
                                       mom, birthdate, number.littermates,
                                       litrank, mortality.source, death.date,
                                       weaned, clan, park)), by = "id")
        
    ## d) Join luma_data to tbl_dart_hy
>>>>>>> bcc99c27b6282d228a325396aaae17b00cac78d0
        # A Left join of 'luma_data' with 'tbl_dart_hy', making an updated
        # 'luma_data' dataframe. Parent tables are linked on 'kay.code' and 
        # 'darting.date'  
        luma_data <- luma_data %>%
          left_join(select(tbl_dart_hy, -one_of("id")),
            by = c("kay.code", "darting.date"))

    ## f) Reorder Age variable for graphing
        luma_data$age <- factor(luma_data$Age, levels = c("cub", "subadult",
                                "adult"))
    
    ## g) Manual data clean up
        # cash is in tblHyenas twice so need to remove erroroneous entry
        luma_data <- luma_data %>%
          filter (!grepl("cash", id) | !grepl("14-apr-08", first.seen))
   
    
  
 
###############################################################################
##############       6. Descriptive Statistics (Univariate)      ##############
###############################################################################               
     

  ### 6.1 Outcome Univariate 
    ## a) Descriptive Stats Outcome
      # calculate the mean, median and standard deviation of % methylation
      # and the adjusted % methylation values
        
        univar_meth <- luma_data %>%
          #group_by(plate_rxn_ID, sample_ID) %>%
          summarize (n = n(),
                     avg = round (mean (methylation, na.rm = T), 2),
                     median =  round (quantile (methylation, c(.5), na.rm = T),
                                      2),
                     stdev = round (sd (methylation, na.rm = T), 2),
                     n_adjust = sum(!is.na(meth_adjust)),
                     avg_adjust = round (mean (meth_adjust, 
                                               na.rm = T), 2),
                     median_adjust =  round (quantile (meth_adjust, 
                                                       c(.5), na.rm = T), 2),
                     stdev = round (sd (meth_adjust, 
                                        na.rm = T), 2))
                    

    ## b) save the data frame of summary stats out as a pdf into output file
        pdf(paste0(here(),"/output/output_luma_prep/univar_meth.pdf"),
            height = 4, width = 8)
        grid.table(univar_meth)
        dev.off()
        
    ## c) Histogram Outcome (methylation)
        ggplot(data=luma_data, aes(x=methylation, y = ..density..)) + 
          geom_histogram(breaks=seq(60, 100, by = 0.5), 
                         col="black",
                         aes(fill = ..count..)) +
          scale_fill_gradient("Count", low = "light green", high = 
                                "dark blue") +
          geom_density() +
          xlim(c(20,90)) +
          labs(title="Histogram of % Methylation") +
          labs(x="% Methylation", y="Frequency")
    
    ## d) Save Plot
      # use ggsave to save the linearization plot
        ggsave("meth_histogram.pdf", plot = last_plot(), device = NULL, 
               path = paste0(here(), "/output/output_luma_prep"),
               scale = 1, width = 7, height = 5, 
               units = c("in"), dpi = 300, limitsize = TRUE)
        
    ## e) Remove Outliers
      RemoveOutlier <- function (data, nos_sd, sd, mean) {
        low_cut <- mean - (nos_sd * sd) 
        hi_cut <- mean + (nos_sd * sd)
        data <- filter(data, methylation > low_cut & methylation < hi_cut)
      }  
    
    ## f) Run RemoverOutlier function to generate another data set      
      luma_data_no_out <- RemoveOutlier(data = luma_data, nos_sd = 2, 
                                        sd = univar_meth$stdev,
                                        mean = univar_meth$avg)   

  
    ## g) Histogram Outcome (adjust_meth and outliers removed)
      ggplot(data=luma_data_no_out, aes(x=meth_adjust, y = ..density..)) + 
        geom_histogram(breaks=seq(60, 100, by = 0.5), 
                       col="black",
                       aes(fill = ..count..)) +
        scale_fill_gradient("Count", low = "light green", high = 
                              "dark blue") +
        geom_density() +
        xlim(c(60,90)) +
        labs(title="Histogram of % Methylation
             (No Outliers") +
        labs(x="% Methylation", y="Frequency")
      
    ## h) Save Plot
      # use ggsave to save the linearization plot
      ggsave("meth_histogram_no_out.pdf", plot = last_plot(), device = NULL, 
             path =  paste0(here(), "/output/output_luma_prep"),
             scale = 1, width = 7, height = 5, 
             units = c("in"), dpi = 300, limitsize = TRUE)  
      
    
    ## i) Descriptive Stats Outcome
      # calculate the mean, median and standard deviation of % methylation
      # and the adjusted % methylation values
      univar_meth_no_out <- luma_data_no_out %>%
        summarize (n = n(),
                   avg = round (mean (methylation, na.rm = T), 2),
                   median =  round (quantile (methylation, c(.5), na.rm = T),
                                    2),
                   stdev = round (sd (methylation, na.rm = T), 2),
                   n_adjust = sum(!is.na(meth_adjust)),
                   avg_adjust = round (mean (meth_adjust, 
                                             na.rm = T), 2),
                   median_adjust =  round (quantile (meth_adjust, 
                                                     c(.5), na.rm = T), 2),
                   stdev = round (sd (meth_adjust, 
                                      na.rm = T), 2))
      
      
    ## j) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/output_luma_prep/univar_meth_no_out.pdf"),
          height = 4, width = 8)
      grid.table(univar_meth_no_out)
      dev.off()
      

  ### 6.2 Identify Sample Re-Runs                    
                                 
    ## a) Anti-join luma_data_no_out with samp_select
      # Identify Sample's from the orignial sample selection that need to be
      # re-run that failed QAQC or that were outliers
      re_runs <- anti_join(samp_record, luma_data_no_out, by = "sample_id")
        
        
        
###############################################################################
##############        7. Save Intermediate Tables as .csv        ##############
###############################################################################         

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
      csv.file.name.luma <- paste (luma_data_out_path, "luma_data_no_out",
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
        
        
        
        
        
        
        
        
            
