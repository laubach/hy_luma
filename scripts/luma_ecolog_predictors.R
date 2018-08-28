###############################################################################
##############       Spotted Hyena Global DNA Methylation        ##############
##############             LUMA Ecological Predictors            ##############
##############                 By: Zach Laubach                  ##############
##############             last updated: 10 Aug 2018              ##############
###############################################################################


### PURPOSE: This code is desingned to analyze LUMA data and associations of 
          # social and ecological predictors of DNA methylation in hyenas.


  # Code Blocks
    # 1: Configure workspace
    # 2: Import data
    # 3: Data management
    # 4: Univariate analyses
    # 5: Data transformations 
    # 6: Bi-variate analyses 
    # 7: Cub models
    # 8: Subadult models
    # 9: Adult models



###############################################################################
##############             1.  Configure workspace               ##############
###############################################################################

  ### 1.1 clear global environment
    rm(list = ls())
  

  ### 1.2 Install and load packages 
    ## a) Data Manipulation and Descriptive Stats Packages
    
      # Check for tidyverse and install if not already installed
        if (!'tidyverse' %in% installed.packages()[,1]){
          install.packages ('tidyYverse')
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
      
      # Check for lubridate and install if not already installed
        if (!'lubridate' %in% installed.packages()[,1]){
          install.packages ('lubridate')
        }
      # load lubridate packages
        library ('lubridate') 
        
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
    
    ## c) Modeling Packages
      # Check for broom and install if not already installed
        if (!'broom' %in% installed.packages()[,1]){
          install.packages ('broom')
        }
      # load broom packages
        library ('broom')   
        
      # Check for lmerTest and install if not already installed
      #  if (!'lmerTest' %in% installed.packages()[,1]){
      #    install.packages ('lmerTest')
      #  }
      # load lmerTest packages
      #  library ('lmerTest')   
        
        
      # Check for nlme and install if not already installed
        if (!'nlme' %in% installed.packages()[,1]){
          install.packages ('nlme')
        }
      # load nlme packages
        library ('nlme')
      
      # Check for lme4 and install if not already installed
        if (!'lme4' %in% installed.packages()[,1]){
          install.packages ('lme4')
        }
      # load lme4 packages
        library ('lme4')
        
      # Check for effects and install if not already installed
        if (!'effects' %in% installed.packages()[,1]){
          install.packages ('effects')
        }
        # load effects packages
        library ('effects')
        
      # Check for aod and install if not already installed
        if (!'aod' %in% installed.packages()[,1]){
          install.packages ('aod')
        }
      # load aod packages
        library ('aod')
    
  
  ### 1.3 Get Version and Session Info
    R.Version()
    sessionInfo()
    
    # Developed in:   
      # R version 3.5.1 (2018-07-02)
      # Platform: x86_64-apple-darwin15.6.0 (64-bit)
      # Running under: macOS High Sierra 10.13.6
    
    
  ### 1.4 Set working directory 
    setwd(here())

  
  ### 1.5 Set file paths for data importing and exporting
    ## a) The path to cleand LUMA data
      luma_data_path <- paste("~/R/R_wd/fisi/project/hy_GR_global_DNA_meth/",
                              "LUMA/output/", sep = '')
    
    ## b) The path to prey data
      prey_data_path <- paste("~/Git/fisi_lab/hy_prey_density/output/",
                              sep = '')
      
    ## c) The path to maternal care data
      fas_data_path <- paste("/Volumes/Holekamp/CurrentCollaborations/",
                             "maternal_care/R_code/fas/data/",
                             "final_data_with_variables/",
                             sep = '')
    

    
###############################################################################
##############                  2. Import data                   ##############
###############################################################################    
      
  ### 2.1 Import LUMA data files

    ## a) Import LUMA data, which has undergone QAQC (see R script 
      # luma_prep_analysis.R)
      luma_data <- read_csv(paste(luma_data_path,
                                  "luma_data_no_out.csv", sep = ''))
      
    ## b) check that all variables are of appropriate class    
      sapply(luma_data, class)
      
  ### 2.2 Import prey density data       
    # Import prey density data, which was created with the R script
    # 'calc_prey_density'
     
    ## a) Create a list of prey summary file names 
      files <- list.files(paste(prey_data_path), pattern = "*.csv")
      
    ## b) Import prey density files and row bind into a single data frame
      prey_density <- files %>%
        map (function(x) read_csv(file.path(prey_data_path, x))) %>% 
        reduce(rbind)
     
       
  ### 2.3 Import Access fisi backend #TEMPORARILY NOT WORKING
    ## a) read in tidy Access fisi backend tables and save as data frames
    #  source(paste0("/Volumes/Holekamp/code_repository/R/1_output_tidy_tbls/",
    #               "load_tidy_tbls.R"))
      
    ## b) manually load tblFemalerank  
    tblFemalerank <- read_csv(paste0("/Volumes/Holekamp/code_repository/R",
                                    "/1_output_tidy_tbls/tblFemalerank.csv"))
    

    
###############################################################################
##############                 3. Data management                ##############
###############################################################################
  
    
  ### 3.1 Tidy luma_data  
    ## a) view the variable classes to determine which need to be modified
      spec(luma_data)
    
    ## b) fix dates and times functions TEMPORARILY NOT WORKING
    #source (file = paste0("/Volumes/Holekamp/code_repository/R/",
      #                    "4_scripts_source/fix_dates_and_times.R"))
    
    ## c) Convert dates stored as character (e.g. 03-aug-05) to formatted dates
      luma_data$darting.date <- as.Date(luma_data$darting.date, 
                                      format = "%d-%b-%y")
      luma_data$birthdate <- as.Date(luma_data$birthdate, 
                                      format = "%d-%b-%y")
      luma_data$weaned <-  as.Date(luma_data$weaned, 
                                      format = "%d-%b-%y")
      luma_data$first.seen <-  as.Date(luma_data$first.seen, 
                                      format = "%d-%b-%y")
      luma_data$den.grad <-  as.Date(luma_data$den.grad, 
                                      format = "%d-%b-%y")
      luma_data$disappeared <-  as.Date(luma_data$disappeared, 
                                      format = "%d-%b-%y")
      luma_data$death.date <-  as.Date(luma_data$death.date, 
                                        format = "%d-%b-%y")
      
    
    ## d) Create an estimated age in months by subracting birthdate from
      # darting.date using lubridate
      luma_data <- luma_data %>%
        mutate(age.mon = interval(birthdate, darting.date) %/% months(1))
      
    ## e) Cobmine age columns to fill in NA
      # replace NA in age.months column where there is a value in 
      # estimated.age.mo column
      luma_data$age.months <- ifelse(is.na(luma_data$age.months),
                                     luma_data$estimated.age.mo,
                                     luma_data$age.months)
      # drop extra EstimatedAgeMo column
      luma_data <- luma_data %>%
        select (- c(estimated.age.mo)) 
      
      # replace NA in age.mon column where there is a value in 
      # age.months column
      luma_data$age.mon <- ifelse(is.na(luma_data$age.mon),
                                   luma_data$age.months,
                                   luma_data$age.mon)
      
    ## f) fix age transcription inconsistencies (e.g. all 'sub' to 'subadult') 
      luma_data <- luma_data %>%
        mutate(age = replace(age, str_detect(age, "s"), "subadult")) 
      
    ## g) Create a categorical age variable based age.mon 
      # darting.date using lubridate
      luma_data <- luma_data %>%
        mutate(age.cat = case_when(age.mon <= 12 ~ c("cub"),
                                          age.mon > 12 & 
                                            age.mon <=24 ~ c("subadult"),
                                          age.mon > 24 ~ c("adult")))
      
    ## h) replace NA in age.cat column where there is a value in 
      # age column
      # first convert age to character
      luma_data$age <- as.character(luma_data$age)
      luma_data$age.cat <- ifelse(is.na(luma_data$age.cat),
                                         luma_data$age,
                                         luma_data$age.cat)
      
      # drop extra EstimatedAgeMo column
      luma_data <- luma_data %>%
        select (- c(age)) 
      
    ## f) Create *nominal* factor (with ordered levels)  
      # Set levels (odering) age.cat variable and sets the reference level 
      # to cub makes this
      # NOTE: model output is difference in means btwn reference and each 
      # level; factorial contrasts are differences in group level means
      # stored internally as 1, 2, 3 (equal spacing)
      luma_data <- transform(luma_data, 
                             age.cat = factor(age.cat,
                                              levels = c("cub", 
                                                         "subadult", 
                                                         "adult")))
    ## e) Create *interval* factor (with ordered levels)  
      # Set order and *equally* spaced categorical age variable a
      # NOTE: model output a p for trend
      # polynomial contrasts are linear, quadratic, cubic etc. estimates
      # stored internally as 1, 2, 3 etc. (equal spacing)
      luma_data <- transform(luma_data, 
                             age.inter = ordered(age.cat,
                                                 levels = c("cub",
                                                            "subadult",
                                                            "adult")))
     
    ## g) Create *ordinal* factor (with ordered levels)  
      # Set order and *unequally* spaced categorical age variable 
      # (e.g. median/mean of each category)
      # NOTE: model output a p for trend; assumes linear monotonic association
      # estimates, which is a change y with respect to each unit change in x
      # is scaled to x and stored internally as numeric value
      luma_data <- luma_data %>%
        mutate(age.ordin = as.numeric(case_when(age.mon <= 12 ~ c(10.65),
                                   age.mon > 12 & 
                                     age.mon <=24 ~ c(17.49),
                                   age.mon > 24 ~ c(60.34))))

    ## h) Re-code sex as a nominal factor  
      luma_data$sex <- as.factor(luma_data$sex)

    ## i) extract the year for date of interest (here Birthdate) using lubridate
      # and make a new variable
      luma_data$rank_year <- year(as.Date(luma_data$birthdate, 
                                          format="%Y-%m-%d"))
      
    ## j) Make a new variable hum.distrb
      # create a 3-level ordinal factor indicating human population size 
      # (hi - talek hyena born after 2000, med - fig tree born after 2000,
      # and lo - all other hynenas for which clan is known and birthdate)
      luma_data <- luma_data  %>%
        mutate(hum.disturb = case_when(luma_data$clan %in% c("talek") &
                                         luma_data$rank_year > 2000 ~ c("hi"),
                                       (luma_data$clan %in% c("fig.tree") &
                                          luma_data$rank_year > 2000 ~ 
                                          c("med")),
                                       luma_data$clan %in% c("talek", 
                                                             "fig.tree") &
                                         luma_data$rank_year <= 2000 |
                                         luma_data$clan %in% c("serena","happy", 
                                                               "mara")
                                       ~ c("low")))
                                    
                                    
    ## j) Re-order hum.distrb levels and re-code variable as categorical
      luma_data <- transform(luma_data,
                             hum.distrb = factor(hum.distrb,
                                              levels = c("low", "med", "hi")))
      

  
    ## k) Make a new variable intra.lit.rank
      #create blank variable to store data 
    #   luma_data$intra.lit.rank <- NA
    
    ## l) Loop through number.littermates and litrank variables to create
      # a three-level factor 
    #  for (i in 1:nrow(luma_data)) {
    #    if((!is.na (luma_data$number.littermates[i]))  # check no NAs
    #       && luma_data$number.littermates[i] == 0) { # check 0 littermates
    #      luma_data$intra.lit.rank[i] <- "single"}
    #    else if((!is.na (luma_data$litrank[i]))  # check no NAs
    #            && luma_data$litrank[i] == 1) { # check twin's rank 1
    #      luma_data$intra.lit.rank[i] <- "twin_hi"}
    #    else if((!is.na (luma_data$litrank[i]))  # check no NAs
    #            && luma_data$litrank[i] == 2) { # check twin's rank 2
    #      luma_data$intra.lit.rank[i] <- "twin_lo"}
    #  }

    ## k) Make a new 2 level nomial variable lit.size
      luma_data <- luma_data  %>%
        mutate(lit.size = case_when(number.littermates == 0 ~ c("single"),
                                          litrank == 1 ~ c("twin"),
                                          litrank == 2 ~ c("twin")))

    ## l) Recod as nominal factor and set level (order)
      luma_data <- transform(luma_data,
                             lit.size = factor(lit.size, 
                                               levels = c("single","twin")))

   
  ### 3.2 Tidy tblFemalerank
    ## a) Pattern recognize numbers from Year variable and copy; gets rid of
      # unwanted text characters
      tblFemalerank$rank_year <- as.numeric(regmatches(tblFemalerank$year, 
                                           gregexpr("[[:digit:]]+",
                                                    tblFemalerank$year)))
      
    ## b) rename 'id' variable as 'mom' 
      tblFemalerank <- rename_(tblFemalerank, "mom" = "id")
      ### 3.3 Join Data Sets: Append female rank covariates to LUMA data
      
      
      
  ### 3.3 Left join tblFemalerank to luma_data  
      # append to luma_data each hyena's mom's rank from the year that 
      # they were born 
      luma_data2 <- sqldf("SELECT
                         luma_data.*           
                         , tblFemalerank. absrank, strank
                         FROM luma_data      
                         LEFT JOIN tblFemalerank       
                         ON tblFemalerank.mom = luma_data.mom
                         AND tblFemalerank.rank_year = luma_data.rank_year") 
      
  
  ### 3.4 Tidy joined table 
    ## a) rename 'absrank' variable as 'mom.absrank' 
      luma_data <- rename_(luma_data, "mom.absrank" = "absrank")
      
    ## b) rename 'strank' variable as 'mom.strank' 
      luma_data <- rename_(luma_data, "mom.strank" = "strank") 
      
    ## c) convert mom's strank to numeric
      luma_data$mom.strank <- as.numeric(luma_data$mom.strank)
      
    ## f) Create quartiles of maternal rank
      # use quantile function to cut strank into 4 levels (approximately same
      # number in each level)
      luma_data <- within(luma_data, mom.strank.quart.order
                          <- as.integer(cut(mom.strank, 
                                            quantile(mom.strank, probs=0:4/4,
                                                     na.rm = T), 
                                            include.lowest = T)))
      # NOTE: Could also use dplyr
      #  luma_data <- luma_data %>% 
      #    mutate(strank.quart.order = ntile(strank, 4))
      

    ## g) Create a nominal factor and rename and re-order the levels to 
      # sets the reference level for lowest rank 
      luma_data <- transform (luma_data, 
                              mom.strank.quart = 
                                factor(mom.strank.quart.order,
                                       levels = c(1, 2, 3, 4),
                                       labels= c("Q1 (lowest)",
                                                 "Q2","Q3",
                                                 "Q4 (highest)")))
  
  ### 3.5 Reduce and group luma_data
    # Create new reduced and grouped data sets to be used in some downstream
    # anlayses

    ## a) Group rows with same hyena ID and within an age group
        luma_data_group <- luma_data %>% 
          filter (!is.na(methylation))%>% # check/remove rows where meth NA
          group_by (id, age.cat) %>% # set grouping same ID within same cat age
          summarise (age.reps = sum(!is.na(methylation)), # n per ID w/in age 
                                                          # class
                     methylation = mean(methylation), # avg methylation
                                                      # value for repeat ID 
                                                      # within same age cat
                     meth_adjust = mean(meth_adjust), # avg adjust meth
                     kay.code = first(kay.code),
                     darting.date = as.Date(mean(darting.date)), # for multiple
                      # darting dates in same age category, take avg of dates
                     sex = first(sex),
                     #age = first(age),
                     age.months = mean(age.months),   # avg age in months
                     reproductive.status = first(reproductive.status),
                     mass = mean(mass),
                     glucose.green = mean(glucose.green),
                     glucose.red = mean(glucose.red),
                     pcv = mean(pcv),
                     wbc = mean(wbc),
                     rbc = mean(rbc),
                     total.solids = mean(total.solids),
                     first.seen = as.Date(first(first.seen)),
                     den.grad = as.Date(first(den.grad)),
                     disappeared = as.Date(first(disappeared)),
                     mom = first(mom),
                     birthdate = as.Date(first(birthdate)),
                     number.littermates = first(number.littermates),
                     litrank = first(litrank),
                     mortality.source = first(mortality.source),
                     death.date = as.Date(first(death.date)),
                     weaned = as.Date(first(weaned)),
                     clan = as.factor(first(clan)),
                     park = as.factor(first(park)),
                     mom.absrank = as.factor(first(mom.absrank)),
                     mom.strank = as.numeric(first(mom.strank)),
                     mom.strank.quart = as.factor(first(mom.strank.quart)),
                     mom.strank.quart.order = (first(mom.strank.quart.order)),
                     hum.distrb = (first(hum.distrb)),
                     lit.size = (first(lit.size)),
                     age.mon = (mean(age.mon)),
                     rank_year = (first(rank_year)),
                     age.ordin = first(age.ordin))
                     
    ## b) Ungroup the data frame
        # dplyr retains grouping after creation of data frame that uses 
        # group_by
          luma_data_group <- ungroup(luma_data_group)
        
    
  ### 3.6 Clean prey density data and combine with LUMA data
    ## a) Select a subset of the prey_density data by column names
        prey_density <- prey_density %>%
          select(grep("tot", names(prey_density)), # contains 'tot'
                 grep("num", names(prey_density)), # contains 'num'
                 grep("^id$", names(prey_density))) %>% # exact match 'ID'
          select (-c(number.littermates)) %>%
          rename ("total.birth.3" = "total.birth-3") %>%
          rename ("total.3.6" = "total.3-6") %>%
          rename ("total.6.9" = "total.6-9")
          
    ## b) Left join prey_density to luma_data   
        luma_data <- left_join(luma_data,
                               prey_density, by = "id")
        
    ## c) Left join prey_density to luma_data_group   
        luma_data_group <- left_join(luma_data_group,
                              prey_density, by = "id")
  
        
             
###############################################################################
##############               4. UniVariate analyses              ##############
###############################################################################      

  ### 4.1 Methylation (outcome) univariate
    ## a) Descriptive stats methylation
      # calculate the mean, median and standard deviation of % methylation
      # NOTE: uses luma_data_group data frame
        univar_meth <- luma_data_group %>%
          summarize (n = sum(!is.na(methylation)),
                     avg = round (mean (methylation, na.rm = T), 2),
                     median =  round (quantile (methylation, c(.5), na.rm = T),
                                      2),
                     stdev = round (sd (methylation, na.rm = T), 2),
                     avg_adjust = round (mean (meth_adjust, 
                                               na.rm = T), 2),
                     median_adjust =  round (quantile (meth_adjust, 
                                                       c(.5), na.rm = T), 2),
                     stdev_adjust = round (sd (meth_adjust, 
                                        na.rm = T), 2))
     
    ## b) save the data frame of summary stats out as a pdf into output file
        pdf(paste0(here(),"/output/output_luma_ecolog/univar_meth.pdf"),
            height = 2, width = 8)
        grid.table(univar_meth)
        dev.off()
    
  
  ### 4.2 Descriptive stats sex
    ## a) Sex ratio summary 
        sex_ratio_summary <- luma_data_group %>%
          #group_by (id) %>%
          group_by (sex) %>%
          summarise(n=n_distinct(id)) %>%
          mutate(freq = n / sum(n))
       
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/output_luma_ecolog/sex_ratio_summary.pdf"), 
                 height = 3, width = 5)
      grid.table(sex_ratio_summary)
      dev.off() 
     
      
  ### 4.3 Descriptive stats age    
    ## a) Age summary  
      age_var_summary <- luma_data_group %>%
        group_by (age.cat) %>%
        summarise (n.age = n(),
                   avg.age = round (mean(age.months, na.rm = T), 2),
                   med.age = round (median(age.months, na.rm = T), 2),
                   stdev.age = round (sd(age.months, na.rm = T), 2)) %>%
        mutate (freq.age =  (n.age / sum(n.age)))
          
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/output_luma_ecolog/age_var_summary.pdf"), 
          height = 4, width = 5)
      grid.table(age_var_summary)
      dev.off() 
    
      
  ### 4.4 Descriptive stats mom rank    
    ## a) Mom rank summary  
      mom_rank_summary <- luma_data_group %>%
        group_by (mom.strank.quart) %>%
        summarise (n.mom.rank = n(),
                   avg.mom.rank = round (mean(mom.strank, na.rm = T),2),
                   stdev.mom.rank = round (sd(mom.strank, na.rm = T), 2)) %>%
        mutate (freq.mom.rank =  (n.mom.rank / sum(n.mom.rank)))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/output_luma_ecolog/mom_rank_summary.pdf"), 
          height = 4, width = 8)
      grid.table(mom_rank_summary)
      dev.off() 
      
      
  ### 4.5 Descriptive stats intra_lit_rank
    ## a) Intra litter rank summary 
      intra_lit_rank_summary <- luma_data_group %>%
        #group_by (id) %>%
        group_by (lit.size) %>%
        summarise(n=n_distinct(id)) %>%
        mutate(freq = n / sum(n))
      
      ## b) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/output_luma_ecolog/",
                 "intra_lit_rank_summary.pdf"), 
          height = 3, width = 5)
      grid.table(intra_lit_rank_summary)
      dev.off() 
   
        
  ### 4.6 Descriptive stats human pop size
      ## a) Human pop size (proxy) summary 
      hum_pop_summary <- luma_data_group %>%
        #group_by (id) %>%
        group_by (hum.distrb) %>%
        summarise(n=n_distinct(id)) %>%
        mutate(freq = n / sum(n))
      
      ## b) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/output_luma_ecolog/",
                 "hum_pop_summary.pdf"), 
          height = 3, width = 5)
      grid.table(hum_pop_summary)
      dev.off() 
      
      
  ### 4.7 Descriptive stats prey density            
    ## a) Prey density summary 
      prey_peri_summary <- luma_data_group %>%
        summarise (dev.period = print("peri.concpt"),
                   n = sum(!is.na(total.peri.concpt)),
                   avg = round (mean(total.peri.concpt, 
                                               na.rm = T),2),
                   stdev = round (sd(total.peri.concpt, 
                                               na.rm = T), 2))
      
      prey_gest_summary <- luma_data_group %>%
        summarise (dev.period = print("gest"),
                   n = sum(!is.na(total.gest)),
                   avg = round (mean(total.gest, 
                                  na.rm = T),2),
                   stdev = round (sd(total.gest, 
                                  na.rm = T), 2))
      
      prey_birth.3_summary <- luma_data_group %>%
        summarise (dev.period = print("birth.3"),
                   n = sum(!is.na(total.birth.3)),
                   avg = round (mean(total.birth.3,
                                     na.rm = T),2),
                   stdev = round (sd(total.birth.3,
                                     na.rm = T), 2))
                   
                   
      prey_3.6_summary <- luma_data_group %>%
        summarise (dev.period = print("3.6"),
                   n = sum(!is.na(total.3.6)),               
                   avg = round (mean(total.3.6,
                                     na.rm = T),2),
                   stdev = round (sd(total.3.6,
                                     na.rm = T), 2))
                   
      prey_6.9_summary <- luma_data_group %>%
        summarise (dev.period = print("6.9"),
                   n = sum(!is.na(total.6.9)),              
                   avg = round (mean(total.6.9,
                                     na.rm = T),2),
                   stdev = round (sd(total.6.9,
                                     na.rm = T), 2))
      
    ## b) combine the prey density descriptive stats into a single data frame
      prey_var_summary <- rbind(prey_peri_summary, prey_gest_summary,
                                prey_birth.3_summary, prey_3.6_summary,
                                prey_6.9_summary)


    ## c) save the data frame of summary stats out as a pdf into output file
      pdf("output/output_luma_ecolog/prey_var_summary.pdf", height = 3, 
          width = 5)
      grid.table(prey_var_summary)
      dev.off()     
      
     

###############################################################################
##############             5.  Data transformations              ##############
###############################################################################
     
      
    ### 5.1 Center and Transform Predictive Variables
      ## a) center function based on column means    
        center_fxn <- function(x) {
          xcenter = colMeans(x, na.rm = T)
          x - rep(xcenter, rep.int(nrow(x), ncol(x)))
        }

      ## b) make a list of variable names to center, value = T is necessary
        # or column positions will be returned
        vars_to_center <- grep("total", names(luma_data_group), value = T)
        # drop total solids variable
        vars_to_center <- vars_to_center[vars_to_center != "total.solids"]
     
      ## c) Z-score standardize
        # Use 'scale' function to z-score standardization function based on 
        # subtracting mean from each x and dividing 
        # by 1 sd in the select columns; raw data are replaced with z-scores 
        # Here Prey densities have been z-score transformed
        # Z-score standardize luma_data
        # standardize luma_data
        luma_data[ ,c(vars_to_center)] <- 
          scale(luma_data[ ,c(vars_to_center)])
        # standardize luma_data_group
        luma_data_group[ ,c(vars_to_center)] <- 
          scale(luma_data_group[ ,c(vars_to_center)]) 
      
      
      
###############################################################################
##############              6.  Bi-variate analyses              ##############
###############################################################################  

  ### 6.1 Bivariate statistics methylation by sex
    ## a) Summary stats methylation by sex 
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg
      meth_by_sex <- luma_data_group %>%
        group_by (sex) %>%
        summarise (n.id = n(),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/output_luma_ecolog/meth_by_sex.pdf"),
          height = 6, width = 7)
      grid.table(meth_by_sex)
      dev.off()
    
    ## c) Plot Methylation by Sex
      # NOTE: uses luma_data_group
      ggplot(data = subset(luma_data_group, !is.na(x = sex)), 
             aes(x = sex, y = methylation, color = sex)) + 
        geom_boxplot() +
        theme(text = element_text(size=20))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        labs(title = "Percent Global DNA 
Mehtylation by Sex") +
        ylab("% Global DNA Methylation") +
        xlab("Sex")
      
    ## d) Save Plot
      # use ggsave to save the plot
      ggsave("meth_by_sex_plot.pdf", plot = last_plot(), device = NULL, 
             path = paste0(here(),"/output/output_luma_ecolog"), 
                           scale = 1, width = 7, height = 5,
                           units = c("in"), dpi = 300, limitsize = TRUE)  
      
    ## e) Bivariate Regression Methylatin by Sex
      # NOTE: here mulitple modeling options are coded and commented out;
      # can be used to compare different types of model fitting
      
      # Ordinary least squares method
      #sex.lm <- lm (methylation ~ sex, subset(luma_data_group, 
      #                                       !is.na(x = sex))) 
  
      # Maximum likelihood method
      #sex.glm <- glm(methylation ~ sex , 
      #               subset(luma_data_group,
      #                      !is.na(x = sex)),family = gaussian) 
     
      # uses 'nmle' package, which will provided p-value estimates
      sex.lme <- lme(methylation ~ sex, random =~1|id, 
                     subset(luma_data_group,!is.na(x = sex)))
      
      # Satterthwaite approximation of DF using lmerTest 
      #sex.sat <- lmerTest::lmer(methylation ~ sex + (1|id),  
      #                data = subset(luma_data_group, !is.na(x = sex)))
       
      #sex.lmer <- lme4::lmer(methylation ~ sex + (1|id),
      #                 data = subset(luma_data_group, !is.na(x = sex)))
   
      #summary(sex.lm)     # print model summary, effects and SE
      #confint(sex.lm)     # print 95% CIs for parameter estimates 
      #summary(sex.glm) 
      #confint(sex.glm)
      summary(sex.lme) 
      intervals(sex.lme, which = "fixed")
      #summary(sex.sat)  
      #confint(sex.sat)  
      #summary(sex.lmer)  
      #confint(sex.lmer)
      
      
  ### 6.2 Bivariate Statistics Methylation by Age 
    ## a) Summary stats methylation by age 
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg
      meth_by_age <- luma_data_group %>%
        group_by (age.cat) %>%
        summarise (n.id = n(),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/output_luma_ecolog/meth_by_age.pdf"), 
          height = 6, width = 7)
      grid.table(meth_by_age)
      dev.off()
    
    ## c) Plot Methylation by Age
      # graph of the raw data for percent global DNA methylaiton by maternal 
      # rank stratified by age 
      ggplot(data = subset(luma_data_group, !is.na(x = age.cat)), 
             aes(x = age.cat, y = methylation, color = age.cat)) + 
        geom_boxplot() +
        theme(text = element_text(size=20))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        labs(title = "Percent Global DNA 
Mehtylation by Age") +
        ylab("% Global DNA Methylation") +
        xlab("Categorical Age")
      
    ## d) Save Plot
      # use ggsave to save the plot
      ggsave("meth_by_age_plot.pdf", plot = last_plot(), device = NULL, 
             path = paste0(here(),"/output/output_luma_ecolog"),
             scale = 1, width = 7, height = 5, 
             units = c("in"), dpi = 300, limitsize = TRUE)  
  
    ## e) Bivariate Regression Methylatino by Age
      # uses 'nmle' package, which will provided p-value estimates
      age.lme <- lme(methylation ~ age.cat, random =~1|id, 
                     subset(luma_data_group,!is.na(x = age.cat)))
      
      summary(age.lme)            #  print model summary, effects and SE
      intervals(age.lme, 
                which = "fixed")  # print 95% CIs for parameter estimates
      anova.lme(age.lme)          # generate p-value from Wald test
     
      # Satterthwaite approximation of DF using lmerTest 
      #age.lmer <- lmer(methylation ~ Age + (1|ID), 
      #              data = subset(luma_data_group, !is.na(x = Age)))
   
      #summary(age.sat)  # print model summary, effects and SE
      #confint(age.sat)  # print 95% CIs for parameter estimates
      #anova(age.sat) # type 3 ANOVA with Satterthwaite DF
      #means.age.sat <- Effect("Age", age.mod) # use Effects to see group means


  ### 6.3 Bivariate Statistics Methylation by Rank
    ## a) Summary stats methylation by age 
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg
      meth_by_mom_rank <- luma_data_group %>%
        group_by (mom.strank.quart) %>%
        summarise (n.id = n(),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))

    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/output_luma_ecolog/meth_by_mom_rank.pdf"), 
          height = 6, width = 7)
      grid.table(meth_by_mom_rank)
      dev.off()
    
    ## c) Plot mehtylation by rank
      # graph of the raw data for percent global DNA methylaiton by maternal 
      # rank
      ggplot(data = subset(luma_data_group, !is.na(x = mom.strank.quart)),
                           aes(x = mom.strank.quart, y = methylation,
                                    color = mom.strank.quart)) + 
        geom_boxplot() +
        theme(text = element_text(size=20))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        labs(title = "Percent Global DNA 
Mehtylation by Maternal Rank") +
        ylab("% Global DNA Methylation") +
        xlab("Maternal Rank Quartiles")
    
    ## d) Save Plot
      # use ggsave to save the plot
      ggsave("meth_by_mom_rank_plot.pdf", plot = last_plot(), device = NULL, 
             path = paste0(here(),"/output/output_luma_ecolog"),
             scale = 1, width = 7, height = 5, 
             units = c("in"), dpi = 300, limitsize = TRUE)  
      
    ## e) Bivariate Regression Methylatino by maternal rank 
      # uses 'nmle' package, which will provided p-value estimates
      mom.rank.lme <- lme(methylation ~ mom.strank.quart, random =~1|id, 
                     subset(luma_data_group,!is.na(x = mom.strank.quart)))
      
      summary(mom.rank.lme)       #  print model summary, effects and SE
      intervals(mom.rank.lme, 
                which = "fixed")  # print 95% CIs for parameter estimates
      anova(mom.rank.lme)         # generate p-value from Wald test
  

  ### 6.4 Bivariate Statistics Methylation by intra litter rank
    ## a) Summary stats methylation by lit.size 
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg
      meth_by_intra_lit_rank <- luma_data_group %>%
        group_by (lit.size) %>%
        summarise (n.id = n(),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/output_luma_ecolog/",
                 "meth_by_intra_lit_rank.pdf"), 
          height = 6, width = 7)
      grid.table( meth_by_intra_lit_rank)
      dev.off()
      
    ## c) Plot mehtylation by intra litter rank
      # graph of the raw data for percent global DNA methylaiton by intra
      # litter rank
      ggplot(data = subset(luma_data_group, !is.na(x = lit.size)),
             aes(x = lit.size, y = methylation,
                 color = lit.size)) + 
        geom_boxplot() +
        theme(text = element_text(size=20))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        labs(title = "Percent Global DNA 
Mehtylation by Intra Litter Rank") +
        ylab("% Global DNA Methylation") +
        xlab("Intra Litter Rank")
      
      ## d) Save Plot
      # use ggsave to save the plot
      ggsave("meth_by_intra_lit_rank_plot.pdf", plot = last_plot(), 
             device = NULL, 
             path = paste0(here(),"/output/output_luma_ecolog"),
             scale = 1, width = 7, height = 5, 
             units = c("in"), dpi = 300, limitsize = TRUE)  
      
      ## e) Bivariate Regression Methylatino by maternal rank 
      # uses 'nmle' package, which will provided p-value estimates
      lit.size.lme <- lme(methylation ~ lit.size, random =~1|id, 
                          subset(luma_data_group,!is.na(x = lit.size)))
      
      summary(lit.size.lme)   #  print model summary, effects and SE
      intervals(lit.size.lme, 
                which = "fixed")    # print 95% CIs for parameter estimates
      anova(lit.size.lme)     # generate p-value from Wald test
      
      
  ### 6.5 Bivariate Statistics Methylation by human population size
    ## a) Summary stats methylation by hum.distrb
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg
      meth_by_hum_pop <- luma_data_group %>%
        group_by (hum.distrb) %>%
        summarise (n.id = n(),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/output_luma_ecolog/",
                 "meth_by_hum_pop.pdf"), 
          height = 6, width = 7)
      grid.table(meth_by_hum_pop)
      dev.off()
      
    ## c) Plot mehtylation by human population size
      # graph of the raw data for percent global DNA methylaiton by human 
      # populaiton size
      ggplot(data = subset(luma_data_group, !is.na(x = hum.distrb)),
             aes(x = hum.distrb, y = methylation,
                 color = hum.distrb)) + 
        geom_boxplot() +
        theme(text = element_text(size=20))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        labs(title = "Percent Global DNA 
Mehtylation by Human Population Size") +
        ylab("% Global DNA Methylation") +
        xlab("Intra Litter Rank")
      
    ## d) Save Plot
      # use ggsave to save the plot
      ggsave("meth_by_hum_pop_plot.pdf", plot = last_plot(), 
             device = NULL, 
             path = paste0(here(),"/output/output_luma_ecolog"),
             scale = 1, width = 7, height = 5, 
             units = c("in"), dpi = 300, limitsize = TRUE)  
      
    ## e) Bivariate Regression Methylatino by human population size
      # uses 'nmle' package, which will provided p-value estimates
      hum.distrb.lme <- lme(methylation ~ hum.distrb, random =~1|id, 
                                subset(luma_data_group,!is.na(x = hum.distrb)))
      
      summary(hum.distrb.lme)        #  print model summary, effects and SE
      intervals(hum.distrb.lme, 
                which = "fixed")  # print 95% CIs for parameter estimates
      anova(hum.distrb.lme)          # generate p-value from Wald test
      
    
  ### 6.6 Bivariate statistics methylation by periconceptional prey density 
    ## a) Summary stats methylation by total.peri.concpt
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg
      #meth_by_peri_concpt <- luma_data_group %>%
      #  summarise (n.id = n(),
      #             avg = round (mean(methylation, na.rm = T), 2),
      #             median =  round (quantile(methylation, c(.5), na.rm = T), 2),
      #             sd = round (sd(methylation, na.rm = T), 2))
      
      ## b) save the data frame of summary stats out as a pdf into output file
      #pdf(paste0(here(),"/output/output_luma_ecolog/",
      #           "meth_by_peri_concpt.pdf"), 
      #    height = 6, width = 7)
      #grid.table(meth_by_peri_concpt)
      #dev.off() 
      
    ## b) Bivariate regression methylatino by by pericoceptional total prey
      # density
      # uses 'nmle' package, which will provided p-value estimates
      peri.concpt.lme <- lme(methylation ~ total.peri.concpt, random =~1|id, 
                         subset(luma_data_group,!is.na(x = total.peri.concpt)))
      
      summary(peri.concpt.lme)        #  print model summary, effects and SE
      intervals(peri.concpt.lme, 
                which = "fixed")  # print 95% CIs for parameter estimates

       
  ### 6.7 Bivariate statistics methylation by gestational prey density
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg
      
    ## a) Bivariate regression methylatino by gestation total prey
      # density
      # uses 'nmle' package, which will provided p-value estimates
      total.gest.lme <- lme(methylation ~ total.gest, random =~1|id, 
                             subset(luma_data_group,!is.na(x = total.gest)))
      
      summary(total.gest.lme)        #  print model summary, effects and SE
      intervals(total.gest.lme, 
                which = "fixed")  # print 95% CIs for parameter estimates
      
      
  ### 6.8 Bivariate statistics methylation by birth to 3 months prey density
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg
      
    ## a) Bivariate regression methylatino by birth to 3 months total prey
      # density
      # uses 'nmle' package, which will provided p-value estimates
      total.birth.3.lme <- lme(methylation ~ total.birth.3, random =~1|id, 
                            subset(luma_data_group,
                                   !is.na(x = total.birth.3)))
      
      summary(total.birth.3.lme)        #  print model summary, effects and SE
      intervals(total.birth.3.lme, 
                which = "fixed")  # print 95% CIs for parameter estimates
      
  
  ### 6.9 Bivariate statistics methylation by 3 to 6 months prey density
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg
      
    ## a) Bivariate regression methylatino by 3 to 6 months total prey
      # density
      # uses 'nmle' package, which will provided p-value estimates
      total.3.6.lme <- lme(methylation ~ total.3.6, random =~1|id, 
                               subset(luma_data_group,
                                      !is.na(x = total.3.6)))
      
      summary(total.3.6.lme)        #  print model summary, effects and SE
      intervals(total.3.6.lme, 
                which = "fixed")  # print 95% CIs for parameter estimates
      
      
  ### 6.10 Bivariate statistics methylation by 6 to 9 months prey density
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg
      
    ## a) Bivariate regression methylatino by 6 to 9 months total prey
      # density
      # uses 'nmle' package, which will provided p-value estimates
      total.6.9.lme <- lme(methylation ~ total.6.9, random =~1|id, 
                               subset(luma_data_group,!is.na(x = total.6.9)))
      
      summary(total.6.9.lme)        #  print model summary, effects and SE
      intervals(total.6.9.lme, 
                which = "fixed")  # print 95% CIs for parameter estimates

    
  ### 6.11 View methylation by rank within age strata
    ## a) Graph of the raw data for percent global DNA methylaiton by maternal rank
      # stratified by age
      ggplot(subset(luma_data_group, !is.na(x = age.cat)),
             aes(x = mom.strank.quart.order, y = methylation, 
                 color = age.cat )) +
        geom_point(shape = 1) +
        theme(text = element_text(size=20))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        geom_smooth(method = lm, se = F) + # Add linear regression best fit lines
        labs(title = "Percent Global DNA Mehtylation by 
             Maternal Rank by Age",
             fill = "age") +
        theme(plot.title = element_text(hjust = 0.5))+
        ylab("% Global DNA Methylation") +
        xlab("Maternal Rank")
      
    ## b) Save Plot
      # use ggsave to save the linearization plot
      ggsave("meth_by_rank_by_age_lm_plot.pdf", plot = last_plot(), 
             device = NULL,
             path = paste0(here(),"/output/output_luma_ecolog"), 
             scale = 1, width = 7, height = 5,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## c) Box-plots of methylation by rank within age categories  
      ggplot(subset(luma_data_group, 
                    !is.na(x = mom.strank.quart)&!is.na(x=age.cat)),
             aes(y = methylation, x = factor(mom.strank.quart), 
                 fill = age.cat))+
        geom_boxplot() +
        theme(text = element_text(size=20))+
        labs(title = "Percent Global DNA Mehtylation by 
             Maternal Rank by Age",
             fill = "age") +
        theme(plot.title = element_text(hjust = 0.5))+
        ylab("% Global DNA Methylation") +
        xlab("Maternal Rank")
      
    ## d) Save Plot
      # use ggsave to save the linearization plot
      ggsave("meth_by_rank_by_age_boxplot.pdf", plot = last_plot(), 
             device = NULL,
             path = paste0(here(),"/output/output_luma_ecolog"),
             scale = 1, width = 7,height = 5,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
      
  ### 6.12 Stratify data by life history variables
    ## a) Cub subset that includes both females and males
      luma_data_cub <- luma_data_group %>%
        filter(grepl('^cub$', age.cat))
      
      #**hack**# Fill in two missing maternal rank quartiles
      # change mom.absrank to an integer
      luma_data_cub$mom.absrank <- as.integer(as.character
                                              (luma_data_cub$mom.absrank))
      # fill in missing mom.strank.quart info based on mom.absrank
      luma_data_cub$mom.strank.quart <- ifelse(
        (is.na(luma_data_cub$mom.strank.quart) & 
           luma_data_cub$mom.absrank > 18),
        "Q1 (lowest)", as.character(luma_data_cub$mom.strank.quart))
      
      # re-order to sets the reference level for high
      luma_data_cub <- transform (luma_data_cub, 
                                  mom.strank.quart = factor(mom.strank.quart,
                                              levels = c("Q1 (lowest)",
                                                         "Q2","Q3",
                                                         "Q4 (highest)")))
      
    ## b) Subadult subset that includes both females and males
      luma_data_sub <- luma_data_group %>%
        dplyr::filter(grepl('sub', age.cat))
      
    ## c) Adult subset that includes both females and males
      luma_data_adult <- luma_data_group %>%
        dplyr::filter(grepl('^adult$', age.cat))  
      
      
      
###############################################################################
##############                  7. Cub models                    ##############
###############################################################################      

        
  ### 7.1 Cub model: methylation by mom rank
    ## a) Check within strata descritpive stats
      luma_data_cub %>%
        group_by (mom.strank.quart) %>%
        summarise (n.id = n(),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Mehtylation by maternal rank quartile
      cub.mom.rank.mod <- glm(methylation ~ mom.strank.quart + sex + age.mon,
                         data = luma_data_cub)
      
    ## c) Parameter estimates
      summary(cub.mom.rank.mod)  # print model summary, effects and SE
      confint(cub.mom.rank.mod)  # print 95% CIs for parameter estimates
    
    ## d) Extract mom.strank.quart estimates and 
      cub.rank.ef <- effect("mom.strank.quart", cub.mom.rank.mod)
      summary(cub.rank.ef)
      # Save effects as data frame
      cub.rank.ef.table <- as.data.frame(cub.rank.ef)
      # Set the reference level to Q1
      cub.rank.ef.table <- transform( cub.rank.ef.table,
                          mom.strank.quart = factor(mom.strank.quart,
                                                    levels = c("Q1 (lowest)", 
                                                              "Q2", "Q3", 
                                                              "Q4 (highest)")))
      
    ## e) Graph cub.mom.rank effects
      ggplot(cub.rank.ef.table, aes(x = mom.strank.quart, y = fit)) +
        geom_point() +
        geom_errorbar(aes(ymin= fit-se, ymax= fit+se), width=0.4) +
        theme(text = element_text(size=20))+
        #theme_bw(base_size=12) +
        labs(title = "Cub Percent Global DNA Mehtylation 
  by Maternal Rank") +
        theme(plot.title = element_text(hjust = 0.5))+
        ylab("% Global DNA Methylation ± SE") +
        xlab("Maternal Rank")
    
    ## g) Save Plot
      # use ggsave to save the linearization plot
      ggsave("mat_rank_cub_mod_beta.pdf", plot = last_plot(), device = NULL,
             path = "./output/output_luma_ecolog", scale = 1, width = 7,
             height = 5.5,
             units = c("in"), dpi = 300, limitsize = TRUE)
 
    ## h) Do a post-hoc test to determine if Q2, Q3, and Q4 differ
      pairwise.t.test(luma_data_cub$methylation,
                      luma_data_cub$mom.strank.quart,
                      p.adj = "none")
      TukeyHSD(aov(methylation ~ mom.strank.quart + sex,
                   data = luma_data_cub))
    
    ## i) Combine quartiles 2-4 into a single category
      luma_data_cub$mom.strank.quart.comb <- as.factor(
      ifelse(luma_data_cub$mom.strank.quart.order <= 1,
             "Q1 (lowest)", "Q2-Q4 (highest)"))
    
    
    ## j)  Mehtylation by Maternal Rank Quartiles binned
      cub.rank.mod2 <- glm(methylation ~ mom.strank.quart.comb + sex + 
                             age.mon,
                         data = luma_data_cub)
      
      summary(cub.rank.mod2)  # print model summary, effects and SE
      confint(cub.rank.mod2)  # print 95% CIs for parameter estimates
      
    
  ### 7.2 Cub model: methylation by intra litter rank      
    ## a) Check within strata descritpive stats
    luma_data_cub %>%
        group_by (lit.size) %>%
        summarise (n.id = n(),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Mehtylation by intra litter rank
      cub.intra.rank.mod <- glm(methylation ~ lit.size + sex + 
                                  age.mon, data = luma_data_cub)
    ## c) Parameter estimates
    summary(cub.intra.rank.mod)  # print model summary, effects and SE
    confint(cub.intra.rank.mod)  # print 95% CIs for parameter estimates  
    
    ## d) Combine quartiles tiwn hi and lo into a single category
    luma_data_cub$lit.size.comb <- as.factor(
      ifelse(luma_data_cub$lit.size == "single",
             "single", "twin"))
    
    ## e)  Mehtylation by litter size
    cub.intra.rank.mod2 <- glm(methylation ~ lit.size.comb  + sex + 
                           age.mon,
                         data = luma_data_cub)
    
    summary(cub.intra.rank.mod2)  # print model summary, effects and SE
    confint(cub.intra.rank.mod2)  # print 95% CIs for parameter estimates
    
  ### 7.3 Cub model: methylation by human population size     
    ## a) Check within strata descritpive stats
    luma_data_cub %>%
      group_by (hum.distrb) %>%
      summarise (n.id = n(),
                 avg = round (mean(methylation, na.rm = T), 2),
                 median =  round (quantile(methylation, c(.5), na.rm = T), 2),
                 sd = round (sd(methylation, na.rm = T), 2))
    
    ## b) Mehtylation by intra litter rank
    cub.hum.distrb.mod <- glm(methylation ~ hum.distrb + sex + 
                                age.mon, data = luma_data_cub)
    ## b) Parameter estimates
    summary(cub.hum.distrb.mod)  # print model summary, effects and SE
    confint(cub.hum.distrb.mod)  # print 95% CIs for parameter estimates  
  
    
  ### 7.4 Cub model: methylation by prey density     
    ## a) Mehtylation by periconceptional prey density
      cub.peri.prey.mod <- glm(methylation ~ total.peri.concpt + sex + 
                                age.mon, data = luma_data_cub)
    ## b) Parameter estimates
      summary(cub.peri.prey.mod)  # print model summary, effects and SE
      confint(cub.peri.prey.mod)  # print 95% CIs for parameter estimates
      
    ## c) Mehtylation by gestational prey density
      cub.gest.prey.mod <- glm(methylation ~ total.gest + total.peri.concpt +
                                 sex +  age.mon, data = luma_data_cub)
    ## d) Parameter estimates
      summary(cub.gest.prey.mod)  # print model summary, effects and SE
      confint(cub.gest.prey.mod)  # print 95% CIs for parameter estimates  
      
    ## e) Mehtylation by birth to 3 months prey density
      cub.birth.3.prey.mod <- glm(methylation ~ total.birth.3 + total.gest + 
                                    total.peri.concpt + sex + age.mon, 
                                  data = luma_data_cub)
    ## f) Parameter estimates
      summary(cub.birth.3.prey.mod)  # print model summary, effects and SE
      confint(cub.birth.3.prey.mod)  # print 95% CIs for parameter estimates
      
    ## g) Mehtylation by 3 to 6 months prey density
      cub.3.6.prey.mod <- glm(methylation ~ total.3.6 + total.birth.3 + 
                                total.gest + total.peri.concpt +
                                  + sex + age.mon, data = luma_data_cub)
    ## h) Parameter estimates
      summary(cub.3.6.prey.mod)  # print model summary, effects and SE
      confint(cub.3.6.prey.mod)  # print 95% CIs for parameter estimates
      
    ## i) Mehtylation by 3 to 6 months prey density
      cub.6.9.prey.mod <- glm(methylation ~ total.6.9 + total.3.6 +
                              total.birth.3 + total.gest + total.peri.concpt 
                              + sex + age.mon, data = luma_data_cub)
    
    ## j) Parameter estimates
      summary(cub.6.9.prey.mod)  # print model summary, effects and SE
      confint(cub.6.9.prey.mod)  # print 95% CIs for parameter estimates
  
        
#******************************************************************************#  
  ### 7.5 Cub interaction model: mom rank by prey density   
    ## a) Mehtylation by periconceptional prey density
      cub.peri.x.mom.rank.mod <- glm(methylation ~ 
                                       mom.strank.quart.comb * 
                                       total.peri.concpt + 
                                       sex + age.mon, data = luma_data_cub)
      # Parameter estimates
      summary(cub.peri.x.mom.rank.mod)  # print model summary, effects and SE
      confint(cub.peri.x.mom.rank.mod)  # print 95% CIs for parameter estimates
    
    ## b) Mehtylation by gestation prey density
      cub.gest.x.mom.rank.mod <- glm(methylation ~ 
                                       mom.strank.quart.comb * 
                                       total.gest + total.peri.concpt +
                                       sex + age.mon, data = luma_data_cub)
      # Parameter estimates
      summary(cub.gest.x.mom.rank.mod)  # print model summary, effects and SE
      confint(cub.gest.x.mom.rank.mod)  # print 95% CIs for parameter estimates
      
    ## c) Mehtylation by birth to 3 month prey density
      cub.birth.3.x.mom.rank.mod <- glm(methylation ~ 
                                       mom.strank.quart.comb * 
                                       total.birth.3 + total.peri.concpt + 
                                       total.gest +
                                       sex + age.mon, data = luma_data_cub)
      # Parameter estimates
      summary(cub.birth.3.x.mom.rank.mod) # print model summary, effects and SE
      confint(cub.birth.3.x.mom.rank.mod)# print 95% CIs for parameter estimates
    
    ## d) Stratified analysis
      # subset data based hi vs lo than median birth.3 prey density
      luma_data_cub_birth.3_hi<- luma_data_cub %>%
        filter(total.birth.3 > median(total.birth.3, na.rm = T))
      
      luma_data_cub_birth.3_lo<- luma_data_cub %>%
        filter(total.birth.3 <= median(total.birth.3, na.rm = T))
    
      # re-run models on stratified data
      # birth.3 prey hi
      cub.birth.3.hi.mom.rank.mod <- glm(methylation ~ 
                                          mom.strank.quart.comb +
                                          sex + age.mon, 
                                          data = luma_data_cub_birth.3_hi)
      # Parameter estimates
      summary(cub.birth.3.hi.mom.rank.mod) 
      confint( cub.birth.3.hi.mom.rank.mod)
      
      # re-run models on stratified data
      # birth.3 prey lo
      cub.birth.3.lo.mom.rank.mod <- glm(methylation ~ 
                                           mom.strank.quart.comb +
                                           sex + age.mon, 
                                         data = luma_data_cub_birth.3_lo)
      # Parameter estimates
      summary(cub.birth.3.lo.mom.rank.mod) 
      confint( cub.birth.3.lo.mom.rank.mod)
   
      
    ## e) Mehtylation by 3 to 6 month prey density
      cub.3.6.x.mom.rank.mod <- glm(methylation ~ 
                                          mom.strank.quart.comb * 
                                          total.3.6 + total.peri.concpt + 
                                          total.gest + total.birth.3 +
                                          sex + age.mon, data = luma_data_cub)
      # Parameter estimates
      summary(cub.3.6.x.mom.rank.mod) # print model summary, effects and SE
      confint(cub.3.6.x.mom.rank.mod) # print 95% CIs for parameter estimates
      
    ## f) Mehtylation by 6 to 9 month prey density
      cub.6.9.x.mom.rank.mod <- glm(methylation ~ 
                                      mom.strank.quart.comb * 
                                      total.6.9 + total.peri.concpt + 
                                      total.gest + total.birth.3 + total.3.6 +
                                      sex + age.mon, data = luma_data_cub)
      # Parameter estimates
      summary(cub.6.9.x.mom.rank.mod) # print model summary, effects and SE
      confint(cub.6.9.x.mom.rank.mod) # print 95% CIs for parameter estimates
  
        
  ### 7.6 Cub interaction model: mom rank by human pop   
    ## a) Maternal rank by human population size
      hum.distrb.x.mom.rank.mod <- glm(methylation ~ 
                                       mom.strank.quart.comb * hum.distrb + 
                                       sex + age.mon, data = luma_data_cub)
    ## b) Parameter estimates
      summary(hum.distrb.x.mom.rank.mod)  # print model summary, effects and SE
      confint(hum.distrb.x.mom.rank.mod)  # print 95% CIs for parameter estimates
    
    ## c) Cub subset that human pop hi
      luma_data_cub_hum_hi<- luma_data_cub %>%
        dplyr::filter(grepl('hi', hum.distrb))
    
      
    ## d)  Mehtylation by Maternal Rank Quartiles binned hi human pop
      cub.rank.hum.hi<- glm(methylation ~ mom.strank.quart.comb + sex + 
                             age.mon,
                           data = luma_data_cub_hum_hi)
      
      summary(cub.rank.hum.hi)  # print model summary, effects and SE
      confint(cub.rank.hum.hi)  # print 95% CIs for parameter estimates
    
    ## e) Cub subset that human pop lo
      luma_data_cub_hum_lo<- luma_data_cub %>%
        filter(grepl('lo', hum.distrb))
      
    ## f)  Mehtylation by Maternal Rank Quartiles binned lo human pop
      cub.rank.hum.lo<- glm(methylation ~ mom.strank.quart.comb + sex + 
                              age.mon,
                            data = luma_data_cub_hum_lo)
      
      summary(cub.rank.hum.lo)  # print model summary, effects and SE
      confint(cub.rank.hum.lo)  # print 95% CIs for parameter estimates


  
###############################################################################
##############                 8. Subadult models                ##############
###############################################################################     
    
      
  ### 8.1 Subadult model: methylation by mom rank
    ## a) Check within strata descritpive stats
      luma_data_sub %>%
        group_by (mom.strank.quart) %>%
        summarise (n.id = n(),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Mehtylation by maternal rank quartile
      sub.mom.rank.mod <- glm(methylation ~ mom.strank.quart + sex + age.mon,
                              data = luma_data_sub)
      
    ## c) Parameter estimates
      summary(sub.mom.rank.mod)  # print model summary, effects and SE
      confint(sub.mom.rank.mod)  # print 95% CIs for parameter estimates
      
    ## d) Extract mom.strank.quart estimates and 
      sub.rank.ef <- effect("mom.strank.quart", sub.mom.rank.mod)
      summary(sub.rank.ef)
      # Save effects as data frame
      sub.rank.ef.table <- as.data.frame(sub.rank.ef)
      # Set the reference level to Q1
      sub.rank.ef.table <- transform( sub.rank.ef.table,
                         mom.strank.quart = factor(mom.strank.quart,
                                                   levels = c("Q1 (lowest)",
                                                              "Q2", "Q3",
                                                              "Q4 (highest)")))
      
    ## e) Graph sub.mom.rank effects
      ggplot(sub.rank.ef.table, aes(x = mom.strank.quart, y = fit)) +
        geom_point() +
        geom_errorbar(aes(ymin= fit-se, ymax= fit+se), width=0.4) +
        theme(text = element_text(size=20))+
        #theme_bw(base_size=12) +
        labs(title = "Subadult Percent Global DNA Mehtylation 
             by Maternal Rank") +
        theme(plot.title = element_text(hjust = 0.5))+
        ylab("% Global DNA Methylation ± SE") +
        xlab("Maternal Rank")
      
    ## g) Save Plot
      # use ggsave to save the linearization plot
      ggsave("mat_rank_sub_mod_beta.pdf", plot = last_plot(), device = NULL,
             path = "./output/output_luma_ecolog", scale = 1, width = 7,
             height = 5.5,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## h) Do a post-hoc test to determine if Q2, Q3, and Q4 differ
      pairwise.t.test(luma_data_sub$methylation,
                      luma_data_sub$mom.strank.quart,
                      p.adj = "none")
      TukeyHSD(aov(methylation ~ mom.strank.quart + sex,
                   data = luma_data_sub))
      
    ## i) Combine quartiles 2-4 into a single category
      luma_data_sub$mom.strank.quart.comb <- as.factor(
        ifelse(luma_data_sub$mom.strank.quart.order <= 1,
               "Q1 (lowest)", "Q2-Q4 (highest)"))
      
    ## j)  Mehtylation by Maternal Rank Quartiles binned
      sub.rank.mod2 <- glm(methylation ~ mom.strank.quart.comb + sex + 
                             age.mon,
                           data = luma_data_sub)
      
      summary(sub.rank.mod2)  # print model summary, effects and SE
      confint(sub.rank.mod2)  # print 95% CIs for parameter estimates
      
      
  ### 8.2 Sub model: methylation by intra litter rank      
    ## a) Check within strata descritpive stats
      luma_data_sub %>%
        group_by (lit.size) %>%
        summarise (n.id = n(),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Mehtylation by intra litter rank
      sub.intra.rank.mod <- glm(methylation ~ lit.size + sex + 
                                  age.mon, data = luma_data_sub)
    ## c) Parameter estimates
      summary(sub.intra.rank.mod)  # print model summary, effects and SE
      confint(sub.intra.rank.mod)  # print 95% CIs for parameter estimates  
      
    ## d) Combine quartiles tiwn hi and lo into a single category
      luma_data_sub$lit.size.comb <- as.factor(
        ifelse(luma_data_sub$lit.size == "single",
               "single", "twin"))
      
    ## e)  Mehtylation by litter size
      sub.intra.rank.mod2 <- glm(methylation ~ lit.size.comb  + sex + 
                                   age.mon,
                                 data = luma_data_sub)
      # Parameter estimates
      summary(sub.intra.rank.mod2)  # print model summary, effects and SE
      confint(sub.intra.rank.mod2)  # print 95% CIs for parameter estimates
      
  
  ### 8.3 Sub model: methylation by human population size     
    ## a) Check within strata descritpive stats
      luma_data_sub %>%
        group_by (hum.distrb) %>%
        summarise (n.id = n(),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Mehtylation by intra litter rank
      sub.hum.distrb.mod <- glm(methylation ~ hum.distrb + sex + 
                               age.mon, data = luma_data_sub)
    ## c) Parameter estimates
      summary(sub.hum.distrb.mod)  # print model summary, effects and SE
      confint(sub.hum.distrb.mod)  # print 95% CIs for parameter estimates  
      
      
  ### 8.4 sub model: methylation by prey density     
    ## a) Mehtylation by periconceptional prey density
      sub.peri.prey.mod <- glm(methylation ~ total.peri.concpt + sex + 
                                 age.mon, data = luma_data_sub)
    ## b) Parameter estimates
      summary(sub.peri.prey.mod)  # print model summary, effects and SE
      confint(sub.peri.prey.mod)  # print 95% CIs for parameter estimates
      
    ## c) Mehtylation by gestational prey density
      sub.gest.prey.mod <- glm(methylation ~ total.gest + total.peri.concpt +
                                 sex +  age.mon, data = luma_data_sub)
    ## d) Parameter estimates
      summary(sub.gest.prey.mod)  # print model summary, effects and SE
      confint(sub.gest.prey.mod)  # print 95% CIs for parameter estimates  
      
    ## e) Mehtylation by birth to 3 months prey density
      sub.birth.3.prey.mod <- glm(methylation ~ total.birth.3 + total.gest +
                                  total.peri.concpt + sex + age.mon, 
                                  data = luma_data_sub)
    ## f) Parameter estimates
      summary(sub.birth.3.prey.mod)  # print model summary, effects and SE
      confint(sub.birth.3.prey.mod)  # print 95% CIs for parameter estimates
      
      
    ## g) Mehtylation by 3 to 6 months prey density
      sub.3.6.prey.mod <- glm(methylation ~ total.3.6 + total.birth.3 +
                              total.gest + total.peri.concpt + sex + age.mon, 
                              data = luma_data_sub)
    ## h) Parameter estimates
      summary(sub.3.6.prey.mod)  # print model summary, effects and SE
      confint(sub.3.6.prey.mod)  # print 95% CIs for parameter estimates
      
    ## i) Mehtylation by 3 to 6 months prey density
      sub.6.9.prey.mod <- glm(methylation ~ total.6.9 + total.3.6 +
                              total.birth.3 + total.gest + total.peri.concpt + 
                                sex + age.mon, data = luma_data_sub)
    ## j) Parameter estimates
      summary(sub.6.9.prey.mod)  # print model summary, effects and SE
      confint(sub.6.9.prey.mod)  # print 95% CIs for parameter estimates
      
    
#******************************************************************************#
  ### 8.5 Sub interaction  model: litter sizeprey density   
    ## a) Mehtylation by periconceptional prey density
      sub.peri.x.lit.size.mod <- glm(methylation ~ 
                                       lit.size.comb * 
                                       total.peri.concpt + 
                                       sex + age.mon, data = luma_data_sub)
      # Parameter estimates
      summary(sub.peri.x.lit.size.mod)  # print model summary, effects and SE
      confint(sub.peri.x.lit.size.mod)  # print 95% CIs for parameter estimates
      
    ## b) Mehtylation by gestation prey density
      sub.gest.x.lit.size.mod <- glm(methylation ~ 
                                       lit.size.comb * 
                                       total.gest + total.peri.concpt +
                                       sex + age.mon, data = luma_data_sub)
      # Parameter estimates
      summary(sub.gest.x.lit.size.mod)  # print model summary, effects and SE
      confint(sub.gest.x.lit.size.mod)  # print 95% CIs for parameter estimates
      
    ## c) Mehtylation by birth to 3 month prey density
      sub.birth.3.x.lit.size.mod <- glm(methylation ~ 
                                          lit.size.comb * 
                                          total.birth.3 + total.peri.concpt + 
                                          total.gest +
                                          sex + age.mon, data = luma_data_sub)
      # Parameter estimates
      summary(sub.birth.3.x.lit.size.mod) # print model summary, effects and SE
      confint(sub.birth.3.x.lit.size.mod)# print 95% CIs for parameter estimates
      
    ## d) Mehtylation by 3 to 6 month prey density
      sub.3.6.x.lit.size.mod <- glm(methylation ~ 
                                      lit.size.comb * 
                                      total.3.6 + total.peri.concpt + 
                                      total.gest + total.birth.3 +
                                      sex + age.mon, data = luma_data_sub)
      # Parameter estimates
      summary(sub.3.6.x.lit.size.mod) # print model summary, effects and SE
      confint(sub.3.6.x.lit.size.mod) # print 95% CIs for parameter estimates
      
    ## e) Mehtylation by 6 to 9 month prey density
      sub.6.9.x.lit.size.mod <- glm(methylation ~ 
                                      lit.size.comb * 
                                      total.6.9 + total.peri.concpt + 
                                      total.gest + total.birth.3 + total.3.6 +
                                      sex + age.mon, data = luma_data_sub)
      # Parameter estimates
      summary(sub.6.9.x.lit.size.mod) # print model summary, effects and SE
      confint(sub.6.9.x.lit.size.mod) # print 95% CIs for parameter estimates
      
      
  ### 8.6 sub model: litter size by human pop   
    ## a) Mehtylation by periconceptional prey density
      hum.distrb.x.lit.size.mod <- glm(methylation ~ 
                                      lit.size.comb * hum.distrb + 
                                      sex + age.mon, data = luma_data_sub)
      # Parameter estimates
      summary(hum.distrb.x.lit.size.mod)  # print model summary, effects and SE
      confint(hum.distrb.x.lit.size.mod)  # print 95% CIs for parameter estimates
      

      
###############################################################################
##############                  9. Adult models                  ##############
############################################################################### 
      
  ### 9.1 Adult model: methylation by mom rank
    ## a) Check within strata descritpive stats
      luma_data_adult %>%
        group_by (mom.strank.quart) %>%
        summarise (n.id = n(),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Mehtylation by maternal rank quartile
      adult.mom.rank.mod <- glm(methylation ~ mom.strank.quart + sex + 
                                  age.mon, data = luma_data_adult)
      
    ## c) Parameter estimates
      summary(adult.mom.rank.mod)  # print model summary, effects and SE
      confint(adult.mom.rank.mod)  # print 95% CIs for parameter estimates
      
    ## d) Extract mom.strank.quart estimates and 
      adult.rank.ef <- effect("mom.strank.quart", adult.mom.rank.mod)
      summary(adult.rank.ef)
      # Save effects as data frame
      adult.rank.ef.table <- as.data.frame(adult.rank.ef)
      # Set the reference level to Q1
      adult.rank.ef.table <- transform(adult.rank.ef.table,
                          mom.strank.quart = factor(mom.strank.quart,
                                                    levels = c("Q1 (lowest)",
                                                               "Q2", "Q3",
                                                               "Q4 (highest)")))
      
    ## e) Graph adult.mom.rank effects
      ggplot(adult.rank.ef.table, aes(x = mom.strank.quart, y = fit)) +
        geom_point() +
        geom_errorbar(aes(ymin= fit-se, ymax= fit+se), width=0.4) +
        theme(text = element_text(size=20))+
        #theme_bw(base_size=12) +
        labs(title = "Adult Percent Global DNA Mehtylation 
             by Maternal Rank") +
        theme(plot.title = element_text(hjust = 0.5))+
        ylab("% Global DNA Methylation ± SE") +
        xlab("Maternal Rank")
      
    ## g) Save Plot
      # use ggsave to save the linearization plot
      ggsave("mat_rank_adult_mod_beta.pdf", plot = last_plot(), device = NULL,
             path = "./output/output_luma_ecolog", scale = 1, width = 7,
             height = 5.5,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## h) Do a post-hoc test to determine if Q2, Q3, and Q4 differ
      pairwise.t.test(luma_data_adult$methylation,
                      luma_data_adult$mom.strank.quart,
                      p.adj = "none")
      TukeyHSD(aov(methylation ~ mom.strank.quart + sex,
                   data = luma_data_adult))
      
    ## i) Combine quartiles 2-4 into a single category
      luma_data_adult$mom.strank.quart.comb <- as.factor(
        ifelse(luma_data_adult$mom.strank.quart.order <= 1,
               "Q1 (lowest)", "Q2-Q4 (highest)"))
      
      
    ## j)  Mehtylation by Maternal Rank Quartiles binned
      adult.rank.mod2 <- glm(methylation ~ mom.strank.quart.comb + sex + 
                             age.mon,
                           data = luma_data_adult)
      
      summary(adult.rank.mod2)  # print model summary, effects and SE
      confint(adult.rank.mod2)  # print 95% CIs for parameter estimates
      
      
  ### 9.2 Adult model: methylation by intra litter rank      
    ## a) Check within strata descritpive stats
      luma_data_adult %>%
        group_by (lit.size) %>%
        summarise (n.id = n(),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Mehtylation by intra litter rank
      adult.intra.rank.mod <- glm(methylation ~ lit.size + sex + 
                                  age.mon, data = luma_data_adult)
    ## c) Parameter estimates
      summary(adult.intra.rank.mod)  # print model summary, effects and SE
      confint(adult.intra.rank.mod)  # print 95% CIs for parameter estimates  
      
    ## d) Combine quartiles tiwn hi and lo into a single category
      luma_data_adult$lit.size.comb <- as.factor(
        ifelse(luma_data_adult$lit.size == "single",
               "single", "twin"))
      
    ## e)  Mehtylation by litter size
      adult.intra.rank.mod2 <- glm(methylation ~ lit.size.comb  + sex + 
                                   age.mon,
                                 data = luma_data_adult)
      
      summary(adult.intra.rank.mod2)  # print model summary, effects and SE
      confint(adult.intra.rank.mod2)  # print 95% CIs for parameter estimates
      
  
  ### 9.3 Adult model: methylation by human population size     
    ## a) Check within strata descritpive stats
      luma_data_adult %>%
        group_by (hum.distrb) %>%
        summarise (n.id = n(),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Mehtylation by intra litter rank
      adult.hum.distrb.mod <- glm(methylation ~ hum.distrb + sex + 
                               age.mon, data = luma_data_adult)
    ## c) Parameter estimates
      summary(adult.hum.distrb.mod)  # print model summary, effects and SE
      confint(adult.hum.distrb.mod)  # print 95% CIs for parameter estimates  
      
      
  ### 9.4 Adult model: methylation by prey density     
    ## a) Mehtylation by periconceptional prey density
      adult.peri.prey.mod <- glm(methylation ~ total.peri.concpt + sex + 
                                 age.mon, data = luma_data_adult)
    ## b) Parameter estimates
      summary(adult.peri.prey.mod)  # print model summary, effects and SE
      confint(adult.peri.prey.mod)  # print 95% CIs for parameter estimates
      
    ## c) Mehtylation by gestational prey density
      adult.gest.prey.mod <- glm(methylation ~ total.gest + total.peri.concpt +
                                 sex +  age.mon, data = luma_data_adult)
    ## d) Parameter estimates
      summary(adult.gest.prey.mod)  # print model summary, effects and SE
      confint(adult.gest.prey.mod)  # print 95% CIs for parameter estimates  
      
    ## e) Mehtylation by birth to 3 months prey density
      adult.birth.3.prey.mod <- glm(methylation ~ total.birth.3 + total.gest + 
                                      total.peri.concpt + sex + age.mon, 
                                    data = luma_data_adult)
    ## f) Parameter estimates
      summary(adult.birth.3.prey.mod)  # print model summary, effects and SE
      confint(adult.birth.3.prey.mod)  # print 95% CIs for parameter estimates
      
    ## g) Mehtylation by 3 to 6 months prey density
      adult.3.6.prey.mod <- glm(methylation ~ total.3.6 + total.birth.3 +
                                  total.gest + total.peri.concpt + 
                                  sex + age.mon, data = luma_data_adult)
    ## h) Parameter estimates
      summary(adult.3.6.prey.mod)  # print model summary, effects and SE
      confint(adult.3.6.prey.mod)  # print 95% CIs for parameter estimates
      
    ## i) Mehtylation by 3 to 6 months prey density
      adult.6.9.prey.mod <- glm(methylation ~ total.6.9 + total.3.6 + 
                                total.birth.3 + total.gest + total.peri.concpt 
                                + sex + age.mon, data = luma_data_adult)
      
    ## j) Parameter estimates
      summary(adult.6.9.prey.mod)  # print model summary, effects and SE
      confint(adult.6.9.prey.mod)  # print 95% CIs for parameter estimates    
 
           
#******************************************************************************#      
  ### 9.5 Adult interaction  model: mom rank by prey density   
    ## a) Mehtylation by periconceptional prey density
      adult.peri.x.mom.rank.mod <- glm(methylation ~ 
                                       mom.strank.quart.comb * 
                                       total.peri.concpt + 
                                       sex + age.mon, data = luma_data_adult)
      # Parameter estimates
      summary(adult.peri.x.mom.rank.mod)# print model summary, effects and SE
      confint(adult.peri.x.mom.rank.mod)# print 95% CIs for parameter estimates
      
    ## b) Mehtylation by gestation prey density
      adult.gest.x.mom.rank.mod <- glm(methylation ~ 
                                       mom.strank.quart.comb * 
                                       total.gest + total.peri.concpt +
                                       sex + age.mon, data = luma_data_adult)
      # Parameter estimates
      summary(adult.gest.x.mom.rank.mod)  # print model summary, effects and SE
      confint(adult.gest.x.mom.rank.mod)  # print 95% CIs for parameter estimates
      
    ## c) Mehtylation by birth to 3 month prey density
      adult.birth.3.x.mom.rank.mod <- glm(methylation ~ 
                                          mom.strank.quart.comb * 
                                          total.birth.3 + total.peri.concpt + 
                                          total.gest +
                                          sex + age.mon, data = luma_data_adult)
      # Parameter estimates
      summary(adult.birth.3.x.mom.rank.mod) 
      confint(adult.birth.3.x.mom.rank.mod)
      
    ## d) Mehtylation by 3 to 6 month prey density
      adult.3.6.x.mom.rank.mod <- glm(methylation ~ 
                                      mom.strank.quart.comb * 
                                      total.3.6 + total.peri.concpt + 
                                      total.gest + total.birth.3 +
                                      sex + age.mon, data = luma_data_adult)
      # Parameter estimates
      summary(adult.3.6.x.mom.rank.mod) # print model summary, effects and SE
      confint(adult.3.6.x.mom.rank.mod) # print 95% CIs for parameter estimates
      
    ## e) Mehtylation by 6 to 9 month prey density
      adult.6.9.x.mom.rank.mod <- glm(methylation ~ 
                                      mom.strank.quart.comb * 
                                      total.6.9 + total.peri.concpt + 
                                      total.gest + total.birth.3 + total.3.6 +
                                      sex + age.mon, data = luma_data_adult)
      # Parameter estimates
      summary(adult.6.9.x.mom.rank.mod) # print model summary, effects and SE
      confint(adult.6.9.x.mom.rank.mod) # print 95% CIs for parameter estimates
      
  ### 9.6 adult model: litter size by human pop   
    ## a) Mehtylation by periconceptional prey density
      hum.distrb.x.mom.rank.mod <- glm(methylation ~ 
                                      mom.strank.quart.comb * hum.distrb + 
                                      sex + age.mon, data = luma_data_adult)
      # Parameter estimates
      summary(hum.distrb.x.mom.rank.mod)  # print model summary, effects and SE
      confint(hum.distrb.x.mom.rank.mod)  # print 95% CIs for parameter estimates

           








      
      
      
      
      
   