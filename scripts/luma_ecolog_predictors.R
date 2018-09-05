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
    # 7: Re-tidy data for analyses
    # 8: Cub models
    # 9: Subadult models
    # 10: Adult models



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
        
      # Check for car and install if not already installed
        if (!'car' %in% installed.packages()[,1]){
          install.packages ('car')
        }
      # load car packages
        library ('car')


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
      source(paste0("/Volumes/Holekamp/code_repository/R/1_output_tidy_tbls/",
                   "load_tidy_tbls.R"))
      
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

      
    ## d) Create an estimated age in months by subtracting birthdate from
      # grad.date and weandate using lubridate
      luma_data <- luma_data %>%
        mutate(grad.age.mon = interval(birthdate, den.grad) %/% months(1))
      # and check what mean age at graduation
          luma_data %>%
            summarise(avg.grad.age = round(mean (grad.age.mon,
                                                  na.rm = T), 2),
                      max.grad.age = max(grad.age.mon, na.rm = T))
  
        luma_data <- luma_data %>%
            mutate(wean.age.mon = interval(birthdate, weaned) %/% months(1))
      # and check what mean age at graduation
        luma_data %>%
          summarise(avg.wean.age = round(mean(wean.age.mon,
                                               na.rm = T), 2),
                    max.wean.age = max(wean.age.mon, na.rm = T))    
      
        
    ## e) Create an estimated age in months by subtracting birthdate from
      # darting.date using lubridate
      luma_data <- luma_data %>%
        mutate(age.mon = interval(birthdate, darting.date) %/% months(1))
      
    ## f) Cobmine age columns to fill in NA
      # replace NA in age.months column where there is a value in 
      # estimated.age.mo column
      luma_data$age.months <- ifelse(is.na(luma_data$age.months),
                                     luma_data$estimated.age.mo,
                                     luma_data$age.months)
     
      # replace NA in age.mon column where there is a value in 
      # age.months column
      luma_data$age.mon <- ifelse(is.na(luma_data$age.mon),
                                   luma_data$age.months,
                                   luma_data$age.mon)
      
    ## g) fix age transcription inconsistencies (e.g. all 'sub' to 'subadult') 
      luma_data <- luma_data %>%
        mutate(age = replace(age, str_detect(age, "s"), "subadult")) 
      
    ## h) Create a categorical age variable based age.mon 
      # darting.date using lubridate
      # based on Holekamp and Smale 1998
      luma_data <- luma_data %>%
        mutate(age.cat = case_when(sex == "m" & age.mon <= 12 ~ c("cub"),
                                   sex == "m" & age.mon > 12 & 
                                            age.mon <=24 ~ c("subadult"),
                                   sex == "m" & age.mon > 24 ~ c("adult"),
                                   sex == "f" & age.mon <= 12 ~ c("cub"),
                                   sex == "f" & age.mon > 12 & 
                                     age.mon <=36 ~ c("subadult"),
                                   sex == "f" & age.mon > 36 ~ c("adult")))
      
    ## i) replace NA in age.cat column where there is a value in 
      # age column
      # first convert age to character
      luma_data$age <- as.character(luma_data$age)
      luma_data$age.cat <- ifelse(is.na(luma_data$age.cat),
                                         luma_data$age,
                                         luma_data$age.cat)
      
    ## j) drop old age columns
      luma_data <- luma_data %>%
        select (- c(age)) 
      
      luma_data <- luma_data %>%
        select (- c(age.months))
      
      luma_data <- luma_data %>%
        select (- c(estimated.age.mo)) 
   
    ## k) extract the year for date of interest (here Birthdate) using lubridate
      # and make a new variable
      luma_data$rank_year <- year(as.Date(luma_data$birthdate, 
                                          format="%Y-%m-%d"))
      
    ## l) Make a new variable hum.pres
      # create a 3-level ordinal factor indicating human pastoralist presence 
      # (mid - talek hyena born after 2000, mid - fig tree born after 2000,
      # and lo - all other hynenas for which clan is known and birthdate)
      # based Green et. al 2017
      luma_data <- luma_data  %>%
        mutate(hum.pres = case_when(luma_data$clan %in% c("talek") &
                                        luma_data$rank_year > 2008
                                      ~ c("hi"),
                                      luma_data$clan %in% c("fig.tree") &
                                        luma_data$rank_year > 2008
                                      ~ c("hi"),
                                      luma_data$clan %in% c("talek") &
                                        luma_data$rank_year >= 2000 & 
                                        luma_data$rank_year <= 2008 
                                      ~ c("med"),
                                      luma_data$clan %in% c("fig.tree") &
                                        luma_data$rank_year >= 2000 &
                                        luma_data$rank_year <= 2008 
                                      ~ c("med"),
                                      luma_data$clan %in% c("talek", 
                                                            "fig.tree") &
                                        luma_data$rank_year < 2000 
                                      ~ c("low"),
                                      luma_data$clan %in% c("serena.s", 
                                                            "serena.s",
                                                            "happy.zebra",
                                                            "mara.river")
                                      ~ c("low")))
          
    ## m) Make a new 2 level nomial variable lit.size
      luma_data <- luma_data  %>%
        mutate(lit.size = case_when(number.littermates == 0 ~ c("single"),
                                          litrank == 1 ~ c("twin"),
                                          litrank == 2 ~ c("twin")))

    ## n) extract the year for sample collection data using lubridate
      # and make a new variable  
      luma_data$samp_year <- year(as.Date(luma_data$darting.date, 
                                          format="%Y-%m-%d"))
      
      
  ### 3.2 Tidy tblFemalerank
    ## a) Pattern recognize numbers from Year variable and copy; gets rid of
      # unwanted text characters
      tblFemalerank$rank_year <- as.numeric(regmatches(tblFemalerank$year, 
                                           gregexpr("[[:digit:]]+",
                                                    tblFemalerank$year)))
      
    ## b) rename 'id' variable as 'mom' 
      tblFemalerank <- rename_(tblFemalerank, "mom" = "id")
     

  ### 3.3 Left join tblFemalerank to luma_data  
    ## a) append to luma_data each hyena's mom's rank from the year that 
      # they were born 
      luma_data <- sqldf("SELECT
                         luma_data.*           
                         , tblFemalerank. absrank, strank
                         FROM luma_data      
                         LEFT JOIN tblFemalerank       
                         ON tblFemalerank.mom = luma_data.mom
                         AND tblFemalerank.rank_year = luma_data.rank_year") 
      
    ## b) Manual data clean up
      # mchl's mom, ele, has two ranks that change in the year 2012
      # use rank from jul-jan because mchl born in nov. and discard other row
      luma_data <- luma_data %>%
        filter (!grepl("mchl", id) | !grepl("^6$", absrank))
  
      
  ### 3.4 Tidy joined table 
    ## a) rename 'absrank' variable as 'mom.absrank' 
      luma_data <- rename_(luma_data, "mom.absrank" = "absrank")
      
    ## b) rename 'strank' variable as 'mom.strank' 
      luma_data <- rename_(luma_data, "mom.strank" = "strank") 
      
    ## c) convert mom's strank to numeric
      luma_data$mom.strank <- as.numeric(luma_data$mom.strank)
      
    ## d) Create quartiles of maternal rank
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

    ## e) Create a nominal factor and rename and re-order the levels to 
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
                     hum.pres = (first(hum.pres)),
                     lit.size = (first(lit.size)),
                     age.mon = (mean(age.mon)),
                     rank_year = (first(rank_year)),
                     samp_year = (first(samp_year)),
                     hum.pres = (first(hum.pres)))
                    

                     
    ## b) Ungroup the data frame
        # dplyr retains grouping after creation of data frame that uses 
        # group_by
          luma_data_group <- ungroup(luma_data_group)
        
    
  ### 3.6 Clean prey density data and combine with LUMA data
    ## a) Select a subset of the prey_density data by column names
          
          prey_density <- prey_density %>%
            select(grep("thomsons", names(prey_density)), # contains 'thomsons'
                   grep("topi", names(prey_density)), # contains 'topi'
                   grep("gnu", names(prey_density)), # contains 'gnu'
                   grep("zebra", names(prey_density)), # contains 'zebra'
                   grep("num", names(prey_density)), # contains 'num'
                   grep("^id$", names(prey_density))) %>% # exact match 'id'
            select (-c(number.littermates)) %>%
            select (-starts_with("thom.topi.")) %>%
            select (-starts_with("gnu.zebra."))
          
#        prey_density <- prey_density %>%
#          select(grep("tot", names(prey_density)), # contains 'tot'
#                 grep("thomsons", names(prey_density)), # contains 'thomsons'
#                 grep("topi", names(prey_density)), # contains 'topi'
#                 grep("gnu", names(prey_density)), # contains 'gnu'
#                 grep("zebra", names(prey_density)), # contains 'zebra'
#                 grep("num", names(prey_density)), # contains 'num'
#                 grep("^id$", names(prey_density))) %>% # exact match 'ID'
#          select (-c(number.littermates)) %>%
#          rename ("total.birth.3" = "total.birth-3") %>%
#          rename ("total.3.6" = "total.3-6") %>%
#          rename ("total.6.9" = "total.6-9")
          
    ## b) Left join prey_density to luma_data   
        luma_data <- left_join(luma_data,
                               prey_density, by = "id")
        
    ## c) Left join prey_density to luma_data_group   
        luma_data_group <- left_join(luma_data_group,
                              prey_density, by = "id")
  
  
  ## 3.7 Re-code variables to the appropriate class
    ## a) Re-code *nominal* factor (with ordered levels)  
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
        
        luma_data_group <- transform(luma_data_group, 
                               age.cat = factor(age.cat,
                                                levels = c("cub", 
                                                           "subadult", 
                                                           "adult")))
        
    ## b) Re-code *interval* factor (with ordered levels)  
        # Set order and *equally* spaced categorical age variable a
        # NOTE: model output a p for trend
        # polynomial contrasts are linear, quadratic, cubic etc. estimates
        # stored internally as 1, 2, 3 etc. (equal spacing)
        luma_data <- transform(luma_data, 
                               age.inter = ordered(age.cat,
                                                   levels = c("cub",
                                                              "subadult",
                                                              "adult")))
        
    ## c) Re-code *ordinal* factor (with ordered levels)  
        # Set order and *unequally* spaced categorical age variable 
        # (e.g. median of each category)
        # NOTE: model output a p for trend; assumes linear monotonic association
        # estimates, which is a change y with respect to each unit change in x
        # is scaled to x and stored internally as numeric value
        luma_data <- luma_data %>%
          mutate(age.ordin = as.numeric(case_when(age.mon <= 12 ~ c(11.0),
                                                  age.mon > 12 & 
                                                    age.mon <=24 ~ c(17.5),
                                                  age.mon > 24 ~ c(51.0))))
        
    ## d) Re-code sex as a nominal factor  
        luma_data$sex <- as.factor(luma_data$sex)
        
        luma_data_group$sex <- as.factor(luma_data_group$sex)
        
    ## e) Re-code hum.pres as nominal factor and set level (order)
        luma_data <- transform(luma_data,
                               hum.pres = factor(hum.pres,
                                                   levels = c("low", "med", 
                                                              "hi")))
        luma_data_group <- transform(luma_data_group,
                               hum.pres = factor(hum.pres,
                                                   levels = c("low", "med", 
                                                               "hi")))
        
    ## f) Re-code lit.size as nominal factor and set level (order)
        luma_data <- transform(luma_data,
                               lit.size = factor(lit.size, 
                                                 levels = c("single","twin")))
        
        luma_data_group <- transform(luma_data_group,
                               lit.size = factor(lit.size, 
                                                 levels = c("single","twin")))
        
             
        
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
                                        na.rm = T), 2),
                     med_samp_age = round (quantile (samp_year, 
                                                     c(.5), na.rm = T), 2))
     
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
        summarise (n.age = sum(!is.na(age.cat)),
                   avg.age = round (mean(age.mon, na.rm = T), 2),
                   med.age = round (median(age.mon, na.rm = T), 2),
                   stdev.age = round (sd(age.mon, na.rm = T), 2)) %>%
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
        summarise (n.mom.rank = sum(!is.na(mom.strank.quart)),
                   avg.mom.rank = round (mean(mom.strank, na.rm = T),2),
                   stdev.mom.rank = round (sd(mom.strank, na.rm = T), 2)) %>%
        mutate (freq.mom.rank =  (n.mom.rank / sum(n.mom.rank)))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/output_luma_ecolog/mom_rank_summary.pdf"), 
          height = 4, width = 8)
      grid.table(mom_rank_summary)
      dev.off() 
      
      
  ### 4.5 Descriptive stats litter size
    ## a) Intra litter rank summary 
      lit_size_summary <- luma_data_group %>%
        #group_by (id) %>%
        group_by (lit.size) %>%
        summarise(n = sum(!is.na(lit.size))) %>%
        mutate(freq = n / sum(n))
      
      ## b) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/output_luma_ecolog/",
                 "intra_lit_rank_summary.pdf"), 
          height = 3, width = 5)
      grid.table(lit_size_summary)
      dev.off() 
   
        
  ### 4.6 Descriptive stats human presence
      ## a) Human presence (proxy) summary 
      hum_pres_summary <- luma_data_group %>%
        #group_by (id) %>%
        group_by (hum.pres) %>%
        summarise(n  =sum(!is.na(hum.pres))) %>%
        mutate(freq = n / sum(n))
      
      ## b) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/output_luma_ecolog/",
                 "hum_pres_summary.pdf"), 
          height = 3, width = 5)
      grid.table(hum_pop_summary)
      dev.off() 
      
      
  ### 4.7 Descriptive stats prey density            
    ## a) Prey density summary 
      prey_peri_summary <- luma_data_group %>%
        summarise (dev.period = print("peri.concpt"),
                   n.prey = sum(!is.na(thomsons.peri.concpt)),
                   avg.thomsons = round (mean(thomsons.peri.concpt, 
                                               na.rm = T),2),
                   stdev.thomsons = round (sd(thomsons.peri.concpt, 
                                               na.rm = T), 2),
                   avg.topi = round (mean(topi.peri.concpt, 
                                         na.rm = T),2),
                   stdev.topi = round (sd(topi.peri.concpt, 
                                         na.rm = T), 2),
                   avg.gnu = round (mean(gnu.peri.concpt, 
                                          na.rm = T),2),
                   stdev.gnu = round (sd(gnu.peri.concpt, 
                                          na.rm = T), 2),
                   avg.zebra = round (mean(zebra.peri.concpt, 
                                         na.rm = T),2),
                   stdev.zebra = round (sd(zebra.peri.concpt, 
                                         na.rm = T), 2))
                   
      prey_gest_summary <-  luma_data_group %>%
        summarise (dev.period = print("gest"),
                   n.prey = sum(!is.na(thomsons.gest)),
                   avg.thomsons = round (mean(thomsons.gest, 
                                              na.rm = T),2),
                   stdev.thomsons = round (sd(thomsons.gest, 
                                              na.rm = T), 2),
                   avg.topi = round (mean(topi.gest, 
                                          na.rm = T),2),
                   stdev.topi = round (sd(topi.gest, 
                                          na.rm = T), 2),
                   avg.gnu = round (mean(gnu.gest, 
                                         na.rm = T),2),
                   stdev.gnu = round (sd(gnu.gest, 
                                         na.rm = T), 2),
                   avg.zebra = round (mean(zebra.gest, 
                                           na.rm = T),2),
                   stdev.zebra = round (sd(zebra.gest, 
                                           na.rm = T), 2))
      
      prey_birth.3_summary <- luma_data_group %>%
        summarise (dev.period = print("birth.3"),
                   n.prey = sum(!is.na(thomsons.birth.3)),
                   avg.thomsons = round (mean(thomsons.birth.3, 
                                              na.rm = T),2),
                   stdev.thomsons = round (sd(thomsons.birth.3, 
                                              na.rm = T), 2),
                   avg.topi = round (mean(topi.birth.3, 
                                          na.rm = T),2),
                   stdev.topi = round (sd(topi.birth.3, 
                                          na.rm = T), 2),
                   avg.gnu = round (mean(gnu.birth.3, 
                                         na.rm = T),2),
                   stdev.gnu = round (sd(gnu.birth.3, 
                                         na.rm = T), 2),
                   avg.zebra = round (mean(zebra.birth.3, 
                                           na.rm = T),2),
                   stdev.zebra = round (sd(zebra.birth.3, 
                                           na.rm = T), 2))
                   
      prey_3.6_summary <- luma_data_group %>%
        summarise (dev.period = print("3.6"),
                   n.prey = sum(!is.na(thomsons.3.6)),
                   avg.thomsons = round (mean(thomsons.3.6, 
                                              na.rm = T),2),
                   stdev.thomsons = round (sd(thomsons.3.6, 
                                              na.rm = T), 2),
                   avg.topi = round (mean(topi.3.6, 
                                          na.rm = T),2),
                   stdev.topi = round (sd(topi.3.6, 
                                          na.rm = T), 2),
                   avg.gnu = round (mean(gnu.3.6, 
                                         na.rm = T),2),
                   stdev.gnu = round (sd(gnu.3.6, 
                                         na.rm = T), 2),
                   avg.zebra = round (mean(zebra.3.6, 
                                           na.rm = T),2),
                   stdev.zebra = round (sd(zebra.3.6, 
                                           na.rm = T), 2))
                   
      prey_6.9_summary <- luma_data_group %>%
        summarise (dev.period = print("6.9"),
                   n.prey = sum(!is.na(thomsons.6.9)),
                   avg.thomsons = round (mean(thomsons.6.9, 
                                              na.rm = T),2),
                   stdev.thomsons = round (sd(thomsons.6.9, 
                                              na.rm = T), 2),
                   avg.topi = round (mean(topi.6.9, 
                                          na.rm = T),2),
                   stdev.topi = round (sd(topi.6.9, 
                                          na.rm = T), 2),
                   avg.gnu = round (mean(gnu.6.9, 
                                         na.rm = T),2),
                   stdev.gnu = round (sd(gnu.6.9, 
                                         na.rm = T), 2),
                   avg.zebra = round (mean(zebra.6.9, 
                                           na.rm = T),2),
                   stdev.zebra = round (sd(zebra.6.9, 
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
        center_fxn_cols <- function(x) {
          xcenter = colMeans(x, na.rm = T)
          x - rep(xcenter, rep.int(nrow(x), ncol(x)))
        }
  
    ## b) Center samp_year
      luma_data$samp_year_cnt <- as.numeric(scale(luma_data$samp_year, 
                                                  scale = F))
      luma_data_group$samp_year_cnt <- as.numeric(scale
                                                  (luma_data_group$samp_year,
                                                    scale = F))

    ## c) make a list of variable names to center, value = T is necessary
        # or column positions will be returned
        var_names <- c("total","thomsons", "topi", "gnu", "zebra")
        vars_to_center <- grep(paste(var_names, collapse = "|"), 
                               names(luma_data_group), value = T)
        
        # drop total solids variable
        vars_to_center <- vars_to_center[vars_to_center != "total.solids"]
     
    ## d) Z-score standardize
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
  
  ### 6.1 Bivariate statistics methylation by sample age
    # NOTE: uses luma_data_group; first average over ID within age cat. 
    # and then take avg 
    # uses 'nmle' package, which will provided p-value estimates
      samp.year.lme <- lme(methylation ~ samp_year_cnt, random =~1|id, 
                            subset(luma_data_group,!is.na(x = samp_year_cnt)))
        
      summary(samp.year.lme)    #  print model summary, effects and SE
      intervals(samp.year.lme, 
                  which = "fixed")  # print 95% CIs for parameter estimates
      
      
  ### 6.2 Bivariate statistics methylation by sex
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
    ## a) Graph of the raw data for percent global DNA methylaiton by maternal 
      # rank stratified by age
      ggplot(data = subset(luma_data_group, !is.na(x = mom.strank)),
             aes(x = mom.strank, y = methylation)) +
        geom_point(shape = 1) +
        theme(text = element_text(size=20))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        geom_smooth(method = loess, se = F) + # Add smooth curve best fit lines
        labs(title = "Percent Global DNA Mehtylation by 
             Maternal Rank",
             fill = "age") +
        theme(plot.title = element_text(hjust = 0.5))+
        ylab("% Global DNA Methylation") +
        xlab("Maternal Rank")
      
    ## b) Save Plot
      # use ggsave to save the linearization plot
      ggsave("meth_by_rank_loess_plot.pdf", plot = last_plot(), 
             device = NULL,
             path = paste0(here(),"/output/output_luma_ecolog"), 
             scale = 1, width = 7, height = 5,
             units = c("in"), dpi = 300, limitsize = TRUE)
    
    ## c) Graph of the raw data for percent global DNA methylaiton by maternal 
      # rank stratified by age because age seems to be important (see above)
      ggplot(data = subset(luma_data_group, !is.na(x = mom.strank)),
             aes(x = mom.strank, y = methylation, color = age.cat)) +
        geom_point(shape = 1) +
        theme(text = element_text(size=20))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        geom_smooth(method = loess, se = F) + # Add smooth curve best fit lines
        labs(title = "Percent Global DNA Mehtylation by 
             Maternal Rank",
             fill = "age") +
        theme(plot.title = element_text(hjust = 0.5))+
        ylab("% Global DNA Methylation") +
        xlab("Maternal Rank")
      
      ## d) Save Plot
      # use ggsave to save the linearization plot
      ggsave("meth_by_rank_by_age_loess_plot.pdf", plot = last_plot(), 
             device = NULL,
             path = paste0(here(),"/output/output_luma_ecolog"), 
             scale = 1, width = 7, height = 5,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## e) Summary stats methylation by age 
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg
      meth_by_mom_rank <- luma_data_group %>%
        group_by (mom.strank.quart) %>%
        summarise (n.id = n(),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))

    ## f) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/output_luma_ecolog/meth_by_mom_rank.pdf"), 
          height = 6, width = 7)
      grid.table(meth_by_mom_rank)
      dev.off()
    
    ## g) Plot mehtylation by rank
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
    
    ## h) Save Plot
      # use ggsave to save the plot
      ggsave("meth_by_mom_rank_plot.pdf", plot = last_plot(), device = NULL, 
             path = paste0(here(),"/output/output_luma_ecolog"),
             scale = 1, width = 7, height = 5, 
             units = c("in"), dpi = 300, limitsize = TRUE)  
      
    ## i) Bivariate Regression Methylatino by maternal rank 
      # uses 'nmle' package, which will provided p-value estimates
      mom.rank.lme <- lme(methylation ~ mom.strank.quart, random =~1|id, 
                     subset(luma_data_group,!is.na(x = mom.strank.quart)))
      
      summary(mom.rank.lme)       #  print model summary, effects and SE
      intervals(mom.rank.lme, 
                which = "fixed")  # print 95% CIs for parameter estimates
      anova(mom.rank.lme)         # generate p-value from Wald test
  


  ### 6.4 Bivariate statistics methylation by litter size
    ## a) Summary stats methylation by lit.size 
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg
      meth_by_lit_size <- luma_data_group %>%
        group_by (lit.size) %>%
        summarise (n.id = n(),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/output_luma_ecolog/",
                 "meth_by_lit_size.pdf"), 
          height = 6, width = 7)
      grid.table(meth_by_lit_size)
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
Mehtylation by Litter Size") +
        ylab("% Global DNA Methylation") +
        xlab("Litter Size")
      
      ## d) Save Plot
      # use ggsave to save the plot
      ggsave("meth_by_lit_size_plot.pdf", plot = last_plot(), 
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
    ## a) Summary stats methylation by hum.pres
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg
      meth_by_hum_pop <- luma_data_group %>%
        group_by (hum.pres) %>%
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
      ggplot(data = subset(luma_data_group, !is.na(x = hum.pres)),
             aes(x = hum.pres, y = methylation,
                 color = hum.pres)) + 
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
      hum.pres.lme <- lme(methylation ~ hum.pres, random =~1|id, 
                                subset(luma_data_group,!is.na(x = hum.pres)))
      
      summary(hum.pres.lme)        #  print model summary, effects and SE
      intervals(hum.pres.lme, 
                which = "fixed")  # print 95% CIs for parameter estimates
      anova(hum.pres.lme)          # generate p-value from Wald test
      
    
  ### 6.6 Bivariate statistics methylation by periconceptional prey density
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg 
      # uses 'nmle' package, which will provided p-value estimates
      thomsons.peri.concpt.lme <- lme(methylation ~ thomsons.peri.concpt, 
                                      random =~1|id, 
                               subset(luma_data_group,
                                      !is.na(x = thomsons.peri.concpt)))
      
      summary(thomsons.peri.concpt.lme)   # print model summary
      intervals(thomsons.peri.concpt.lme, 
                which = "fixed")  # print 95% CIs
      
      topi.peri.concpt.lme <- lme(methylation ~ topi.peri.concpt, 
                                      random =~1|id, 
                                      subset(luma_data_group,
                                             !is.na(x = topi.peri.concpt)))
      
      summary(topi.peri.concpt.lme)   # print model summary
      intervals(topi.peri.concpt.lme, 
                which = "fixed")  # print 95% CIs
      
      gnu.peri.concpt.lme <- lme(methylation ~ gnu.peri.concpt, 
                                      random =~1|id, 
                                      subset(luma_data_group,
                                             !is.na(x = gnu.peri.concpt)))
      
      summary(gnu.peri.concpt.lme)   # print model summary
      intervals(gnu.peri.concpt.lme, 
                which = "fixed")  # print 95% CIs
      
      zebra.peri.concpt.lme <- lme(methylation ~ zebra.peri.concpt, 
                                      random =~1|id, 
                                      subset(luma_data_group,
                                             !is.na(x = zebra.peri.concpt)))
      
      summary(zebra.peri.concpt.lme)   # print model summary
      intervals(zebra.peri.concpt.lme, 
                which = "fixed")  # print 95% CIs
    
  ### 6.7 Bivariate statistics methylation by gestational prey density
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg 
      # uses 'nmle' package, which will provided p-value estimates
      thomsons.gest.lme <- lme(methylation ~ thomsons.gest, 
                                      random =~1|id, 
                                      subset(luma_data_group,
                                             !is.na(x = thomsons.gest)))
      
      summary(thomsons.gest.lme)   # print model summary
      intervals(thomsons.gest.lme, 
                which = "fixed")  # print 95% CIs
      
      topi.gest.lme <- lme(methylation ~ topi.gest, 
                                  random =~1|id, 
                                  subset(luma_data_group,
                                         !is.na(x = topi.gest)))
      
      summary(topi.gest.lme)   # print model summary
      intervals(topi.gest.lme, 
                which = "fixed")  # print 95% CIs
      
      gnu.gest.lme <- lme(methylation ~ gnu.gest, 
                                 random =~1|id, 
                                 subset(luma_data_group,
                                        !is.na(x = gnu.gest)))
      
      summary(gnu.gest.lme)   # print model summary
      intervals(gnu.gest.lme, 
                which = "fixed")  # print 95% CIs
      
      zebra.gest.lme <- lme(methylation ~ zebra.gest, 
                                   random =~1|id, 
                                   subset(luma_data_group,
                                          !is.na(x = zebra.gest)))
      
      summary(zebra.gest.lme)   # print model summary
      intervals(zebra.gest.lme, 
                which = "fixed")  # print 95% CIs
      
      
  ### 6.8 Bivariate statistics methylation by birth to 3 months prey density
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg
      thomsons.birth.3.lme <- lme(methylation ~ thomsons.birth.3, 
                                      random =~1|id, 
                                      subset(luma_data_group,
                                             !is.na(x = thomsons.birth.3)))
      
      summary(thomsons.birth.3.lme)   # print model summary
      intervals(thomsons.birth.3.lme, 
                which = "fixed")  # print 95% CIs
      
      topi.birth.3.lme <- lme(methylation ~ topi.birth.3, 
                                  random =~1|id, 
                                  subset(luma_data_group,
                                         !is.na(x = topi.birth.3)))
      
      summary(topi.birth.3.lme)   # print model summary
      intervals(topi.birth.3.lme, 
                which = "fixed")  # print 95% CIs
      
      gnu.birth.3.lme <- lme(methylation ~ gnu.birth.3, 
                                 random =~1|id, 
                                 subset(luma_data_group,
                                        !is.na(x = gnu.birth.3)))
      
      summary(gnu.birth.3.lme)   # print model summary
      intervals(gnu.birth.3.lme, 
                which = "fixed")  # print 95% CIs
      
      zebra.birth.3.lme <- lme(methylation ~ zebra.birth.3, 
                                   random =~1|id, 
                                   subset(luma_data_group,
                                          !is.na(x = zebra.birth.3)))
      
      summary(zebra.birth.3.lme)   # print model summary
      intervals(zebra.birth.3.lme, 
                which = "fixed")  # print 95% CIs
      
  
  ### 6.9 Bivariate statistics methylation by 3 to 6 months prey density
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg 
      # uses 'nmle' package, which will provided p-value estimates
      thomsons.3.6.lme <- lme(methylation ~ thomsons.3.6, 
                                        random =~1|id, 
                                        subset(luma_data_group,
                                               !is.na(x = thomsons.3.6)))
      
      summary(thomsons.3.6.lme)   # print model summary
      intervals(thomsons.3.6.lme, 
                which = "fixed")  # print 95% CIs
      
      topi.3.6.lme <- lme(methylation ~ topi.3.6, 
                                  random =~1|id, 
                                  subset(luma_data_group,
                                         !is.na(x = topi.3.6)))
      
      summary(topi.3.6.lme)   # print model summary
      intervals(topi.3.6.lme, 
                which = "fixed")  # print 95% CIs
      
      gnu.3.6.lme <- lme(methylation ~ gnu.3.6, 
                                 random =~1|id, 
                                 subset(luma_data_group,
                                        !is.na(x = gnu.3.6)))
      
      summary(gnu.3.6.lme)   # print model summary
      intervals(gnu.3.6.lme, 
                which = "fixed")  # print 95% CIs
      
      zebra.3.6.lme <- lme(methylation ~ zebra.3.6, 
                                   random =~1|id, 
                                   subset(luma_data_group,
                                          !is.na(x = zebra.3.6)))
      
      summary(zebra.3.6.lme)   # print model summary
      intervals(zebra.3.6.lme, 
                which = "fixed")  # print 95% CIs
      
  ### 6.10 Bivariate statistics methylation by 6 to 9 months prey density
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg 
      # uses 'nmle' package, which will provided p-value estimates
      thomsons.6.9.lme <- lme(methylation ~ thomsons.6.9, 
                                      random =~1|id, 
                                      subset(luma_data_group,
                                             !is.na(x = thomsons.6.9)))
      
      summary(thomsons.6.9.lme)   # print model summary
      intervals(thomsons.6.9.lme, 
                which = "fixed")  # print 95% CIs
      
      topi.6.9.lme <- lme(methylation ~ topi.6.9, 
                                  random =~1|id, 
                                  subset(luma_data_group,
                                         !is.na(x = topi.6.9)))
      
      summary(topi.6.9.lme)   # print model summary
      intervals(topi.6.9.lme, 
                which = "fixed")  # print 95% CIs
      
      gnu.6.9.lme <- lme(methylation ~ gnu.6.9, 
                                 random =~1|id, 
                                 subset(luma_data_group,
                                        !is.na(x = gnu.6.9)))
      
      summary(gnu.6.9.lme)   # print model summary
      intervals(gnu.6.9.lme, 
                which = "fixed")  # print 95% CIs
      
      zebra.6.9.lme <- lme(methylation ~ zebra.6.9, 
                                   random =~1|id, 
                                   subset(luma_data_group,
                                          !is.na(x = zebra.6.9)))
      
      summary(zebra.6.9.lme)   # print model summary
      intervals(zebra.6.9.lme, 
                which = "fixed")  # print 95% CIs
 
      
      
###############################################################################
##############            7. Re-tidy daa for analyses            ##############
###############################################################################       
      
         
  ### 7.1 View methylation by rank within age strata
    ## a) Graph of the raw data for percent global DNA methylaiton by maternal rank
      # stratified by age
      ggplot(subset(luma_data_group, !is.na(x = mom.strank)),
             aes(x = mom.strank.quart.order, y = methylation, 
                 color = age.cat )) +
        geom_point(shape = 1) +
        theme(text = element_text(size=20))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        geom_smooth(method = loess, se = F) + # Add linear regression best fit lines
        labs(title = "Percent Global DNA Mehtylation by 
             Maternal Rank by Age",
             fill = "age") +
        theme(plot.title = element_text(hjust = 0.5))+
        ylab("% Global DNA Methylation") +
        xlab("Maternal Rank")
      
    ## b) Save Plot
      # use ggsave to save the linearization plot
      ggsave("meth_by_rank_quart_by_age_loess_plot.pdf", plot = last_plot(), 
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
      
      
  ### 7.2 Stratify data by life history variables
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
##############                  8. Cub models                    ##############
###############################################################################      

        
  ### 8.1 Cub model: methylation by sex
    ## a) Check within strata descritpive stats
      luma_data_cub %>%
        group_by (sex) %>%
        summarise (n.id = n(),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Unadjusted: methlyation by sex
      cub.sex.unadj <- glm(methylation ~  sex, data = luma_data_cub)
      
    ## c) Parameter estimates
      summary(cub.sex.unadj)  # model parameter estimates
      confint(cub.sex.unadj)  # 95% CIs 

    ## d) Adjusted: methlyation by sex
      cub.sex.adj <- glm(methylation ~  sex + age.mon +
                                samp_year_cnt,
                              data = luma_data_cub)
      
    ## e) Parameter estimates
      summary(cub.sex.adj)  # model parameter estimates
      confint(cub.sex.adj)  # 95% CIs 
    

  ### 8.2 Cub model: methylation by age in months
    ## a) Check within strata descritpive stats
      luma_data_cub %>%
        #group_by () %>%
        summarise (n.id = sum(!is.na(age.mon)),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Unadjusted: methlyation by age.mon
      cub.age.mon.unadj <- glm(methylation ~  age.mon, data = luma_data_cub)
      
    ## c) Parameter estimates
      summary(cub.age.mon.unadj)  # model parameter estimates
      confint(cub.age.mon.unadj)  # 95% CIs 
      
    ## d) Adjusted: methlyation by age.mon
      cub.age.mon.adj <- glm(methylation ~  age.mon + sex +
                               samp_year_cnt,
                             data = luma_data_cub)
      
    ## e) Parameter estimates
      summary(cub.age.mon.adj)  # model parameter estimates
      confint(cub.age.mon.adj)  # 95% CIs 
      

  ### 8.3 Cub model: methylation by mom rank quartiles    
    ## a) Check within strata descritpive stats
      luma_data_cub %>%
        group_by (mom.strank.quart) %>%
        summarise (n.id = sum(!is.na(age.mon)),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Unadjusted: methlyation by mom.strank.quart
      cub.mom.rank.unadj <- glm(methylation ~  mom.strank.quart, 
                               data = luma_data_cub)
      
    ## c) Parameter estimates
      summary(cub.mom.rank.unadj)  # model parameter estimates
      confint(cub.mom.rank.unadj)  # 95% CIs 
      # use 'car' package to extract a Wald Chi-Sq p-value; test for overall 
      # difference
      Anova(cub.mom.rank.unadj, Type ="II", test = "Wald") # Wald test p
      
    ## d) Adjusted: methlyation by mom.strank.quart
      cub.mom.rank.mod.adj <- glm(methylation ~ mom.strank.quart + sex + 
                                    age.mon + samp_year_cnt,
                              data = luma_data_cub)
      
    ## e) Parameter estimates
      summary(cub.mom.rank.adj)  # print model summary, effects and SE
      confint(cub.mom.rank.adj)  # print 95% CIs for parameter estimates
      Anova(cub.mom.rank.adj, Type ="II", test = "Wald") # Wald test p
      
      
    ## f) Extract mom.strank.quart estimates and 
      cub.rank.ef <- effect("mom.strank.quart", cub.mom.rank.adj)
      summary(cub.rank.ef)
      # Save effects as data frame
      cub.rank.ef.table <- as.data.frame(cub.rank.ef)
      # Set the reference level to Q1
      cub.rank.ef.table <- transform( cub.rank.ef.table,
                          mom.strank.quart = factor(mom.strank.quart,
                                                    levels = c("Q1 (lowest)", 
                                                              "Q2", "Q3", 
                                                              "Q4 (highest)")))
      
    ## g) Graph cub.mom.rank effects
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
    
    ## h) Save Plot
      # use ggsave to save the linearization plot
      ggsave("mat_rank_cub_mod_beta.pdf", plot = last_plot(), device = NULL,
             path = "./output/output_luma_ecolog", scale = 1, width = 7,
             height = 5.5,
             units = c("in"), dpi = 300, limitsize = TRUE)
 
    ## i) Do a post-hoc test to determine if Q2, Q3, and Q4 differ
      pairwise.t.test(luma_data_cub$methylation,
                      luma_data_cub$mom.strank.quart,
                      p.adj = "none")
      TukeyHSD(aov(methylation ~ mom.strank.quart + sex,
                   data = luma_data_cub))
    
    ## j) Combine quartiles 2-4 into a single category
      luma_data_cub$mom.strank.quart.comb <- as.factor(
      ifelse(luma_data_cub$mom.strank.quart ==  "Q1 (lowest)",
             "Q1 (lowest)", "Q2-Q4 (highest)"))
    
    ## k)  Adjusted: methlyation by mom.strank.quart binned
      cub.rank.adj2 <- glm(methylation ~ mom.strank.quart.comb + sex + 
                             age.mon + samp_year_cnt,
                         data = luma_data_cub)
      
      summary(cub.rank.adj2)  # print model summary, effects and SE
      confint(cub.rank.mod2)  # print 95% CIs for parameter estimates
      Anova(cub.mom.rank.adj2, Type ="II", test = "Wald") # Wald test p
    
      
  ### 8.4 Cub model: methylation by litter size
    ## a) Check within strata descritpive stats
      luma_data_cub %>%
        group_by (lit.size) %>%
        summarise (n.id = sum(!is.na(lit.size)),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Unadjusted: methlyation by age.mon
      cub.lit.size.unadj <- glm(methylation ~  lit.size, data = luma_data_cub)
      
    ## c) Parameter estimates
      summary(cub.lit.size.unadj)  # model parameter estimates
      confint(cub.lit.size.unadj)  # 95% CIs
      Anova(cub.lit.size.unadj, Type ="II", test = "Wald") # Wald test p
      
    ## d) Adjusted: methlyation by age.mon
      cub.lit.size.adj <- glm(methylation ~  lit.size + age.mon + sex +
                               samp_year_cnt,
                             data = luma_data_cub)
      
    ## e) Parameter estimates
      summary(cub.lit.size.adj)  # model parameter estimates
      confint(cub.lit.size.adj)  # 95% CIs 
      Anova(cub.lit.size.adj, Type ="II", test = "Wald") # Wald test p
      
      
  ### 8.5 Cub model: methylation by human presence proxy
      ## a) Check within strata descritpive stats
      luma_data_cub %>%
        group_by (hum.pres) %>%
        summarise (n.id = sum(!is.na(lit.size)),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
      ## b) Unadjusted: methlyation by age.mon
      cub.hum.pres.unadj <- glm(methylation ~  hum.pres, data = luma_data_cub)
      
      ## c) Parameter estimates
      summary(cub.hum.pres.unadj)  # model parameter estimates
      confint(cub.hum.pres.unadj)  # 95% CIs
      Anova(cub.hum.pres.unadj, Type ="II", test = "Wald") # Wald test p
      
      ## d) Adjusted: methlyation by age.mon
      cub.hum.pres.adj <- glm(methylation ~  hum.pres + age.mon + sex +
                                samp_year_cnt,
                              data = luma_data_cub)
      
      ## e) Parameter estimates
      summary(cub.hum.pres.adj)  # model parameter estimates
      confint(cub.hum.pres.adj)  # 95% CIs 
      Anova(cub.hum.pres.adj, Type ="II", test = "Wald") # Wald test p
    

  ### 8.6 Cub model: methylation by periconceptional prey density     
    ## a) Unadjusted: methlyation by periconceptional thomsons density
      cub.peri.thomsons.unadj <- glm(methylation ~ thomsons.peri.concpt, 
                               data = luma_data_cub)
      
    ## b) Parameter estimates
      summary(cub.peri.thomsons.unadj)  # model parameter estimates
      confint(cub.peri.thomsons.unadj)  # 95% CIs 
      
    ## c) Adjusted: methlyation by periconceptional thomsons density
      cub.peri.thomsons.adj <- glm(methylation ~ thomsons.peri.concpt + sex + 
                                      age.mon + samp_year_cnt, 
                                    data = luma_data_cub)
      
    ## d) Parameter estimates
      summary(cub.peri.thomsons.adj)  # model parameter estimates
      confint(cub.peri.thomsons.adj)  # 95% CIs 
      
    ## e) Unadjusted: methlyation by periconceptional topi density
      cub.peri.topi.unadj <- glm(methylation ~ topi.peri.concpt, 
                                     data = luma_data_cub)
      
    ## f) Parameter estimates
      summary(cub.peri.topi.unadj)  # model parameter estimates
      confint(cub.peri.topi.unadj)  # 95% CIs 
      
    ## g) Adjusted: methlyation by periconceptional topi density
      cub.peri.topi.adj <- glm(methylation ~ topi.peri.concpt + sex + 
                                     age.mon + samp_year_cnt, 
                                   data = luma_data_cub)
      
    ## h) Parameter estimates
      summary(cub.peri.topi.adj)  # model parameter estimates
      confint(cub.peri.topi.adj)  # 95% CIs 
      
    ## i) Unadjusted: methlyation by periconceptional gnu density
      cub.peri.gnu.unadj <- glm(methylation ~ gnu.peri.concpt, 
                                     data = luma_data_cub)
      
    ## j) Parameter estimates
      summary(cub.peri.gnu.unadj)  # model parameter estimates
      confint(cub.peri.gnu.unadj)  # 95% CIs 
      
    ## k) Adjusted: methlyation by periconceptional gnu density
      cub.peri.gnu.adj <- glm(methylation ~ gnu.peri.concpt + sex + 
                                     age.mon + samp_year_cnt, 
                                   data = luma_data_cub)
      
    ## l) Parameter estimates
      summary(cub.peri.gnu.adj)  # model parameter estimates
      confint(cub.peri.gnu.adj)  # 95% CIs 
      
    ## m) Unadjusted: methlyation by periconceptional zebra density
      cub.peri.zebra.unadj <- glm(methylation ~ zebra.peri.concpt, 
                                     data = luma_data_cub)
      
    ## n) Parameter estimates
      summary(cub.peri.zebra.unadj)  # model parameter estimates
      confint(cub.peri.zebra.unadj)  # 95% CIs 
      
    ## o) Adjusted: methlyation by periconceptional zebra density
      cub.peri.zebra.adj <- glm(methylation ~ zebra.peri.concpt + sex + 
                                     age.mon + samp_year_cnt, 
                                   data = luma_data_cub)
      
    ## p) Parameter estimates
      summary(cub.peri.zebra.adj)  # model parameter estimates
      confint(cub.peri.zebra.adj)  # 95% CIs 
  
      
  ### 8.7 Cub model: methylation by gestational prey density     
    ## a) Unadjusted: methlyation by gestational thomsons density
      cub.gest.thomsons.unadj <- glm(methylation ~ thomsons.gest, 
                                     data = luma_data_cub)
      
    ## b) Parameter estimates
      summary(cub.gest.thomsons.unadj)  # model parameter estimates
      confint(cub.gest.thomsons.unadj)  # 95% CIs 
      
    ## c) Adjusted: methlyation by gestational thomsons density
      cub.gest.thomsons.adj <- glm(methylation ~ thomsons.gest + 
                                     thomsons.peri.concpt + sex + 
                                     age.mon + samp_year_cnt, 
                                   data = luma_data_cub)
      
    ## d) Parameter estimates
      summary(cub.gest.thomsons.adj)  # model parameter estimates
      confint(cub.gest.thomsons.adj)  # 95% CIs 
      
    ## e) Unadjusted: methlyation by gestational topi density
      cub.gest.topi.unadj <- glm(methylation ~ topi.peri.concpt, 
                                 data = luma_data_cub)
      
    ## f) Parameter estimates
      summary(cub.gest.topi.unadj)  # model parameter estimates
      confint(cub.gest.topi.unadj)  # 95% CIs 
      
    ## g) Adjusted: methlyation by gestational topi density
      cub.gest.topi.adj <- glm(methylation ~ topi.gest + topi.peri.concpt + sex + 
                                 age.mon + samp_year_cnt, 
                               data = luma_data_cub)
      
    ## h) Parameter estimates
      summary(cub.gest.topi.adj)  # model parameter estimates
      confint(cub.gest.topi.adj)  # 95% CIs 
      
    ## i) Unadjusted: methlyation by gestational gnu density
      cub.gest.gnu.unadj <- glm(methylation ~ gnu.gest, 
                                data = luma_data_cub)
      
    ## j) Parameter estimates
      summary(cub.gest.gnu.unadj)  # model parameter estimates
      confint(cub.gest.gnu.unadj)  # 95% CIs 
      
    ## k) Adjusted: methlyation by gestational gnu density
      cub.gest.gnu.adj <- glm(methylation ~ gnu.gest + gnu.peri.concpt + sex + 
                                age.mon + samp_year_cnt, 
                              data = luma_data_cub)
      
      ## l) Parameter estimates
      summary(cub.gest.gnu.adj)  # model parameter estimates
      confint(cub.gest.gnu.adj)  # 95% CIs 
      
    ## m) Unadjusted: methlyation by gestational zebra density
      cub.gest.zebra.unadj <- glm(methylation ~ zebra.gest,
                                  data = luma_data_cub)
      
    ## n) Parameter estimates
      summary(cub.gest.zebra.unadj)  # model parameter estimates
      confint(cub.gest.zebra.unadj)  # 95% CIs 
      
    ## o) Adjusted: methlyation by gestational zebra density
      cub.gest.zebra.adj <- glm(methylation ~ zebra.gest + 
                                  zebra.peri.concpt + sex + 
                                  age.mon + samp_year_cnt, 
                                data = luma_data_cub)
      
    ## p) Parameter estimates
      summary(cub.gest.zebra.adj)  # model parameter estimates
      confint(cub.gest.zebra.adj)  # 95% CIs 
      
      
  ### 8.8 Cub model: methylation by birth to 3 months prey density     
    ## a) Unadjusted: methlyation by birth to 3 months thomsons density
      cub.birth.3.thomsons.unadj <- glm(methylation ~ thomsons.birth.3, 
                                     data = luma_data_cub)
      
    ## b) Parameter estimates
      summary(cub.birth.3.thomsons.unadj)  # model parameter estimates
      confint(cub.birth.3.thomsons.unadj)  # 95% CIs 
      
    ## c) Adjusted: methlyation by birth to 3 months thomsons density
      cub.birth.3.thomsons.adj <- glm(methylation ~ thomsons.3.6 + 
                                        thomsons.gest + thomsons.peri.concpt + 
                                        sex + age.mon + samp_year_cnt, 
                                   data = luma_data_cub)
      
    ## d) Parameter estimates
      summary(cub.birth.3.thomsons.adj)  # model parameter estimates
      confint(cub.birth.3.thomsons.adj)  # 95% CIs 
      
    ## e) Unadjusted: methlyation by birth to 3 months topi density
      cub.birth.3.topi.unadj <- glm(methylation ~ topi.birth.3, 
                                 data = luma_data_cub)
      
    ## f) Parameter estimates
      summary(cub.birth.3.topi.unadj)  # model parameter estimates
      confint(cub.birth.3.topi.unadj)  # 95% CIs 
      
    ## g) Adjusted: methlyation by birth to 3 months topi density
      cub.birth.3.topi.adj <- glm(methylation ~ topi.birth.3 + topi.gest + 
                                    topi.peri.concpt + sex + 
                                 age.mon + samp_year_cnt, 
                               data = luma_data_cub)
      
    ## h) Parameter estimates
      summary(cub.birth.3.topi.adj)  # model parameter estimates
      confint(cub.birth.3.topi.adj)  # 95% CIs 
      
    ## i) Unadjusted: methlyation by birth to 3 months gnu density
      cub.birth.3.gnu.unadj <- glm(methylation ~ gnu.birth.3, 
                                data = luma_data_cub)
      
    ## j) Parameter estimates
      summary(cub.birth.3.gnu.unadj)  # model parameter estimates
      confint(cub.birth.3.gnu.unadj)  # 95% CIs 
      
    ## k) Adjusted: methlyation by birth to 3 months gnu density
      cub.birth.3.gnu.adj <- glm(methylation ~ gnu.birth.3 + gnu.gest + 
                                   gnu.peri.concpt + sex + 
                                age.mon + samp_year_cnt, 
                              data = luma_data_cub)
      
    ## l) Parameter estimates
      summary(cub.birth.3.gnu.adj)  # model parameter estimates
      confint(cub.birth.3.gnu.adj)  # 95% CIs 
      
    ## m) Unadjusted: methlyation by birth to 3 months zebra density
      cub.birth.3.zebra.unadj <- glm(methylation ~ zebra.birth.3, 
                                  data = luma_data_cub)
      
    ## n) Parameter estimates
      summary(cub.birth.3.zebra.unadj)  # model parameter estimates
      confint(cub.birth.3.zebra.unadj)  # 95% CIs 
      
    ## o) Adjusted: methlyation by birth to 3 months zebra density
      cub.birth.3.zebra.adj <- glm(methylation ~ zebra.birth.3 + zebra.gest + 
                                     zebra.peri.concpt + sex + 
                                  age.mon + samp_year_cnt, 
                                data = luma_data_cub)
      
    ## p) Parameter estimates
      summary(cub.birth.3.zebra.adj)  # model parameter estimates
      confint(cub.birth.3.zebra.adj)  # 95% CIs 
      
      
  ### 8.9 Cub model: methylation by 3 to 6 months prey density     
    ## a) Unadjusted: methlyation by 3 to 6 months thomsons density
      cub.3.6.thomsons.unadj <- glm(methylation ~ thomsons.3.6, 
                                     data = luma_data_cub)
      
    ## b) Parameter estimates
      summary(cub.3.6.thomsons.unadj)  # model parameter estimates
      confint(cub.3.6.thomsons.unadj)  # 95% CIs 
      
    ## c) Adjusted: methlyation by 3 to 6 months thomsons density
      cub.3.6.thomsons.adj <- glm(methylation ~ thomsons.3.6 + 
                                    thomsons.birth.3 + thomsons.gest + 
                                    thomsons.peri.concpt + sex + 
                                     age.mon + samp_year_cnt, 
                                   data = luma_data_cub)
      
    ## d) Parameter estimates
      summary(cub.3.6.thomsons.adj)  # model parameter estimates
      confint(cub.3.6.thomsons.adj)  # 95% CIs 
      
    ## e) Unadjusted: methlyation by 3 to 6 months topi density
      cub.3.6.topi.unadj <- glm(methylation ~ topi.3.6, 
                                 data = luma_data_cub)
      
    ## f) Parameter estimates
      summary(cub.3.6.topi.unadj)  # model parameter estimates
      confint(cub.3.6.topi.unadj)  # 95% CIs 
      
    ## g) Adjusted: methlyation by 3 to 6 months topi density
      cub.3.6.topi.adj <- glm(methylation ~ topi.3.6 + topi.birth.3 + 
                                topi.gest + topi.peri.concpt + sex + 
                                 age.mon + samp_year_cnt, 
                               data = luma_data_cub)
      
    ## h) Parameter estimates
      summary(cub.3.6.topi.adj)  # model parameter estimates
      confint(cub.3.6.topi.adj)  # 95% CIs 
      
    ## i) Unadjusted: methlyation by 3 to 6 months gnu density
      cub.3.6.gnu.unadj <- glm(methylation ~ gnu.3.6, 
                                data = luma_data_cub)
      
    ## j) Parameter estimates
      summary(cub.3.6.gnu.unadj)  # model parameter estimates
      confint(cub.3.6.gnu.unadj)  # 95% CIs 
      
    ## k) Adjusted: methlyation by 3 to 6 months gnu density
      cub.3.6.gnu.adj <- glm(methylation ~ gnu.3.6 + gnu.birth.3 + 
                               gnu.gest + gnu.peri.concpt + sex + 
                                age.mon + samp_year_cnt, 
                              data = luma_data_cub)
      
    ## l) Parameter estimates
      summary(cub.3.6.gnu.adj)  # model parameter estimates
      confint(cub.3.6.gnu.adj)  # 95% CIs 
      
    ## m) Unadjusted: methlyation by 3 to 6 months zebra density
      cub.3.6.zebra.unadj <- glm(methylation ~ zebra.3.6, 
                                  data = luma_data_cub)
      
    ## n) Parameter estimates
      summary(cub.3.6.zebra.unadj)  # model parameter estimates
      confint(cub.3.6.zebra.unadj)  # 95% CIs 
      
    ## o) Adjusted: methlyation by 3 to 6 months zebra density
      cub.3.6.zebra.adj <- glm(methylation ~ zebra.3.6 + zebra.birth.3 + 
                                 zebra.gest + zebra.peri.concpt + sex + 
                                  age.mon + samp_year_cnt, 
                                data = luma_data_cub)
      
    ## p) Parameter estimates
      summary(cub.3.6.zebra.adj)  # model parameter estimates
      confint(cub.3.6.zebra.adj)  # 95% CIs 
      
      
  ### 8.10 Cub model: methylation by 6 to 9 months prey density     
    ## a) Unadjusted: methlyation by 6 to 9 months thomsons density
      cub.6.9.thomsons.unadj <- glm(methylation ~ thomsons.6.9, 
                                     data = luma_data_cub)
      
    ## b) Parameter estimates
      summary(cub.6.9.thomsons.unadj)  # model parameter estimates
      confint(cub.6.9.thomsons.unadj)  # 95% CIs 
      
    ## c) Adjusted: methlyation by 6 to 9 months thomsons density
      cub.6.9.thomsons.adj <- glm(methylation ~ thomsons.6.9 + 
                                    thomsons.3.6 + thomsons.birth.3 + 
                                    thomsons.gest + thomsons.peri.concpt + 
                                    sex + age.mon + samp_year_cnt, 
                                   data = luma_data_cub)
      
    ## d) Parameter estimates
      summary(cub.6.9.thomsons.adj)  # model parameter estimates
      confint(cub.6.9.thomsons.adj)  # 95% CIs 
      
    ## e) Unadjusted: methlyation by 6 to 9 months topi density
      cub.6.9.topi.unadj <- glm(methylation ~ topi.6.9, 
                                 data = luma_data_cub)
      
    ## f) Parameter estimates
      summary(cub.6.9.topi.unadj)  # model parameter estimates
      confint(cub.6.9.topi.unadj)  # 95% CIs 
      
    ## g) Adjusted: methlyation by 6 to 9 months topi density
      cub.6.9.topi.adj <- glm(methylation ~ topi.6.9 + topi.3.6 + 
                                topi.birth.3 + topi.gest + topi.peri.concpt + 
                                sex + age.mon + samp_year_cnt, 
                               data = luma_data_cub)
      
    ## h) Parameter estimates
      summary(cub.6.9.topi.adj)  # model parameter estimates
      confint(cub.6.9.topi.adj)  # 95% CIs 
      
    ## i) Unadjusted: methlyation by 6 to 9 months gnu density
      cub.6.9.gnu.unadj <- glm(methylation ~ gnu.6.9, 
                                data = luma_data_cub)
      
    ## j) Parameter estimates
      summary(cub.6.9.gnu.unadj)  # model parameter estimates
      confint(cub.6.9.gnu.unadj)  # 95% CIs 
      
    ## k) Adjusted: methlyation by 6 to 9 months gnu density
      cub.6.9.gnu.adj <- glm(methylation ~ gnu.6.9 + gnu.3.6 + 
                               gnu.birth.3 + gnu.gest + gnu.peri.concpt + 
                               sex + age.mon + samp_year_cnt, 
                              data = luma_data_cub)
      
    ## l) Parameter estimates
      summary(cub.6.9.gnu.adj)  # model parameter estimates
      confint(cub.6.9.gnu.adj)  # 95% CIs 
      
    ## m) Unadjusted: methlyation by 6 to 9 months zebra density
      cub.6.9.zebra.unadj <- glm(methylation ~ zebra.6.9, 
                                  data = luma_data_cub)
      
    ## n) Parameter estimates
      summary(cub.6.9.zebra.unadj)  # model parameter estimates
      confint(cub.6.9.zebra.unadj)  # 95% CIs 
      
    ## o) Adjusted: methlyation by 6 to 9 months zebra density
      cub.6.9.zebra.adj <- glm(methylation ~ zebra.6.9 + zebra.3.6 + 
                                 zebra.birth.3 + zebra.gest + 
                                 zebra.peri.concpt + sex + 
                                  age.mon + samp_year_cnt, 
                                data = luma_data_cub)
      
    ## p) Parameter estimates
      summary(cub.6.9.zebra.adj)  # model parameter estimates
      confint(cub.6.9.zebra.adj)  # 95% CIs 

      
  ### 8.11 Cub model: mutual adjustment by significatn predictors 
    ## a) Adjusted: methlyation by mom rank and birth to 3 months gnu density
      cub.mutual.adj <- glm(methylation ~ mom.strank.quart + gnu.birth.3 + 
                               gnu.gest + gnu.peri.concpt + sex + 
                                   age.mon + samp_year_cnt, 
                                 data = luma_data_cub)
      
    ## b) Parameter estimates
      summary(cub.mutual.adj)  # model parameter estimates
      confint(cub.mutual.adj)  # 95% CIs 
      
      
      
###############################################################################
##############                 9. Subadult models                ##############
###############################################################################     

  ### 9.1 sub model: methylation by sex
    ## a) Check within strata descritpive stats
      luma_data_sub %>%
        group_by (sex) %>%
        summarise (n.id = n(),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Unadjusted: methlyation by sex
      sub.sex.unadj <- glm(methylation ~  sex, data = luma_data_sub)
      
    ## c) Parameter estimates
      summary(sub.sex.unadj)  # model parameter estimates
      confint(sub.sex.unadj)  # 95% CIs 
      
    ## d) Adjusted: methlyation by sex
      sub.sex.adj <- glm(methylation ~  sex + age.mon +
                           samp_year_cnt,
                         data = luma_data_sub)
      
      ## e) Parameter estimates
      summary(sub.sex.adj)  # model parameter estimates
      confint(sub.sex.adj)  # 95% CIs 
      
      
  ### 9.2 sub model: methylation by age in months
    ## a) Check within strata descritpive stats
      luma_data_sub %>%
        #group_by () %>%
        summarise (n.id = sum(!is.na(age.mon)),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Unadjusted: methlyation by age.mon
      sub.age.mon.unadj <- glm(methylation ~  age.mon, data = luma_data_sub)
      
    ## c) Parameter estimates
      summary(sub.age.mon.unadj)  # model parameter estimates
      confint(sub.age.mon.unadj)  # 95% CIs 
      
    ## d) Adjusted: methlyation by age.mon
      sub.age.mon.adj <- glm(methylation ~  age.mon + sex +
                               samp_year_cnt,
                             data = luma_data_sub)
      
    ## e) Parameter estimates
      summary(sub.age.mon.adj)  # model parameter estimates
      confint(sub.age.mon.adj)  # 95% CIs 
      
      
  ### 9.3 sub model: methylation by mom rank quartiles    
    ## a) Check within strata descritpive stats
      luma_data_sub %>%
        group_by (mom.strank.quart) %>%
        summarise (n.id = sum(!is.na(age.mon)),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Unadjusted: methlyation by mom.strank.quart
      sub.mom.rank.unadj <- glm(methylation ~  mom.strank.quart, 
                                data = luma_data_sub)
      
    ## c) Parameter estimates
      summary(sub.mom.rank.unadj)  # model parameter estimates
      confint(sub.mom.rank.unadj)  # 95% CIs 
      # use 'car' package to extract a Wald Chi-Sq p-value; test for overall 
      # difference
      Anova(sub.mom.rank.unadj, Type ="II", test = "Wald") # Wald test p
      
    ## d) Adjusted: methlyation by mom.strank.quart
      sub.mom.rank.adj <- glm(methylation ~ mom.strank.quart + sex + 
                                    age.mon + samp_year_cnt,
                                  data = luma_data_sub)
      
    ## e) Parameter estimates
      summary(sub.mom.rank.adj)  # print model summary, effects and SE
      confint(sub.mom.rank.adj)  # print 95% CIs for parameter estimates
      Anova(sub.mom.rank.adj, Type ="II", test = "Wald") # Wald test p

    ## f) Extract mom.strank.quart estimates and 
      sub.rank.ef <- effect("mom.strank.quart", sub.mom.rank.adj)
      summary(sub.rank.ef)
      # Save effects as data frame
      sub.rank.ef.table <- as.data.frame(sub.rank.ef)
      # Set the reference level to Q1
      sub.rank.ef.table <- transform( sub.rank.ef.table,
                                      mom.strank.quart = 
                                        factor(mom.strank.quart,
                                               levels = c("Q1 (lowest)",
                                                          "Q2", "Q3",
                                                          "Q4 (highest)")))
      
    ## g) Graph sub.mom.rank effects
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
      
    ## h) Save Plot
      # use ggsave to save the linearization plot
      ggsave("mat_rank_sub_mod_beta.pdf", plot = last_plot(), device = NULL,
             path = "./output/output_luma_ecolog", scale = 1, width = 7,
             height = 5.5,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## i) Do a post-hoc test to determine if Q2, Q3, and Q4 differ
      pairwise.t.test(luma_data_sub$methylation,
                      luma_data_sub$mom.strank.quart,
                      p.adj = "none")
      TukeyHSD(aov(methylation ~ mom.strank.quart + sex,
                   data = luma_data_sub))
      
    ## j) Combine quartiles 2-4 into a single category
      luma_data_sub$mom.strank.quart.comb <- as.factor(
        ifelse(luma_data_sub$mom.strank.quart ==  "Q1 (lowest)",
               "Q1 (lowest)", "Q2-Q4 (highest)"))
      
      ## k)  Adjusted: methlyation by mom.strank.quart binned
      sub.rank.adj2 <- glm(methylation ~ mom.strank.quart.comb + sex + 
                             age.mon + samp_year_cnt,
                           data = luma_data_sub)
      
      summary(sub.rank.adj2)  # print model summary, effects and SE
      confint(sub.rank.mod2)  # print 95% CIs for parameter estimates
      Anova(sub.mom.rank.adj2, Type ="II", test = "Wald") # Wald test p
      
      
  ### 9.4 sub model: methylation by litter size
    ## a) Check within strata descritpive stats
      luma_data_sub %>%
        group_by (lit.size) %>%

        summarise (n.id = sum(!is.na(lit.size)),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Unadjusted: methlyation by age.mon
      sub.lit.size.unadj <- glm(methylation ~  lit.size, data = luma_data_sub)
      
    ## c) Parameter estimates
      summary(sub.lit.size.unadj)  # model parameter estimates
      confint(sub.lit.size.unadj)  # 95% CIs
      Anova(sub.lit.size.unadj, Type ="II", test = "Wald") # Wald test p
      
    ## d) Adjusted: methlyation by age.mon
      sub.lit.size.adj <- glm(methylation ~  lit.size + age.mon + sex +
                                samp_year_cnt,
                              data = luma_data_sub)
      
    ## e) Parameter estimates
      summary(sub.lit.size.adj)  # model parameter estimates
      confint(sub.lit.size.adj)  # 95% CIs 
      Anova(sub.lit.size.adj, Type ="II", test = "Wald") # Wald test p
    
      
  ### 9.5 sub model: methylation by human presence proxy
    ## a) Check within strata descritpive stats
      luma_data_sub %>%
        group_by (hum.pres) %>%
        summarise (n.id = sum(!is.na(lit.size)),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Unadjusted: methlyation by age.mon
      sub.hum.pres.unadj <- glm(methylation ~  hum.pres, data = luma_data_sub)
      
    ## c) Parameter estimates
      summary(sub.hum.pres.unadj)  # model parameter estimates
      confint(sub.hum.pres.unadj)  # 95% CIs
      Anova(sub.hum.pres.unadj, Type ="II", test = "Wald") # Wald test p
      
    ## d) Adjusted: methlyation by age.mon
      sub.hum.pres.adj <- glm(methylation ~  hum.pres + age.mon + sex +
                                samp_year_cnt,
                              data = luma_data_sub)
      
    ## e) Parameter estimates
      summary(sub.hum.pres.adj)  # model parameter estimates
      confint(sub.hum.pres.adj)  # 95% CIs 
      Anova(sub.hum.pres.adj, Type ="II", test = "Wald") # Wald test p
    
      
  ### 9.6 sub model: methylation by periconceptional prey density     
    ## a) Unadjusted: methlyation by periconceptional thomsons density
      sub.peri.thomsons.unadj <- glm(methylation ~ thomsons.peri.concpt, 
                                     data = luma_data_sub)
      
    ## b) Parameter estimates
      summary(sub.peri.thomsons.unadj)  # model parameter estimates
      confint(sub.peri.thomsons.unadj)  # 95% CIs 
      
    ## c) Adjusted: methlyation by periconceptional thomsons density
      sub.peri.thomsons.adj <- glm(methylation ~ thomsons.peri.concpt + sex + 
                                     age.mon + samp_year_cnt, 
                                   data = luma_data_sub)
      
    ## d) Parameter estimates
      summary(sub.peri.thomsons.adj)  # model parameter estimates
      confint(sub.peri.thomsons.adj)  # 95% CIs 
      
    ## e) Unadjusted: methlyation by periconceptional topi density
      sub.peri.topi.unadj <- glm(methylation ~ topi.peri.concpt, 
                                 data = luma_data_sub)
      
    ## f) Parameter estimates
      summary(sub.peri.topi.unadj)  # model parameter estimates
      confint(sub.peri.topi.unadj)  # 95% CIs 
      
    ## g) Adjusted: methlyation by periconceptional topi density
      sub.peri.topi.adj <- glm(methylation ~ topi.peri.concpt + sex + 
                                 age.mon + samp_year_cnt, 
                               data = luma_data_sub)
      
    ## h) Parameter estimates
      summary(sub.peri.topi.adj)  # model parameter estimates
      confint(sub.peri.topi.adj)  # 95% CIs 
      
    ## i) Unadjusted: methlyation by periconceptional gnu density
      sub.peri.gnu.unadj <- glm(methylation ~ gnu.peri.concpt, 
                                data = luma_data_sub)
      
    ## j) Parameter estimates
      summary(sub.peri.gnu.unadj)  # model parameter estimates
      confint(sub.peri.gnu.unadj)  # 95% CIs 
      
    ## k) Adjusted: methlyation by periconceptional gnu density
      sub.peri.gnu.adj <- glm(methylation ~ gnu.peri.concpt + sex + 
                                age.mon + samp_year_cnt, 
                              data = luma_data_sub)
      
    ## l) Parameter estimates
      summary(sub.peri.gnu.adj)  # model parameter estimates
      confint(sub.peri.gnu.adj)  # 95% CIs 
      
    ## m) Unadjusted: methlyation by periconceptional zebra density
      sub.peri.zebra.unadj <- glm(methylation ~ zebra.peri.concpt, 
                                  data = luma_data_sub)
      
    ## n) Parameter estimates
      summary(sub.peri.zebra.unadj)  # model parameter estimates
      confint(sub.peri.zebra.unadj)  # 95% CIs 
      
    ## o) Adjusted: methlyation by periconceptional zebra density
      sub.peri.zebra.adj <- glm(methylation ~ zebra.peri.concpt + sex + 
                                  age.mon + samp_year_cnt, 
                                data = luma_data_sub)
      
    ## p) Parameter estimates
      summary(sub.peri.zebra.adj)  # model parameter estimates
      confint(sub.peri.zebra.adj)  # 95% CIs 
      
      
      
  ### 9.7 sub model: methylation by gestational prey density     
    ## a) Unadjusted: methlyation by gestational thomsons density
      sub.gest.thomsons.unadj <- glm(methylation ~ thomsons.gest, 
                                     data = luma_data_sub)
      
    ## b) Parameter estimates
      summary(sub.gest.thomsons.unadj)  # model parameter estimates
      confint(sub.gest.thomsons.unadj)  # 95% CIs 
      
    ## c) Adjusted: methlyation by gestational thomsons density
      sub.gest.thomsons.adj <- glm(methylation ~ thomsons.gest + 
                                     thomsons.peri.concpt + sex + 
                                     age.mon + samp_year_cnt, 
                                   data = luma_data_sub)
      
    ## d) Parameter estimates
      summary(sub.gest.thomsons.adj)  # model parameter estimates
      confint(sub.gest.thomsons.adj)  # 95% CIs 
      
    ## e) Unadjusted: methlyation by gestational topi density
      sub.gest.topi.unadj <- glm(methylation ~ topi.peri.concpt, 
                                 data = luma_data_sub)
      
    ## f) Parameter estimates
      summary(sub.gest.topi.unadj)  # model parameter estimates
      confint(sub.gest.topi.unadj)  # 95% CIs 
      
    ## g) Adjusted: methlyation by gestational topi density
      sub.gest.topi.adj <- glm(methylation ~ topi.gest + topi.peri.concpt + sex + 
                                 age.mon + samp_year_cnt, 
                               data = luma_data_sub)
    ## h) Parameter estimates
      summary(sub.gest.topi.adj)  # model parameter estimates
      confint(sub.gest.topi.adj)  # 95% CIs 
      
    ## i) Unadjusted: methlyation by gestational gnu density
      sub.gest.gnu.unadj <- glm(methylation ~ gnu.gest, 
                                data = luma_data_sub)
      
    ## j) Parameter estimates
      summary(sub.gest.gnu.unadj)  # model parameter estimates
      confint(sub.gest.gnu.unadj)  # 95% CIs 
      

    ## k) Adjusted: methlyation by gestational gnu density
      sub.gest.gnu.adj <- glm(methylation ~ gnu.gest + gnu.peri.concpt + sex + 
                                age.mon + samp_year_cnt, 
                              data = luma_data_sub)
      
    ## l) Parameter estimates
      summary(sub.gest.gnu.adj)  # model parameter estimates
      confint(sub.gest.gnu.adj)  # 95% CIs 
      
    ## m) Unadjusted: methlyation by gestational zebra density
      sub.gest.zebra.unadj <- glm(methylation ~ zebra.gest,
                                  data = luma_data_sub)
      
    ## n) Parameter estimates
      summary(sub.gest.zebra.unadj)  # model parameter estimates
      confint(sub.gest.zebra.unadj)  # 95% CIs 
      
    ## o) Adjusted: methlyation by gestational zebra density
      sub.gest.zebra.adj <- glm(methylation ~ zebra.gest + 
                                  zebra.peri.concpt + sex + 
                                  age.mon + samp_year_cnt, 
                                data = luma_data_sub)
      
    ## p) Parameter estimates
      summary(sub.gest.zebra.adj)  # model parameter estimates
      confint(sub.gest.zebra.adj)  # 95% CIs 
      
      
  ### 9.8 sub model: methylation by birth to 3 months prey density     
    ## a) Unadjusted: methlyation by birth to 3 months thomsons density
      sub.birth.3.thomsons.unadj <- glm(methylation ~ thomsons.birth.3, 
                                        data = luma_data_sub)
      
    ## b) Parameter estimates
      summary(sub.birth.3.thomsons.unadj)  # model parameter estimates
      confint(sub.birth.3.thomsons.unadj)  # 95% CIs 
      
    ## c) Adjusted: methlyation by birth to 3 months thomsons density
      sub.birth.3.thomsons.adj <- glm(methylation ~ thomsons.3.6 + 
                                        thomsons.gest + thomsons.peri.concpt + 
                                        sex + age.mon + samp_year_cnt, 
                                      data = luma_data_sub)
      
    ## d) Parameter estimates
      summary(sub.birth.3.thomsons.adj)  # model parameter estimates
      confint(sub.birth.3.thomsons.adj)  # 95% CIs 
      
    ## e) Unadjusted: methlyation by birth to 3 months topi density
      sub.birth.3.topi.unadj <- glm(methylation ~ topi.birth.3, 
                                    data = luma_data_sub)
      
    ## f) Parameter estimates
      summary(sub.birth.3.topi.unadj)  # model parameter estimates
      confint(sub.birth.3.topi.unadj)  # 95% CIs 
      
    ## g) Adjusted: methlyation by birth to 3 months topi density
      sub.birth.3.topi.adj <- glm(methylation ~ topi.birth.3 + topi.gest + 
                                    topi.peri.concpt + sex + 
                                    age.mon + samp_year_cnt, 
                                  data = luma_data_sub)
      
    ## h) Parameter estimates
      summary(sub.birth.3.topi.adj)  # model parameter estimates
      confint(sub.birth.3.topi.adj)  # 95% CIs 
      
    ## i) Unadjusted: methlyation by birth to 3 months gnu density
      sub.birth.3.gnu.unadj <- glm(methylation ~ gnu.birth.3, 
                                   data = luma_data_sub)
      
    ## j) Parameter estimates
      summary(sub.birth.3.gnu.unadj)  # model parameter estimates
      confint(sub.birth.3.gnu.unadj)  # 95% CIs 
      
    ## k) Adjusted: methlyation by birth to 3 months gnu density
      sub.birth.3.gnu.adj <- glm(methylation ~ gnu.birth.3 + gnu.gest + 
                                   gnu.peri.concpt + sex + 
                                   age.mon + samp_year_cnt, 
                                 data = luma_data_sub)
      
    ## l) Parameter estimates
      summary(sub.birth.3.gnu.adj)  # model parameter estimates
      confint(sub.birth.3.gnu.adj)  # 95% CIs 
      
    ## m) Unadjusted: methlyation by birth to 3 months zebra density
      sub.birth.3.zebra.unadj <- glm(methylation ~ zebra.birth.3, 
                                     data = luma_data_sub)
      
    ## n) Parameter estimates
      summary(sub.birth.3.zebra.unadj)  # model parameter estimates
      confint(sub.birth.3.zebra.unadj)  # 95% CIs 
      
    ## o) Adjusted: methlyation by birth to 3 months zebra density
      sub.birth.3.zebra.adj <- glm(methylation ~ zebra.birth.3 + zebra.gest + 
                                     zebra.peri.concpt + sex + 
                                     age.mon + samp_year_cnt, 
                                   data = luma_data_sub)
      
    ## p) Parameter estimates
      summary(sub.birth.3.zebra.adj)  # model parameter estimates
      confint(sub.birth.3.zebra.adj)  # 95% CIs 
      
      
      
  ### 8.9 sub model: methylation by 3 to 6 months prey density     
    ## a) Unadjusted: methlyation by 3 to 6 months thomsons density
      sub.3.6.thomsons.unadj <- glm(methylation ~ thomsons.3.6, 
                                    data = luma_data_sub)
      
    ## b) Parameter estimates
      summary(sub.3.6.thomsons.unadj)  # model parameter estimates
      confint(sub.3.6.thomsons.unadj)  # 95% CIs 
      
    ## c) Adjusted: methlyation by 3 to 6 months thomsons density
      sub.3.6.thomsons.adj <- glm(methylation ~ thomsons.3.6 + 
                                    thomsons.birth.3 + thomsons.gest + 
                                    thomsons.peri.concpt + sex + 
                                    age.mon + samp_year_cnt, 
                                  data = luma_data_sub)
      
    ## d) Parameter estimates
      summary(sub.3.6.thomsons.adj)  # model parameter estimates
      confint(sub.3.6.thomsons.adj)  # 95% CIs 
      
    ## e) Unadjusted: methlyation by 3 to 6 months topi density
      sub.3.6.topi.unadj <- glm(methylation ~ topi.3.6, 
                                data = luma_data_sub)
      
    ## f) Parameter estimates
      summary(sub.3.6.topi.unadj)  # model parameter estimates
      confint(sub.3.6.topi.unadj)  # 95% CIs 
      
    ## g) Adjusted: methlyation by 3 to 6 months topi density
      sub.3.6.topi.adj <- glm(methylation ~ topi.3.6 + topi.birth.3 + 
                                topi.gest + topi.peri.concpt + sex + 
                                age.mon + samp_year_cnt, 
                              data = luma_data_sub)
      
    ## h) Parameter estimates
      summary(sub.3.6.topi.adj)  # model parameter estimates
      confint(sub.3.6.topi.adj)  # 95% CIs 
      
    ## i) Unadjusted: methlyation by 3 to 6 months gnu density
      sub.3.6.gnu.unadj <- glm(methylation ~ gnu.3.6, 
                               data = luma_data_sub)
      
    ## j) Parameter estimates
      summary(sub.3.6.gnu.unadj)  # model parameter estimates
      confint(sub.3.6.gnu.unadj)  # 95% CIs 
      
    ## k) Adjusted: methlyation by 3 to 6 months gnu density
      sub.3.6.gnu.adj <- glm(methylation ~ gnu.3.6 + gnu.birth.3 + 
                               gnu.gest + gnu.peri.concpt + sex + 
                               age.mon + samp_year_cnt, 
                             data = luma_data_sub)
      
    ## l) Parameter estimates
      summary(sub.3.6.gnu.adj)  # model parameter estimates
      confint(sub.3.6.gnu.adj)  # 95% CIs 
      
    ## m) Unadjusted: methlyation by 3 to 6 months zebra density
      sub.3.6.zebra.unadj <- glm(methylation ~ zebra.3.6, 
                                 data = luma_data_sub)
      
    ## n) Parameter estimates
      summary(sub.3.6.zebra.unadj)  # model parameter estimates
      confint(sub.3.6.zebra.unadj)  # 95% CIs 
      
    ## o) Adjusted: methlyation by 3 to 6 months zebra density
      sub.3.6.zebra.adj <- glm(methylation ~ zebra.3.6 + zebra.birth.3 + 
                                 zebra.gest + zebra.peri.concpt + sex + 
                                 age.mon + samp_year_cnt, 
                               data = luma_data_sub)
      
    ## p) Parameter estimates
      summary(sub.3.6.zebra.adj)  # model parameter estimates
      confint(sub.3.6.zebra.adj)  # 95% CIs 
      
      
      
  ### 9.10 sub model: methylation by 6 to 9 months prey density     
    ## a) Unadjusted: methlyation by 6 to 9 months thomsons density
      sub.6.9.thomsons.unadj <- glm(methylation ~ thomsons.6.9, 
                                    data = luma_data_sub)
      
    ## b) Parameter estimates
      summary(sub.6.9.thomsons.unadj)  # model parameter estimates
      confint(sub.6.9.thomsons.unadj)  # 95% CIs 
      
    ## c) Adjusted: methlyation by 6 to 9 months thomsons density
      sub.6.9.thomsons.adj <- glm(methylation ~ thomsons.6.9 + 
                                    thomsons.3.6 + thomsons.birth.3 + 
                                    thomsons.gest + thomsons.peri.concpt + 
                                    sex + age.mon + samp_year_cnt, 
                                  data = luma_data_sub)
      
    ## d) Parameter estimates
      summary(sub.6.9.thomsons.adj)  # model parameter estimates
      confint(sub.6.9.thomsons.adj)  # 95% CIs 
      
    ## e) Unadjusted: methlyation by 6 to 9 months topi density
      sub.6.9.topi.unadj <- glm(methylation ~ topi.6.9, 
                                data = luma_data_sub)
      
    ## f) Parameter estimates
      summary(sub.6.9.topi.unadj)  # model parameter estimates
      confint(sub.6.9.topi.unadj)  # 95% CIs 
      
    ## g) Adjusted: methlyation by 6 to 9 months topi density
      sub.6.9.topi.adj <- glm(methylation ~ topi.6.9 + topi.3.6 + 
                                topi.birth.3 + topi.gest + topi.peri.concpt + 
                                sex + age.mon + samp_year_cnt, 
                              data = luma_data_sub)
      
    ## h) Parameter estimates
      summary(sub.6.9.topi.adj)  # model parameter estimates
      confint(sub.6.9.topi.adj)  # 95% CIs 
      
    ## i) Unadjusted: methlyation by 6 to 9 months gnu density
      sub.6.9.gnu.unadj <- glm(methylation ~ gnu.6.9, 
                               data = luma_data_sub)
      
    ## j) Parameter estimates
      summary(sub.6.9.gnu.unadj)  # model parameter estimates
      confint(sub.6.9.gnu.unadj)  # 95% CIs 
      
    ## k) Adjusted: methlyation by 6 to 9 months gnu density
      sub.6.9.gnu.adj <- glm(methylation ~ gnu.6.9 + gnu.3.6 + 
                               gnu.birth.3 + gnu.gest + gnu.peri.concpt + 
                               sex + age.mon + samp_year_cnt, 
                             data = luma_data_sub)
      
    ## l) Parameter estimates
      summary(sub.6.9.gnu.adj)  # model parameter estimates
      confint(sub.6.9.gnu.adj)  # 95% CIs 
      
    ## m) Unadjusted: methlyation by 6 to 9 months zebra density
      sub.6.9.zebra.unadj <- glm(methylation ~ zebra.6.9, 
                                 data = luma_data_sub)
      
    ## n) Parameter estimates
      summary(sub.6.9.zebra.unadj)  # model parameter estimates
      confint(sub.6.9.zebra.unadj)  # 95% CIs 
      
    ## o) Adjusted: methlyation by 6 to 9 months zebra density
      sub.6.9.zebra.adj <- glm(methylation ~ zebra.6.9 + zebra.3.6 + 
                                 zebra.birth.3 + zebra.gest + 
                                 zebra.peri.concpt + sex + 
                                 age.mon + samp_year_cnt, 
                               data = luma_data_sub)
      
    ## p) Parameter estimates
      summary(sub.6.9.zebra.adj)  # model parameter estimates
      confint(sub.6.9.zebra.adj)  # 95% CIs 


      
###############################################################################
##############                  10. Adult models                  ##############
############################################################################### 
      
  ### 10.1 adult model: methylation by sex
    ## a) Check within strata descritpive stats
      luma_data_adult %>%
        group_by (sex) %>%
        summarise (n.id = n(),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))

    ## b) Unadjusted: methlyation by sex
      adult.sex.unadj <- glm(methylation ~  sex, data = luma_data_adult)
      
    ## c) Parameter estimates
      summary(adult.sex.unadj)  # model parameter estimates
      confint(adult.sex.unadj)  # 95% CIs 
      
    ## d) Adjusted: methlyation by sex
      adult.sex.adj <- glm(methylation ~  sex + age.mon +
                           samp_year_cnt,
                         data = luma_data_adult)
      
    ## e) Parameter estimates
      summary(adult.sex.adj)  # model parameter estimates
      confint(adult.sex.adj)  # 95% CIs 
      
      
  ### 10.2 adult model: methylation by age in months
    ## a) Check within strata descritpive stats
      luma_data_adult %>%
        #group_by () %>%
        summarise (n.id = sum(!is.na(age.mon)),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Unadjusted: methlyation by age.mon
      adult.age.mon.unadj <- glm(methylation ~  age.mon, 
                                 data = luma_data_adult)

    ## c) Parameter estimates
      summary(adult.age.mon.unadj)  # model parameter estimates
      confint(adult.age.mon.unadj)  # 95% CIs 
      
    ## d) Adjusted: methlyation by age.mon
      adult.age.mon.adj <- glm(methylation ~  age.mon + sex +
                               samp_year_cnt,
                             data = luma_data_adult)
      
    ## e) Parameter estimates
      summary(adult.age.mon.adj)  # model parameter estimates
      confint(adult.age.mon.adj)  # 95% CIs 
      
      
  ### 10.3 adult model: methylation by mom rank quartiles    
    ## a) Check within strata descritpive stats
      luma_data_adult %>%
        group_by (mom.strank.quart) %>%
        summarise (n.id = sum(!is.na(age.mon)),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Unadjusted: methlyation by mom.strank.quart
      adult.mom.rank.unadj <- glm(methylation ~  mom.strank.quart, 
                                data = luma_data_adult)
      
    ## c) Parameter estimates
      summary(adult.mom.rank.unadj)  # model parameter estimates
      confint(adult.mom.rank.unadj)  # 95% CIs 
      # use 'car' package to extract a Wald Chi-Sq p-value; test for overall 
      # difference
      Anova(adult.mom.rank.unadj, Type ="II", test = "Wald") # Wald test p
      
    ## d) Adjusted: methlyation by mom.strank.quart
      adult.mom.rank.adj <- glm(methylation ~ mom.strank.quart + sex + 
                                age.mon + samp_year_cnt,
                              data = luma_data_adult)
      
    ## e) Parameter estimates
      summary(adult.mom.rank.adj)  # print model summary, effects and SE
      confint(adult.mom.rank.adj)  # print 95% CIs for parameter estimates
      Anova(adult.mom.rank.adj, Type ="II", test = "Wald") # Wald test p
      
    ## f) Extract mom.strank.quart estimates and 
      adult.rank.ef <- effect("mom.strank.quart", adult.mom.rank.adj)
      summary(adult.rank.ef)
      # Save effects as data frame
      adult.rank.ef.table <- as.data.frame(adult.rank.ef)
      # Set the reference level to Q1
      adult.rank.ef.table <- transform( adult.rank.ef.table,
                                      mom.strank.quart = 
                                        factor(mom.strank.quart,
                                               levels = c("Q1 (lowest)",
                                                          "Q2", "Q3",
                                                          "Q4 (highest)")))
      
    ## g) Graph adult.mom.rank effects
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
      
    ## h) Save Plot
      # use ggsave to save the linearization plot
      ggsave("mat_rank_adult_mod_beta.pdf", plot = last_plot(), device = NULL,
             path = "./output/output_luma_ecolog", scale = 1, width = 7,
             height = 5.5,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## i) Do a post-hoc test to determine if Q2, Q3, and Q4 differ
      pairwise.t.test(luma_data_adult$methylation,
                      luma_data_adult$mom.strank.quart,
                      p.adj = "none")
      TukeyHSD(aov(methylation ~ mom.strank.quart + sex,
                   data = luma_data_adult))
      
    ## j) Combine quartiles 2-4 into a single category
      luma_data_adult$mom.strank.quart.comb <- as.factor(
        ifelse(luma_data_adult$mom.strank.quart ==  "Q1 (lowest)",
               "Q1 (lowest)", "Q2-Q4 (highest)"))
      
    ## k)  Adjusted: methlyation by mom.strank.quart binned
      adult.rank.adj2 <- glm(methylation ~ mom.strank.quart.comb + sex + 
                             age.mon + samp_year_cnt,
                           data = luma_data_adult)
      
      summary(adult.rank.adj2)  # print model summary, effects and SE
      confint(adult.rank.mod2)  # print 95% CIs for parameter estimates
      Anova(adult.mom.rank.adj2, Type ="II", test = "Wald") # Wald test p
      
      
  ### 10.4 adult model: methylation by litter size
    ## a) Check within strata descritpive stats
      luma_data_adult %>%
        group_by (lit.size) %>%
        summarise (n.id = sum(!is.na(lit.size)),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Unadjusted: methlyation by age.mon
      adult.lit.size.unadj <- glm(methylation ~  lit.size, 
                                  data = luma_data_adult)
      
    ## c) Parameter estimates
      summary(adult.lit.size.unadj)  # model parameter estimates
      confint(adult.lit.size.unadj)  # 95% CIs
      Anova(adult.lit.size.unadj, Type ="II", test = "Wald") # Wald test p

    ## d) Adjusted: methlyation by age.mon
      adult.lit.size.adj <- glm(methylation ~  lit.size + age.mon + sex +
                                samp_year_cnt,
                              data = luma_data_adult)
      
    ## e) Parameter estimates
      summary(adult.lit.size.adj)  # model parameter estimates
      confint(adult.lit.size.adj)  # 95% CIs 
      Anova(adult.lit.size.adj, Type ="II", test = "Wald") # Wald test p
      
      
  ### 10.5 adult model: methylation by human presence proxy
    ## a) Check within strata descritpive stats
      luma_data_adult %>%
        group_by (hum.pres) %>%
        summarise (n.id = sum(!is.na(lit.size)),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Unadjusted: methlyation by age.mon
      adult.hum.pres.unadj <- glm(methylation ~  hum.pres, 
                                  data = luma_data_adult)
      
    ## c) Parameter estimates
      summary(adult.hum.pres.unadj)  # model parameter estimates
      confint(adult.hum.pres.unadj)  # 95% CIs
      Anova(adult.hum.pres.unadj, Type ="II", test = "Wald") # Wald test p
      
    ## d) Adjusted: methlyation by age.mon
      adult.hum.pres.adj <- glm(methylation ~  hum.pres + age.mon + sex +
                                samp_year_cnt,
                              data = luma_data_adult)
      
    ## e) Parameter estimates
      summary(adult.hum.pres.adj)  # model parameter estimates
      confint(adult.hum.pres.adj)  # 95% CIs 
      Anova(adult.hum.pres.adj, Type ="II", test = "Wald") # Wald test p
      
      
  ### 10.6 adult model: methylation by periconceptional prey density     
    ## a) Unadjusted: methlyation by periconceptional thomsons density
      adult.peri.thomsons.unadj <- glm(methylation ~ thomsons.peri.concpt, 
                                     data = luma_data_adult)
      
    ## b) Parameter estimates
      summary(adult.peri.thomsons.unadj)  # model parameter estimates
      confint(adult.peri.thomsons.unadj)  # 95% CIs 
      

    ## c) Adjusted: methlyation by periconceptional thomsons density
      adult.peri.thomsons.adj <- glm(methylation ~ thomsons.peri.concpt + sex + 
                                     age.mon + samp_year_cnt, 
                                   data = luma_data_adult)

    ## d) Parameter estimates
      summary(adult.peri.thomsons.adj)  # model parameter estimates
      confint(adult.peri.thomsons.adj)  # 95% CIs 
      
    ## e) Unadjusted: methlyation by periconceptional topi density
      adult.peri.topi.unadj <- glm(methylation ~ topi.peri.concpt, 
                                 data = luma_data_adult)
      
    ## f) Parameter estimates
      summary(adult.peri.topi.unadj)  # model parameter estimates
      confint(adult.peri.topi.unadj)  # 95% CIs 
      
    ## g) Adjusted: methlyation by periconceptional topi density
      adult.peri.topi.adj <- glm(methylation ~ topi.peri.concpt + sex + 
                                 age.mon + samp_year_cnt, 
                               data = luma_data_adult)
    ## h) Parameter estimates
      summary(adult.peri.topi.adj)  # model parameter estimates
      confint(adult.peri.topi.adj)  # 95% CIs 
      
    ## i) Unadjusted: methlyation by periconceptional gnu density
      adult.peri.gnu.unadj <- glm(methylation ~ gnu.peri.concpt, 
                                data = luma_data_adult)
      
    ## j) Parameter estimates
      summary(adult.peri.gnu.unadj)  # model parameter estimates
      confint(adult.peri.gnu.unadj)  # 95% CIs 
      
    ## k) Adjusted: methlyation by periconceptional gnu density
      adult.peri.gnu.adj <- glm(methylation ~ gnu.peri.concpt + sex + 
                                age.mon + samp_year_cnt, 
                              data = luma_data_adult)
      
    ## l) Parameter estimates
      summary(adult.peri.gnu.adj)  # model parameter estimates
      confint(adult.peri.gnu.adj)  # 95% CIs 
      
    ## m) Unadjusted: methlyation by periconceptional zebra density
      adult.peri.zebra.unadj <- glm(methylation ~ zebra.peri.concpt, 
                                  data = luma_data_adult)
      
    ## n) Parameter estimates
      summary(adult.peri.zebra.unadj)  # model parameter estimates
      confint(adult.peri.zebra.unadj)  # 95% CIs 
      
    ## o) Adjusted: methlyation by periconceptional zebra density
      adult.peri.zebra.adj <- glm(methylation ~ zebra.peri.concpt + sex + 
                                  age.mon + samp_year_cnt, 
                                data = luma_data_adult)
      
    ## p) Parameter estimates
      summary(adult.peri.zebra.adj)  # model parameter estimates
      confint(adult.peri.zebra.adj)  # 95% CIs 
      
      
  ### 10.7 adult model: methylation by gestational prey density     
    ## a) Unadjusted: methlyation by gestational thomsons density
      adult.gest.thomsons.unadj <- glm(methylation ~ thomsons.gest, 
                                     data = luma_data_adult)
      
    ## b) Parameter estimates
      summary(adult.gest.thomsons.unadj)  # model parameter estimates
      confint(adult.gest.thomsons.unadj)  # 95% CIs 
      
    ## c) Adjusted: methlyation by gestational thomsons density
      adult.gest.thomsons.adj <- glm(methylation ~ thomsons.gest + 
                                     thomsons.peri.concpt + sex + 
                                     age.mon + samp_year_cnt, 
                                   data = luma_data_adult)
      
    ## d) Parameter estimates
      summary(adult.gest.thomsons.adj)  # model parameter estimates
      confint(adult.gest.thomsons.adj)  # 95% CIs 
      
    ## e) Unadjusted: methlyation by gestational topi density
      adult.gest.topi.unadj <- glm(methylation ~ topi.peri.concpt, 
                                 data = luma_data_adult)
      
    ## f) Parameter estimates
      summary(adult.gest.topi.unadj)  # model parameter estimates
      confint(adult.gest.topi.unadj)  # 95% CIs 
      
    ## g) Adjusted: methlyation by gestational topi density
      adult.gest.topi.adj <- glm(methylation ~ topi.gest + topi.peri.concpt + 
                                   sex + age.mon + samp_year_cnt, 
                               data = luma_data_adult)
    ## h) Parameter estimates
      summary(adult.gest.topi.adj)  # model parameter estimates
      confint(adult.gest.topi.adj)  # 95% CIs 
      
    ## i) Unadjusted: methlyation by gestational gnu density
      adult.gest.gnu.unadj <- glm(methylation ~ gnu.gest, 
                                data = luma_data_adult)
      
    ## j) Parameter estimates
      summary(adult.gest.gnu.unadj)  # model parameter estimates
      confint(adult.gest.gnu.unadj)  # 95% CIs 
      
    ## k) Adjusted: methlyation by gestational gnu density
      adult.gest.gnu.adj <- glm(methylation ~ gnu.gest + gnu.peri.concpt + 
                                  sex + age.mon + samp_year_cnt, 
                              data = luma_data_adult)
      
    ## l) Parameter estimates
      summary(adult.gest.gnu.adj)  # model parameter estimates
      confint(adult.gest.gnu.adj)  # 95% CIs 
      
    ## m) Unadjusted: methlyation by gestational zebra density
      adult.gest.zebra.unadj <- glm(methylation ~ zebra.gest,
                                  data = luma_data_adult)
      
    ## n) Parameter estimates
      summary(adult.gest.zebra.unadj)  # model parameter estimates
      confint(adult.gest.zebra.unadj)  # 95% CIs 
      
    ## o) Adjusted: methlyation by gestational zebra density
      adult.gest.zebra.adj <- glm(methylation ~ zebra.gest + 
                                  zebra.peri.concpt + sex + 
                                  age.mon + samp_year_cnt, 
                                data = luma_data_adult)
      
    ## p) Parameter estimates
      summary(adult.gest.zebra.adj)  # model parameter estimates
      confint(adult.gest.zebra.adj)  # 95% CIs 
      
      
  ### 10.8 adult model: methylation by birth to 3 months prey density     
    ## a) Unadjusted: methlyation by birth to 3 months thomsons density
      adult.birth.3.thomsons.unadj <- glm(methylation ~ thomsons.birth.3, 
                                        data = luma_data_adult)
      
    ## b) Parameter estimates
      summary(adult.birth.3.thomsons.unadj)  # model parameter estimates
      confint(adult.birth.3.thomsons.unadj)  # 95% CIs 
      
    ## c) Adjusted: methlyation by birth to 3 months thomsons density
      adult.birth.3.thomsons.adj <- glm(methylation ~ thomsons.3.6 + 
                                        thomsons.gest + thomsons.peri.concpt + 
                                        sex + age.mon + samp_year_cnt, 
                                      data = luma_data_adult)
      
    ## d) Parameter estimates
      summary(adult.birth.3.thomsons.adj)  # model parameter estimates
      confint(adult.birth.3.thomsons.adj)  # 95% CIs 
      
    ## e) Unadjusted: methlyation by birth to 3 months topi density
      adult.birth.3.topi.unadj <- glm(methylation ~ topi.birth.3, 
                                    data = luma_data_adult)
      
    ## f) Parameter estimates
      summary(adult.birth.3.topi.unadj)  # model parameter estimates
      confint(adult.birth.3.topi.unadj)  # 95% CIs 
      
    ## g) Adjusted: methlyation by birth to 3 months topi density
      adult.birth.3.topi.adj <- glm(methylation ~ topi.birth.3 + topi.gest + 
                                    topi.peri.concpt + sex + 
                                    age.mon + samp_year_cnt, 
                                  data = luma_data_adult)
      
    ## h) Parameter estimates
      summary(adult.birth.3.topi.adj)  # model parameter estimates
      confint(adult.birth.3.topi.adj)  # 95% CIs 
      
    ## i) Unadjusted: methlyation by birth to 3 months gnu density
      adult.birth.3.gnu.unadj <- glm(methylation ~ gnu.birth.3, 
                                   data = luma_data_adult)
      
    ## j) Parameter estimates
      summary(adult.birth.3.gnu.unadj)  # model parameter estimates
      confint(adult.birth.3.gnu.unadj)  # 95% CIs 
      
    ## k) Adjusted: methlyation by birth to 3 months gnu density
      adult.birth.3.gnu.adj <- glm(methylation ~ gnu.birth.3 + gnu.gest + 
                                   gnu.peri.concpt + sex + 
                                   age.mon + samp_year_cnt, 
                                 data = luma_data_adult)
      
    ## l) Parameter estimates
      summary(adult.birth.3.gnu.adj)  # model parameter estimates
      confint(adult.birth.3.gnu.adj)  # 95% CIs 
      
    ## m) Unadjusted: methlyation by birth to 3 months zebra density
      adult.birth.3.zebra.unadj <- glm(methylation ~ zebra.birth.3, 
                                     data = luma_data_adult)
      
    ## n) Parameter estimates
      summary(adult.birth.3.zebra.unadj)  # model parameter estimates
      confint(adult.birth.3.zebra.unadj)  # 95% CIs 
      
    ## o) Adjusted: methlyation by birth to 3 months zebra density
      adult.birth.3.zebra.adj <- glm(methylation ~ zebra.birth.3 + zebra.gest + 
                                     zebra.peri.concpt + sex + 
                                     age.mon + samp_year_cnt, 
                                   data = luma_data_adult)
      
    ## p) Parameter estimates
      summary(adult.birth.3.zebra.adj)  # model parameter estimates
      confint(adult.birth.3.zebra.adj)  # 95% CIs 
      
      
      
  ### 10.9 adult model: methylation by 3 to 6 months prey density     
    ## a) Unadjusted: methlyation by 3 to 6 months thomsons density
      adult.3.6.thomsons.unadj <- glm(methylation ~ thomsons.3.6, 
                                    data = luma_data_adult)
      
    ## b) Parameter estimates
      summary(adult.3.6.thomsons.unadj)  # model parameter estimates
      confint(adult.3.6.thomsons.unadj)  # 95% CIs 
      
    ## c) Adjusted: methlyation by 3 to 6 months thomsons density
      adult.3.6.thomsons.adj <- glm(methylation ~ thomsons.3.6 + 
                                    thomsons.birth.3 + thomsons.gest + 
                                    thomsons.peri.concpt + sex + 
                                    age.mon + samp_year_cnt, 
                                  data = luma_data_adult)
      
    ## d) Parameter estimates
      summary(adult.3.6.thomsons.adj)  # model parameter estimates
      confint(adult.3.6.thomsons.adj)  # 95% CIs 
      
    ## e) Unadjusted: methlyation by 3 to 6 months topi density
      adult.3.6.topi.unadj <- glm(methylation ~ topi.3.6, 
                                data = luma_data_adult)
      
    ## f) Parameter estimates
      summary(adult.3.6.topi.unadj)  # model parameter estimates
      confint(adult.3.6.topi.unadj)  # 95% CIs 
      
    ## g) Adjusted: methlyation by 3 to 6 months topi density
      adult.3.6.topi.adj <- glm(methylation ~ topi.3.6 + topi.birth.3 + 
                                topi.gest + topi.peri.concpt + sex + 
                                age.mon + samp_year_cnt, 
                              data = luma_data_adult)
      
    ## h) Parameter estimates
      summary(adult.3.6.topi.adj)  # model parameter estimates
      confint(adult.3.6.topi.adj)  # 95% CIs 
      
    ## i) Unadjusted: methlyation by 3 to 6 months gnu density
      adult.3.6.gnu.unadj <- glm(methylation ~ gnu.3.6, 
                               data = luma_data_adult)
      
    ## j) Parameter estimates
      summary(adult.3.6.gnu.unadj)  # model parameter estimates
      confint(adult.3.6.gnu.unadj)  # 95% CIs 
      
    ## k) Adjusted: methlyation by 3 to 6 months gnu density
      adult.3.6.gnu.adj <- glm(methylation ~ gnu.3.6 + gnu.birth.3 + 
                               gnu.gest + gnu.peri.concpt + sex + 
                               age.mon + samp_year_cnt, 
                             data = luma_data_adult)
      
    ## l) Parameter estimates
      summary(adult.3.6.gnu.adj)  # model parameter estimates
      confint(adult.3.6.gnu.adj)  # 95% CIs 
      
    ## m) Unadjusted: methlyation by 3 to 6 months zebra density
      adult.3.6.zebra.unadj <- glm(methylation ~ zebra.3.6, 
                                 data = luma_data_adult)
      
    ## n) Parameter estimates
      summary(adult.3.6.zebra.unadj)  # model parameter estimates
      confint(adult.3.6.zebra.unadj)  # 95% CIs 
      
    ## o) Adjusted: methlyation by 3 to 6 months zebra density
      adult.3.6.zebra.adj <- glm(methylation ~ zebra.3.6 + zebra.birth.3 + 
                                 zebra.gest + zebra.peri.concpt + sex + 
                                 age.mon + samp_year_cnt, 
                               data = luma_data_adult)
      
    ## p) Parameter estimates
      summary(adult.3.6.zebra.adj)  # model parameter estimates
      confint(adult.3.6.zebra.adj)  # 95% CIs 
      
      
  ### 10.10 adult model: methylation by 6 to 9 months prey density     
    ## a) Unadjusted: methlyation by 6 to 9 months thomsons density
      adult.6.9.thomsons.unadj <- glm(methylation ~ thomsons.6.9, 
                                    data = luma_data_adult)
      
    ## b) Parameter estimates
      summary(adult.6.9.thomsons.unadj)  # model parameter estimates
      confint(adult.6.9.thomsons.unadj)  # 95% CIs 
      
    ## c) Adjusted: methlyation by 6 to 9 months thomsons density
      adult.6.9.thomsons.adj <- glm(methylation ~ thomsons.6.9 + 
                                    thomsons.3.6 + thomsons.birth.3 + 
                                    thomsons.gest + thomsons.peri.concpt + 
                                    sex + age.mon + samp_year_cnt, 
                                  data = luma_data_adult)
      
    ## d) Parameter estimates
      summary(adult.6.9.thomsons.adj)  # model parameter estimates
      confint(adult.6.9.thomsons.adj)  # 95% CIs 
      
    ## e) Unadjusted: methlyation by 6 to 9 months topi density
      adult.6.9.topi.unadj <- glm(methylation ~ topi.6.9, 
                                data = luma_data_adult)
      
    ## f) Parameter estimates
      summary(adult.6.9.topi.unadj)  # model parameter estimates
      confint(adult.6.9.topi.unadj)  # 95% CIs 
      
    ## g) Adjusted: methlyation by 6 to 9 months topi density
      adult.6.9.topi.adj <- glm(methylation ~ topi.6.9 + topi.3.6 + 
                                topi.birth.3 + topi.gest + topi.peri.concpt + 
                                sex + age.mon + samp_year_cnt, 
                              data = luma_data_adult)
      
    ## h) Parameter estimates
      summary(adult.6.9.topi.adj)  # model parameter estimates
      confint(adult.6.9.topi.adj)  # 95% CIs 
      
    ## i) Unadjusted: methlyation by 6 to 9 months gnu density
      adult.6.9.gnu.unadj <- glm(methylation ~ gnu.6.9, 
                               data = luma_data_adult)
      
    ## j) Parameter estimates
      summary(adult.6.9.gnu.unadj)  # model parameter estimates
      confint(adult.6.9.gnu.unadj)  # 95% CIs 
      
    ## k) Adjusted: methlyation by 6 to 9 months gnu density
      adult.6.9.gnu.adj <- glm(methylation ~ gnu.6.9 + gnu.3.6 + 
                               gnu.birth.3 + gnu.gest + gnu.peri.concpt + 
                               sex + age.mon + samp_year_cnt, 
                             data = luma_data_adult)
      
    ## l) Parameter estimates
      summary(adult.6.9.gnu.adj)  # model parameter estimates
      confint(adult.6.9.gnu.adj)  # 95% CIs 
      
    ## m) Unadjusted: methlyation by 6 to 9 months zebra density
      adult.6.9.zebra.unadj <- glm(methylation ~ zebra.6.9, 
                                 data = luma_data_adult)
      
    ## n) Parameter estimates
      summary(adult.6.9.zebra.unadj)  # model parameter estimates
      confint(adult.6.9.zebra.unadj)  # 95% CIs 
      
    ## o) Adjusted: methlyation by 6 to 9 months zebra density
      adult.6.9.zebra.adj <- glm(methylation ~ zebra.6.9 + zebra.3.6 + 
                                 zebra.birth.3 + zebra.gest + 
                                 zebra.peri.concpt + sex + 
                                 age.mon + samp_year_cnt, 
                               data = luma_data_adult)
      
    ## p) Parameter estimates
      summary(adult.6.9.zebra.adj)  # model parameter estimates
      confint(adult.6.9.zebra.adj)  # 95% CIs 
      





      
      
      
      
      
   