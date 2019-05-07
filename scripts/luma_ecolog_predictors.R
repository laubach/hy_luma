###############################################################################
##############       Spotted Hyena Global DNA Methylation        ##############
##############             LUMA Ecological Predictors            ##############
##############                 By: Zach Laubach                  ##############
##############            last updated: 10 Sept 2018             ##############
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
    # 11: Save data tables as .csv



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
        
      # Check for dotwhisker and install if not already installed
        if (!'dotwhisker' %in% installed.packages()[,1]){
          install.packages ('dotwhisker')
        }
      # load dotwhisker packages
        library ('dotwhisker')  
        
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
        # includes vif function
        if (!'car' %in% installed.packages()[,1]){
          install.packages ('car')
        }
      # load car packages
        library ('car')

      # Check for lmtest and install if not already installed
        if (!'lmtest' %in% installed.packages()[,1]){
          install.packages ('lmtest')
        }
      # load lmtest package (used to test heteroskedacity)
        library ('lmtest')
        
      # Check for sandwich and install if not already installed
        if (!'sandwich' %in% installed.packages()[,1]){
          install.packages ('sandwich')
        }
      # load sandwich package (used generate Robust SE)
        library ('sandwich')  
      
        

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
      luma_data_path <- paste("~/R/R_wd/fisi/project/3_hy_GR_global_DNA_meth/",
                              "LUMA/output/", sep = '')
    
    ## b) The path to prey data
      prey_data_path <- paste("~/Git/fisi_lab/hy_prey_density/output/",
                              sep = '')
      
    ## c) The path to static Access Fisi backend
      acess_fisi_data_path <- paste("~/R/R_wd/fisi/project/", 
                                    "3_hy_GR_global_DNA_meth/",
                                    "LUMA/soc_eco_detrmnts_ms/", 
                                    "1_output_tidy_tbls/",
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
     
       
  ### 2.3 Import Access fisi backend # Live backend...swithc to static version
          # after analyses finalized
    ## a) read in tidy Access fisi backend tables and save as data frames
#      source(paste0("/Volumes/Holekamp/code_repository/R/1_output_tidy_tbls/",
#                   "load_tidy_tbls.R"))
      
    ## b) manually load tblFemalerank; the one from merge_databases 
      tblFemalerank <- read_csv(paste0(acess_fisi_data_path,
                                        "tblFemalerank.csv", sep = ''))
    

    
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
        mutate(grad.age.mon = interval(birthdate, den.grad) %/% days(1))
      # and check what mean age at graduation
          luma_data %>%
            summarise(avg.grad.age = round((mean (grad.age.mon,
                                                  na.rm = T)/30.44), 1),
                      max.grad.age = (max(grad.age.mon, na.rm = T)/30.44))
  
        luma_data <- luma_data %>%
            mutate(wean.age.mon = interval(birthdate, weaned) %/% days(1))
      # and check what mean age at graduation
        luma_data %>%
          summarise(avg.wean.age = round((mean(wean.age.mon,
                                               na.rm = T)/30.44), 1),
                    max.wean.age = (max(wean.age.mon, na.rm = T)/30.44)) 
        
    ## NOTE: To convert from age in days to age in months divde by average 
        # number of days per month  (30.44)  
        
    ## e) Create an estimated age in months by subtracting birthdate from
      # darting.date using lubridate
      luma_data <- luma_data %>%
        mutate(age.mon = round((interval(birthdate, 
                                         darting.date) %/% days(1)/ 30.44), 1))
      
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
                                     age.mon <=24 ~ c("subadult"),
                                   sex == "f" & age.mon > 24 ~ c("adult")))
      
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
                                      luma_data$clan %in% c("serena.n", 
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
      
    ## o) mark samples extracted via PAX-gene tubes those collected 2016 
      # onwards
      luma_data <- luma_data  %>%
        mutate(dna_extrct_mth = case_when(samp_year >= 2016 ~ c("pax"),
                                          samp_year < 2016 ~ c("old")))
                        
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
                                                 "Q4")))
     
    #**hack**# Fill in two missing maternal rank quartiles
      # change mom.absrank to an integer
      luma_data$mom.absrank <- as.integer(as.character
                                          (luma_data$mom.absrank))
      # fill in missing mom.strank.quart info based on mom.absrank for serena
      # hyneas
      luma_data$mom.strank.quart <- ifelse(
        (is.na(luma_data$mom.strank.quart) & 
           luma_data$mom.absrank > 18),
        "Q1 (lowest)", as.character(luma_data$mom.strank.quart))  
      
      
  ### 3.5 Left join tblFemalerank to luma_data
    ## a) append to luma_data each hyena's own rank the year that 
      # they were darted
      luma_data <- sqldf("SELECT
                         luma_data.*           
                         , tblFemalerank. absrank, strank
                         FROM luma_data      
                         LEFT JOIN tblFemalerank       
                         ON tblFemalerank.mom = luma_data.id
                         AND tblFemalerank.rank_year = luma_data.samp_year")   
      
    ## b) Create quartiles of a hyena's own rank
      # use use dplyr
        luma_data <- luma_data %>% 
          mutate(strank.quart.order = ntile(strank, 4))
      
      ## e) Create a nominal factor and rename and re-order the levels to 
      # sets the reference level for lowest rank 
      luma_data <- transform (luma_data, 
                             strank.quart = 
                                factor(strank.quart.order,
                                       levels = c(1, 2, 3, 4),
                                       labels= c("Q1 (lowest)",
                                                 "Q2","Q3",
                                                 "Q4")))
      
      
  ### 3.6 Reduce and group luma_data
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
                     age.mon = round(mean(age.mon),1),
                     rank_year = (first(rank_year)),
                     samp_year = (first(samp_year)),
                     hum.pres = (first(hum.pres)),
                     dna_extrct_mth = first(dna_extrct_mth),
                     absrank = as.factor(first(absrank)),
                     strank = as.numeric(first(strank)),
                     strank.quart = as.factor(first(strank.quart)),
                     strank.quart.order = (first(strank.quart.order)))
                    

                     
    ## b) Ungroup the data frame
        # dplyr retains grouping after creation of data frame that uses 
        # group_by
          luma_data_group <- ungroup(luma_data_group)
      
    
  ### 3.7 Clean prey density data and combine with LUMA data
    ## a) Select a subset of the prey_density data by column names
          
          prey_density <- prey_density %>%
            select(grep("prim", names(prey_density)), # contains 'thomsons'
                   #grep("topi", names(prey_density)), # contains 'topi'
                   #grep("gnu", names(prey_density)), # contains 'gnu'
                   #grep("zebra", names(prey_density)), # contains 'zebra'
                   grep("num", names(prey_density)), # contains 'num'
                   grep("^id$", names(prey_density))) %>% # exact match 'id'
            select (-c(number.littermates)) #%>%
            #select (-starts_with("thom.topi.")) %>%
            #select (-starts_with("gnu.zebra."))
          
    ## b) Left join prey_density to luma_data   
        luma_data <- left_join(luma_data,
                               prey_density, by = "id")
        
    ## c) Left join prey_density to luma_data_group   
        luma_data_group <- left_join(luma_data_group,
                              prey_density, by = "id")
  
  
  ## 3.8 Re-code variables to the appropriate class
    ## a) Re-code *nominal* factor (with ordered levels)  
        # Set levels (odering) age.cat variable and sets the reference level 
        # to cub 
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
#        luma_data <- luma_data %>%
#          mutate(age.ordin = as.numeric(case_when(age.mon <= 12 ~ c(11.0),
#                                                  age.mon > 12 & 
#                                                    age.mon <=24 ~ c(17.5),
#                                                  age.mon > 24 ~ c(51.0))))
        
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
        
    ## f) Re-code mom.strank.quart as nominal factor and set level (order)
      luma_data <- transform(luma_data,
                             mom.strank.quart = factor(mom.strank.quart, 
                                               levels = c("Q1 (lowest)",
                                                          "Q2","Q3",
                                                          "Q4")))
        
      luma_data_group <- transform(luma_data_group,
                                   mom.strank.quart = factor(mom.strank.quart,
                                                levels = c("Q1 (lowest)",
                                                           "Q2","Q3",
                                                           "Q4")))
        
             
  ### 3.9 Check the effect of extraction method...literature suggests this
      # matters and looking at data there are striking difference amomg
      # six samples extracted from Pax gene tubes
    ## a) uses 'nmle' package, which will provided p-value estimates
        dna.mtd.lme <- lme(methylation ~ dna_extrct_mth, random =~1|id, 
                           subset(luma_data_group,!is.na(x = dna_extrct_mth)))
        
        summary( dna.mtd.lme)        # model summary
        intervals( dna.mtd.lme, 
                   which = "fixed")  # 95% CIs 
        
    ## b) remove PAX gene extractions from data
        luma_data <- luma_data %>% 
          filter (!grepl("pax", dna_extrct_mth))
        
        luma_data_group <- luma_data_group %>% 
          filter (!grepl("pax", dna_extrct_mth))
        
    ## c) Calculate the number of individual hyeans and repeat samples
        luma_data_hy_id <- luma_data %>% 
          filter (!is.na(methylation))%>% # check/remove rows where meth NA
          group_by (id) %>% # set grouping same ID within same cat age
          summarise (age.reps = sum(!is.na(methylation))) # n per ID w/in age 
        # class   
      
###############################################################################
##############               4. UniVariate analyses              ##############
###############################################################################      

  ### 4.1 Methylation (outcome) univariate
    ## a) Descriptive stats methylation
      # calculate the mean, median and standard deviation of % methylation
      # NOTE: uses luma_data_group data frame
        univar_meth <- luma_data_group %>%
          summarize (n.id = n_distinct(id),
                     n.samp = sum(!is.na(methylation)),
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
      grid.table(hum_pres_summary)
      dev.off() 
      
      
  ### 4.7 Descriptive stats prey density            
    ## a) Prey density summary 
      prey_peri_summary <- luma_data_group %>%
        summarise (dev.period = print("peri.concpt"),
                   n.prey = sum(!is.na(prim.prey.peri.concpt)),
                   avg.prim.prey = round (mean(prim.prey.peri.concpt, 
                                               na.rm = T),2),
                   stdev.prim.prey = round (sd(prim.prey.peri.concpt, 
                                               na.rm = T), 2))
                   
      prey_gest_summary <-  luma_data_group %>%
        summarise (dev.period = print("gest"),
                   n.prey = sum(!is.na(prim.prey.gest)),
                   avg.prim.prey = round (mean(prim.prey.gest, 
                                              na.rm = T),2),
                   stdev.prim.prey = round (sd(prim.prey.gest, 
                                              na.rm = T), 2))
      
      prey_birth.3_summary <- luma_data_group %>%
        summarise (dev.period = print("birth.3"),
                   n.prey = sum(!is.na(prim.prey.birth.3)),
                   avg.prim.prey = round (mean(prim.prey.birth.3, 
                                              na.rm = T),2),
                   stdev.prim.prey = round (sd(prim.prey.birth.3, 
                                              na.rm = T), 2))
                   
      prey_3.6_summary <- luma_data_group %>%
        summarise (dev.period = print("3.6"),
                   n.prey = sum(!is.na(prim.prey.3.6)),
                   avg.prim.prey = round (mean(prim.prey.3.6, 
                                              na.rm = T),2),
                   stdev.prim.prey = round (sd(prim.prey.3.6, 
                                              na.rm = T), 2))
                   
      prey_6.9_summary <- luma_data_group %>%
        summarise (dev.period = print("6.9"),
                   n.prey = sum(!is.na(prim.prey.6.9)),
                   avg.prim.prey = round (mean(prim.prey.6.9, 
                                              na.rm = T),2),
                   stdev.prim.prey = round (sd(prim.prey.6.9, 
                                              na.rm = T), 2))
      
    ## b) combine the prey density descriptive stats into a single data frame
      prey_var_summary <- rbind(prey_peri_summary, prey_gest_summary,
                                prey_birth.3_summary, prey_3.6_summary,
                                prey_6.9_summary)


    ## c) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/output_luma_ecolog/",
                 "hum_pres_summary.pdf"),
      height = 3, width = 5)
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
        #var_names <- c("total","thomsons", "topi", "gnu", "zebra")
        var_names <- c("prim.prey")
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
        
      summary(samp.year.lme)        # model summary
      intervals(samp.year.lme, 
                  which = "fixed")  # 95% CIs 
      
      
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
        labs(title = "Percent Global DNA Mehtylation by Sex") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(legend.position = "none") + # remove legend
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = "white")) +
        # add major axes
        theme(axis.line = element_line(colour = "darkgrey", 
                                        size = 1, linetype = "solid")) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face="bold", color="black", 
                                         size=18, angle=0,
                                         margin = margin(t = 0, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face="bold", color="black", 
                                         size=18, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10))) +
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
   
      #summary(sex.lm)     # model summary
      #confint(sex.lm)     # 95% CIs
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
                   median =  round (quantile(methylation, c(.5), 
                                             na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/output_luma_ecolog/meth_by_age.pdf"), 
          height = 6, width = 7)
      grid.table(meth_by_age)
      dev.off()
    
    ## c) Plot Methylation by Age
      # graph of the raw data for percent global DNA methylaiton by age 
      ggplot(data = subset(luma_data_group, !is.na(x = age.cat)), 
             aes(x = age.cat, y = methylation, color = age.cat)) + 
        geom_boxplot() +
        theme(text = element_text(size=20))+
        labs(title = "Percent Global DNA Mehtylation by Age") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(legend.position = "none") + # remove legend
        #theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = "white")) +
        # add major axes
        #theme(axis.line = element_line(colour = "darkgrey", 
        #                               size = 1, linetype = "solid")) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face="bold", color="black", 
                                         size=18, angle=0,
                                         margin = margin(t = 0, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face="bold", color="black", 
                                         size=18, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10))) +
        ylab("% Global DNA Methylation") +
        xlab("Categorical Age")
      
    ## d) Save Plot
      # use ggsave to save the plot
      ggsave("meth_by_age_plot.pdf", plot = last_plot(), device = NULL, 
             path = paste0(here(),"/output/output_luma_ecolog"),
             scale = 1, width = 7, height = 5, 
             units = c("in"), dpi = 300, limitsize = TRUE)  
  
    ## e) Bivariate regression methylatino by age.cat
      # uses 'nmle' package, which will provided p-value estimates
      age.lme <- lme(methylation ~ age.cat, random =~1|id, 
                     subset(luma_data_group,!is.na(x = age.cat)))
      
      summary(age.lme)            # model summary
      intervals(age.lme, 
                which = "fixed")  # 95% CIs 
      # generate p-value from type II Wald test
      car::Anova(age.lme,Type ="II", test = "Wald") 
      #anova.lme(age.lme)      # generate p-value from Wald test
      
    ## f) Bivariate regression methylatino by age.mon
      # uses 'nmle' package, which will provided p-value estimates
      age.mon.lme <- lme(methylation ~ age.mon, random =~1|id, 
                     subset(luma_data_group,!is.na(x = age.mon)))
      
      summary(age.mon.lme)        # model estimates
      intervals(age.mon.lme, 
                which = "fixed")  # 95% CIs 

  ### 6.3 Bivariate Statistics Methylation by Rank
    ## a) Graph of the raw data for percent global DNA methylaiton by maternal 
      # rank 
      ggplot(data = subset(luma_data_group, !is.na(x = mom.strank)),
             aes(x = mom.strank, y = methylation)) +
        geom_point(shape = 1) +
        theme(text = element_text(size=20))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        geom_smooth(method = loess, se = F) + # Add smooth curve best fit lines
        labs(title = "Percent Global DNA Mehtylation by Maternal Rank",
             fill = "age") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(legend.position = "none") + # remove legend
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = "white")) +
        # add major axes
        theme(axis.line = element_line(colour = "darkgrey", 
                                       size = 1, linetype = "solid")) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face="bold", color="black", 
                                         size=18, angle=0,
                                         margin = margin(t = 0, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face="bold", color="black", 
                                         size=18, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10))) +
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
      
    ## c) Bivariate Regression Methylatino by maternal rank 
      # uses 'nmle' package, which will provided p-value estimates
      mom.rank.lme <- lme(methylation ~ mom.strank.quart, random =~1|id, 
                     subset(luma_data_group,!is.na(x = mom.strank.quart)))
      
      summary(mom.rank.lme)       # model summary
      intervals(mom.rank.lme, 
                which = "fixed")  # 95% CIs
      # generate p-value from type II Wald test
      car::Anova(mom.rank.lme,Type ="II", test = "Wald") 
  

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
             aes(x = lit.size, y = methylation, color = lit.size)) + 
        geom_boxplot() +
        scale_color_manual(values=c("black", "darkgrey")) +
        theme(text = element_text(size=20))+
        labs(title = "%CCGG Methylation by Litter Size") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(legend.position = "none") + # remove legend
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = "white")) +
        # add major axes
        theme(axis.line = element_line(colour = "lightgrey", 
                                       size = 1, linetype = "solid")) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face="bold", color="black", 
                                         size=18, angle=0,
                                         margin = margin(t = 0, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face="bold", color="black", 
                                         size=18, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10))) +
        ylab("%CCGG Methylation") +
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
      
      summary(lit.size.lme)         # model summary
      intervals(lit.size.lme, 
                which = "fixed")    # 95% CIs 
      # generate p-value from type II Wald test
      car::Anova(lit.size.lme,Type ="II", test = "Wald") 
      
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
        scale_color_manual(values=c("lightgrey", "darkgrey", "black")) +
        theme(text = element_text(size=20))+
        labs(title = "%CCGG Methylation by Human Disturbance") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(legend.position = "none") + # remove legend
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = "white")) +
        # add major axes
        theme(axis.line = element_line(colour = "lightgrey", 
                                       size = 1, linetype = "solid")) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face="bold", color="black", 
                                         size=18, angle=0,
                                         margin = margin(t = 0, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face="bold", color="black", 
                                         size=18, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10))) +
        ylab("%CCGG Methylation") +
        xlab("Human Disturbance")
      
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
      
      summary(hum.pres.lme)        # model summary
      intervals(hum.pres.lme, 
                which = "fixed")   # 95% CIs 
      anova(hum.pres.lme)          # generate p-value from Wald test
      # generate p-value from type II Wald test
      car::Anova(hum.pres.lme,Type ="II", test = "Wald") 
    
  ### 6.6 Bivariate statistics methylation by periconceptional prey density
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg 
      # uses 'nmle' package, which will provided p-value estimates
      
      prim.peri.concpt.lme <- lme(methylation ~ prim.prey.peri.concpt, 
                                      random =~1|id, 
                                      subset(luma_data_group,
                                             !is.na(x = prim.prey.peri.concpt)))
      
      summary(prim.peri.concpt.lme)   # model summary
      intervals(prim.peri.concpt.lme, 
                which = "fixed")      # 95% CIs

    
  ### 6.7 Bivariate statistics methylation by gestational prey density
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg 
      # uses 'nmle' package, which will provided p-value estimates
      prim.gest.lme <- lme(methylation ~ prim.prey.gest, 
                                      random =~1|id, 
                                      subset(luma_data_group,
                                             !is.na(x = prim.prey.gest)))
      
      summary(prim.gest.lme)      # model summary
      intervals(prim.gest.lme, 
                which = "fixed")  # 95% CIs
      
      
  ### 6.8 Bivariate statistics methylation by birth to 3 months prey density
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg
      prim.birth.3.lme <- lme(methylation ~ prim.prey.birth.3, 
                                      random =~1|id, 
                                      subset(luma_data_group,
                                             !is.na(x = prim.prey.birth.3)))
      
      summary(prim.birth.3.lme)   # model summary
      intervals(prim.birth.3.lme, 
                which = "fixed")  # 95% CIs
      
  
  ### 6.9 Bivariate statistics methylation by 3 to 6 months prey density
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg 
      # uses 'nmle' package, which will provided p-value estimates
      prim.3.6.lme <- lme(methylation ~ prim.prey.3.6, 
                                        random =~1|id, 
                                        subset(luma_data_group,
                                               !is.na(x = prim.prey.3.6)))
      
      summary(prim.3.6.lme)         # model summary
      intervals(prim.3.6.lme, 
                which = "fixed")    # 95% CIs
      
    
      
  ### 6.10 Bivariate statistics methylation by 6 to 9 months prey density
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg 
      # uses 'nmle' package, which will provided p-value estimates
      prim.6.9.lme <- lme(methylation ~ prim.prey.6.9, 
                                      random =~1|id, 
                                      subset(luma_data_group,
                                             !is.na(x = prim.prey.6.9)))
      
      summary(prim.6.9.lme)       # model summary
      intervals(prim.6.9.lme, 
                which = "fixed")  # 95% CIs
      
     
      
      
###############################################################################
##############  7. Test for effect modification and re-tidy data ##############
###############################################################################       
    
  ### 7.1 Test for effect modification of maternal rank and methylation by sex    
    ## a) Model uses 'nmle' package, which will provided p-value estimates
      rank_by_sex.lme <- lme(methylation ~ sex * mom.strank, random =~1|id, 
                     subset(luma_data_group,!is.na(x = mom.strank)))
      
    ## b) Parameter estimates  
      summary(rank_by_sex.lme) 
      intervals(rank_by_sex.lme, which = "fixed")
      
      
  ### 7.2 Test for effect modification of maternal rank and methylation by age    
    ## a) Model uses 'nmle' package, which will provided p-value estimates
      rank_by_age.lme <- lme(methylation ~ age.cat * mom.strank, random =~1|id, 
                             subset(luma_data_group,!is.na(x = mom.strank)))
      
    ## b) Parameter estimates   
      summary(rank_by_age.lme) 
      intervals(rank_by_age.lme, which = "fixed")    
      # generate p-value from type II Wald test
      car::Anova(rank_by_age.lme,Type ="II", test = "Wald") 
      

  #*** RESUBMISSION 2 *** addressing reviewer's comments
    ## c) Model uses 'nmle' package, which will provided p-value estimates
      rank_by_age_adj.lme <- lme(methylation ~ age.cat * mom.strank + sex, 
                             random =~1|id, 
                             subset(luma_data_group,!is.na(x = mom.strank)))
      
    ## d) Parameter estimates   
      summary(rank_by_age_adj.lme) 
      intervals(rank_by_age_adj.lme, which = "fixed")    
      # generate p-value from type II Wald test
      car::Anova(rank_by_age_adj.lme,Type ="II", test = "Wald") 
      
    ## e) Model uses 'nmle' package, which will provide p-value estimates
      rank_by_cat_age_adj.lme <- lme(methylation ~ age.cat * mom.strank.quart 
                                     + sex, 
                                 random =~1|id, 
                                 subset(luma_data_group,!is.na(x = mom.strank)))
      
    ## f) Parameter estimates   
      summary(rank_by_cat_age_adj.lme) 
      intervals(rank_by_cat_age_adj.lme, which = "fixed")    
      # generate p-value from type II Wald test
      car::Anova(rank_by_cat_age_adj.lme,Type ="II", test = "Wald") 
      
      
  ### 7.3 View methylation by rank within age strata
      
    ## a) Graph of the raw data for percent global DNA methylaiton by maternal 
      # rank stratified by age because age seems to be important (see above)
      ggplot(data = subset(luma_data_group, !is.na(x = mom.strank)),
             aes(x = mom.strank, y = methylation, color = age.cat)) +
        geom_point(shape = 1) +
        theme(text = element_text(size=18))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        geom_smooth(method = loess, se = F) + # Add smooth curve best fit lines
        labs(title = "%CCGG methylation by maternal rank 
             and stratified by age",
             fill = "age") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(legend.position = c(0.007, 0.01),
              legend.justification = c(0, 0), 
              legend.background = element_rect(colour="grey80"),
              legend.title = element_blank()) + # remove legend
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = "white")) +
        # add major axes
        #theme(axis.line = element_line(colour = "lightgrey", 
        #                               size = 1, linetype = "solid")) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face="bold", color="black", 
                                         size=18, angle=0,
                                         margin = margin(t = 0, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face="bold", color="black", 
                                         size=18, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10))) +
        ylab("%CCGG Methylation") +
        xlab("Maternal Rank")
      
    ## b) Save Plot
      # use ggsave to save the linearization plot
      ggsave("meth_by_rank_by_age_loess_plot.pdf", plot = last_plot(), 
             device = NULL,
             path = paste0(here(),"/output/output_luma_ecolog"), 
             scale = 1, width = 7, height = 5,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## c) Summary stats methylation by age 
      # NOTE: uses luma_data_group; first average over ID within age cat. 
      # and then take avg
      meth_by_mom_rank <- luma_data_group %>%
        group_by (mom.strank.quart) %>%
        summarise (n.id = n(),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## d) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/output_luma_ecolog/meth_by_mom_rank.pdf"), 
          height = 6, width = 7)
      grid.table(meth_by_mom_rank)
      dev.off()
      
    ## e) Graph of the raw data for percent global DNA methylaiton by maternal 
      # rank stratified by age
      ggplot(subset(luma_data_group, !is.na(x = mom.strank)),
             aes(x = mom.strank.quart.order, y = methylation, 
                 color = age.cat )) +
        geom_point(shape = 1) +
        geom_smooth(method = loess, se = F) + # Add loess fit lines
        theme(text = element_text(size=20))+
        labs(title = "%CCGG Methylation by 
Maternal Rank by Age",
             fill = "age") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(legend.position = "none") + # remove legend
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = "white")) +
        # add major axes
        theme(axis.line = element_line(colour = "lightgrey", 
                                       size = 1, linetype = "solid")) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face="bold", color="black", 
                                         size=18, angle=0,
                                         margin = margin(t = 0, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face="bold", color="black", 
                                         size=18, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10))) +
        ylab("%CCGG Methylation") +
        xlab("Maternal Rank")
      
    ## f) Save Plot
      # use ggsave to save the linearization plot
      ggsave("meth_by_rank_quart_by_age_loess_plot.pdf", plot = last_plot(), 
             device = NULL,
             path = paste0(here(),"/output/output_luma_ecolog"), 
             scale = 1, width = 7, height = 5,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## g) Box-plots of methylation by rank within age categories  
      ggplot(subset(luma_data_group, 
                    !is.na(x = mom.strank.quart)&!is.na(x=age.cat)),
             aes(y = methylation, x = factor(mom.strank.quart), 
                 fill = age.cat))+
        geom_boxplot() +
        theme(text = element_text(size=20))+
        labs(title = "%CCGG Methylation by Maternal Rank by Age",
             fill = "age") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(legend.position = "none") + # remove legend
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = "white")) +
        # add major axes
        theme(axis.line = element_line(colour = "lightgrey", 
                                       size = 1, linetype = "solid")) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face="bold", color="black", 
                                         size=18, angle=0,
                                         margin = margin(t = 0, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face="bold", color="black", 
                                         size=18, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10))) +
        ylab("%CCGG Methylation") +
        xlab("Maternal Rank")
      
    ## h) Save Plot
      # use ggsave to save the linearization plot
      ggsave("meth_by_rank_by_age_boxplot.pdf", plot = last_plot(), 
             device = NULL,
             path = paste0(here(),"/output/output_luma_ecolog"),
             scale = 1, width = 7,height = 5,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
      
  ### 7.4 Stratify data by life history variables
    ## a) Cub subset that includes both females and males
      luma_data_cub <- luma_data_group %>%
        filter(grepl('^cub$', age.cat))
      
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
      
      plot(cub.mom.rank.unadj) # view the residual and QQ plots
      
    ## d) Adjusted: methlyation by mom.strank.quart
      cub.mom.rank.adj <- glm(methylation ~ mom.strank.quart + sex + 
                                    age.mon,
                              data = luma_data_cub)
      
    ## e) Parameter estimates
      summary(cub.mom.rank.adj)  # print model summary, effects and SE
      confint(cub.mom.rank.adj)  # print 95% CIs for parameter estimates
    
    ## f) Check for heteroskedacity, normality, and outliers    
      plot(cub.mom.rank.adj) # view the residual and QQ plots
      bptest(cub.mom.rank.adj) # heteroskedasticity using the Breusch-Pagan 
                               # test; if significant, then generate test
                               # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(cub.mom.rank.adj, vcov = vcovHC(cub.mom.rank.adj)) 
      
    ## g) Sensitivity: methlyation by mom.strank.quart
#      cub.mom.rank.sens <- glm(methylation ~ mom.strank.quart + sex + 
#                                age.mon + samp_year_cnt,
#                              data = luma_data_cub)
      
    ## h) Parameter estimates
#      summary(cub.mom.rank.sens)  # print model summary, effects and SE
#      confint(cub.mom.rank.sens)  # print 95% CIs for parameter estimates

      
      
    ## i) Extract mom.strank.quart estimates and 
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
      
    ## j) Graph cub.mom.rank effects
      ggplot(cub.rank.ef.table, aes(x = mom.strank.quart, y = fit)) +
        geom_point() +
        geom_errorbar(aes(ymin= fit-se, ymax= fit+se), width=0.4) +
        theme(text = element_text(size=20))+
        labs(title = "Cub %CCGG Methylation 
by Maternal Rank") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(legend.position = "none") + # remove legend
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = "white")) +
        # add major axes
        theme(axis.line = element_line(colour = "lightgrey", 
                                       size = 1, linetype = "solid")) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face="bold", color="black", 
                                         size=18, angle=0,
                                         margin = margin(t = 0, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face="bold", color="black", 
                                         size=18, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10))) +
        ylab("%CCGG Methylation ± SE") +
        xlab("Maternal Rank")
    
    ## k) Save Plot
      # use ggsave to save the linearization plot
      ggsave("mat_rank_cub_mod_beta.pdf", plot = last_plot(), device = NULL,
             path = paste0(here(),"/output/output_luma_ecolog"), 
                           scale = 1, width = 7,
             height = 5.5,
             units = c("in"), dpi = 300, limitsize = TRUE)
 
    ## l) Do a post-hoc test to determine if Q2, Q3, and Q4 differ
      pairwise.t.test(luma_data_cub$methylation,
                      luma_data_cub$mom.strank.quart,
                      p.adj = "none")
      TukeyHSD(aov(methylation ~ mom.strank.quart + sex,
                   data = luma_data_cub))
    
    ## m) Combine quartiles 2-4 into a single category
      luma_data_cub$mom.strank.quart.comb <- as.factor(
      ifelse(luma_data_cub$mom.strank.quart ==  "Q1 (lowest)",
             "Q1 (lowest)", "Q2-Q4 (highest)"))
    
    ## n)  Adjusted: methlyation by mom.strank.quart binned
      cub.rank.adj2 <- glm(methylation ~ mom.strank.quart.comb + sex + 
                             age.mon,
                         data = luma_data_cub)
      
      summary(cub.rank.adj2)  # print model summary, effects and SE
      confint(cub.rank.adj2)  # print 95% CIs for parameter estimates
    
      
  ### 8.4 Cub model: methylation by litter size
    ## a) Check within strata descritpive stats
      luma_data_cub %>%
        group_by (lit.size) %>%
        summarise (n.id = sum(!is.na(lit.size)),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Unadjusted: methlyation by lit.size
      cub.lit.size.unadj <- glm(methylation ~  lit.size, data = luma_data_cub)
      
    ## c) Parameter estimates
      summary(cub.lit.size.unadj)  # model parameter estimates
      confint(cub.lit.size.unadj)  # 95% CIs
      Anova(cub.lit.size.unadj, Type ="II", test = "Wald") # Wald test p
      
      plot(cub.lit.size.unadj) # view the residual and QQ plots
      
    ## d) Adjusted: methlyation by lit.size
      cub.lit.size.adj <- glm(methylation ~  lit.size + age.mon + sex,
                             data = luma_data_cub)
      
    ## e) Parameter estimates
      summary(cub.lit.size.adj)  # model parameter estimates
      confint(cub.lit.size.adj)  # 95% CIs 
      
    ## f) Check for heteroskedacity, normality, and outliers    
      plot(cub.lit.size.adj) # view the residual and QQ plots
      bptest(cub.lit.size.adj) # heteroskedasticity using the Breusch-Pagan 
      # test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(cub.lit.size.adj, vcov = vcovHC(cub.lit.size.adj)) 
    
    ## g) Sensitivity: methlyation by lit.size
#      cub.lit.size.sens <- glm(methylation ~  lit.size + age.mon + sex +
#                                samp_year_cnt,
#                              data = luma_data_cub)
      
    ## h) Parameter estimates
#      summary(cub.lit.size.sens)  # model parameter estimates
#      confint(cub.lit.size.sens)  # 95% CIs 

    
  ### 8.5 Cub model: methylation by human presence proxy
    ## a) Check within strata descritpive stats
      luma_data_cub %>%
        group_by (hum.pres) %>%
        summarise (n.id = sum(!is.na(lit.size)),
                   avg = round (mean(methylation, na.rm = T), 2),
                   median =  round (quantile(methylation, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(methylation, na.rm = T), 2))
      
    ## b) Unadjusted: methlyation by hum.pres
      cub.hum.pres.unadj <- glm(methylation ~  hum.pres, data = luma_data_cub)
      
    ## c) Parameter estimates
      summary(cub.hum.pres.unadj)  # model parameter estimates
      confint(cub.hum.pres.unadj)  # 95% CIs
      Anova(cub.hum.pres.unadj, Type ="II", test = "Wald") # Wald test p
      
      plot(cub.hum.pres.unadj) # view the residual and QQ plots
      
    ## d) Adjusted: methlyation by hum.pres
      cub.hum.pres.adj <- glm(methylation ~  hum.pres + age.mon + sex,
                              data = luma_data_cub)
      
    ## e) Parameter estimates
      summary(cub.hum.pres.adj)  # model parameter estimates
      confint(cub.hum.pres.adj)  # 95% CIs 
      
    ## f) Check for heteroskedacity, normality, and outliers    
      plot(cub.hum.pres.adj) # view the residual and QQ plots
      bptest(cub.hum.pres.adj) # heteroskedasticity using the Breusch-Pagan 
      # test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(cub.lit.size.adj, vcov = vcovHC(cub.lit.size.adj)) 
    
    ## g) Sensitivity: methlyation by hum.pres
      cub.hum.pres.sens <- glm(methylation ~  hum.pres + age.mon + sex +
                                samp_year_cnt,
                              data = luma_data_cub)
      
    ## h) Parameter estimates
      summary(cub.hum.pres.sens)  # model parameter estimates
      confint(cub.hum.pres.sens)  # 95% CIs 
      
    ## i) Calculate VIF
      vif(cub.hum.pres.sens)
      
    ## j) Sensitivity: methlyation by samp_year_cnt
      cub.samp.age.sens <- glm(methylation ~ age.mon + sex +
                                 samp_year_cnt,
                               data = luma_data_cub)
      
    ## k) Parameter estimates
      summary(cub.samp.age.sens)  # model parameter estimates
      confint(cub.samp.age.sens)  # 95% CIs 
      
    ## l) Sensitivity: meth by mat_rank, while controling for samp_year_cnt
      cub.rank.samp.age.sens <- glm(methylation ~ mom.strank.quart + age.mon +
                                      sex + samp_year_cnt,
                               data = luma_data_cub)
      
    ## m) Parameter estimates
      summary(cub.rank.samp.age.sens)  # model parameter estimates
      confint(cub.rank.samp.age.sens)  # 95% CIs 
      
      
      
  ### 8.6 Cub model: methylation by periconceptional prey density     
    ## a) Unadjusted: methlyation by periconceptional prim.prey density
      cub.peri.prim.prey.unadj <- glm(methylation ~ prim.prey.peri.concpt, 
                               data = luma_data_cub)
      
    ## b) Parameter estimates
      summary(cub.peri.prim.prey.unadj)  # model parameter estimates
      confint(cub.peri.prim.prey.unadj)  # 95% CIs 
      
      plot(cub.peri.prim.prey.unadj) # view the residual and QQ plots
      
    ## c) Adjusted: methlyation by periconceptional prim.prey density
      cub.peri.prim.prey.adj <- glm(methylation ~ prim.prey.peri.concpt + sex + 
                                      age.mon, 
                                    data = luma_data_cub)
      
    ## d) Parameter estimates
      summary(cub.peri.prim.prey.adj)  # model parameter estimates
      confint(cub.peri.prim.prey.adj)  # 95% CIs 
      
    ## e) Check for heteroskedacity, normality, and outliers    
      plot(cub.peri.prim.prey.adj) # view the residual and QQ plots
      bptest(cub.peri.prim.prey.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(cub.peri.prim.prey.adj, vcov = vcovHC(cub.peri.prim.prey.adj))  
      
    ## f) Adjusted: methlyation by periconceptional prim.prey density
#      cub.peri.prim.prey.sens <- glm(methylation ~ prim.prey.peri.concpt + sex + 
#                                      age.mon + samp_year_cnt, 
#                                    data = luma_data_cub)
      
    ## g) Parameter estimates
#      summary(cub.peri.prim.prey.sens)  # model parameter estimates
#      confint(cub.peri.prim.prey.sens)  # 95% CIs 
      
      
  ### 8.7 Cub model: methylation by gestational prey density     
    ## a) Unadjusted: methlyation by gestational prim.prey density
      cub.gest.prim.prey.unadj <- glm(methylation ~ prim.prey.gest, 
                                     data = luma_data_cub)
      
    ## b) Parameter estimates
      summary(cub.gest.prim.prey.unadj)  # model parameter estimates
      confint(cub.gest.prim.prey.unadj)  # 95% CIs 
      
      plot(cub.gest.prim.prey.unadj) # view the residual and QQ plots
      
    ## c) Adjusted: methlyation by gestational prim.prey density
      cub.gest.prim.prey.adj <- glm(methylation ~ prim.prey.gest + 
                                     prim.prey.peri.concpt + sex + 
                                     age.mon, 
                                   data = luma_data_cub)
      
    ## d) Parameter estimates
      summary(cub.gest.prim.prey.adj)  # model parameter estimates
      confint(cub.gest.prim.prey.adj)  # 95% CIs 
  
    ## e) Check for heteroskedacity, normality, and outliers    
      plot(cub.gest.prey.adj) # view the residual and QQ plots
      bptest(cub.gest.prim.prey.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(cub.gest.prim.prey.adj, vcov = vcovHC(cub.gest.prim.prey.adj))  
      
    ## f) Sensitivity: methlyation by gestational prim.prey density
#      cub.gest.prim.prey.sens <- glm(methylation ~ prim.prey.gest + 
#                                      prim.prey.peri.concpt + sex + 
#                                      age.mon + samp_year_cnt, 
#                                    data = luma_data_cub)
      
    ## g) Parameter estimates
#      summary(cub.gest.prim.prey.sens)  # model parameter estimates
#      confint(cub.gest.prim.prey.sens)  # 95% CIs  
        
      
  ### 8.8 Cub model: methylation by birth to 3 months prey density     
    ## a) Unadjusted: methlyation by birth to 3 months prim.prey density
      cub.birth.3.prim.prey.unadj <- glm(methylation ~ prim.prey.birth.3, 
                                     data = luma_data_cub)
      
    ## b) Parameter estimates
      summary(cub.birth.3.prim.prey.unadj)  # model parameter estimates
      confint(cub.birth.3.prim.prey.unadj)  # 95% CIs 
      
      plot(cub.birth.3.prim.prey.unadj) # view the residual and QQ plots
      
    ## c) Adjusted: methlyation by birth to 3 months prim.prey density
      cub.birth.3.prim.prey.adj <- glm(methylation ~ prim.prey.birth.3 + 
                                        prim.prey.gest + prim.prey.peri.concpt + 
                                        sex + age.mon, 
                                   data = luma_data_cub)
      
    ## d) Parameter estimates
      summary(cub.birth.3.prim.prey.adj)  # model parameter estimates
      confint(cub.birth.3.prim.prey.adj)  # 95% CIs 
      
    ## e) Check for heteroskedacity, normality, and outliers    
      plot(cub.birth.3.prey.adj) # view the residual and QQ plots
      bptest(cub.birth.3.prim.prey.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(cub.birth.3.prim.prey.adj, 
               vcov = vcovHC(cub.birth.3.prim.prey.adj))
      
    ## f) Sensitivity: methlyation by birth to 3 months prim.prey density
#      cub.birth.3.prim.prey.sens <- glm(methylation ~ prim.prey.birth.3 + 
#                                         prim.prey.gest + prim.prey.peri.concpt + 
#                                         sex + age.mon + samp_year_cnt, 
#                                       data = luma_data_cub)
      
    ## g) Parameter estimates
#      summary(cub.birth.3.prim.prey.sens)  # model parameter estimates
#      confint(cub.birth.3.prim.prey.sens)  # 95% CIs 
      
      
  ### 8.9 Cub model: methylation by 3 to 6 months prey density     
    ## a) Unadjusted: methlyation by 3 to 6 months prim.prey density
      cub.3.6.prim.prey.unadj <- glm(methylation ~ prim.prey.3.6, 
                                     data = luma_data_cub)
      
    ## b) Parameter estimates
      summary(cub.3.6.prim.prey.unadj)  # model parameter estimates
      confint(cub.3.6.prim.prey.unadj)  # 95% CIs 
      
      plot(cub.3.6.prim.prey.unadj) # view the residual and QQ plots
      
    ## c) Adjusted: methlyation by 3 to 6 months prim.prey density
      cub.3.6.prim.prey.adj <- glm(methylation ~ prim.prey.3.6 + 
                                    prim.prey.birth.3 + prim.prey.gest + 
                                    prim.prey.peri.concpt + sex + 
                                     age.mon, 
                                   data = luma_data_cub)
      
    ## d) Parameter estimates
      summary(cub.3.6.prim.prey.adj)  # model parameter estimates
      confint(cub.3.6.prim.prey.adj)  # 95% CIs 
      
    ## e) Check for heteroskedacity, normality, and outliers    
      plot(cub.3.6.prey.adj) # view the residual and QQ plots
      bptest(cub.3.6.prim.prey.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(cub.3.6.prim.prey.adj, vcov = vcovHC(cub.3.6.prim.prey.adj)) 
    
    ## f) Sensitivity: methlyation by 3 to 6 months prim.prey density
#      cub.3.6.prim.prey.sens <- glm(methylation ~ prim.prey.3.6 + 
#                                     prim.prey.birth.3 + prim.prey.gest + 
#                                     prim.prey.peri.concpt + sex + 
#                                     age.mon + samp_year_cnt, 
#                                   data = luma_data_cub)
      
    ## g) Parameter estimates
#      summary(cub.3.6.prim.prey.sens)  # model parameter estimates
#      confint(cub.3.6.prim.prey.sens)  # 95% CIs 
      
      
  ### 8.10 Cub model: methylation by 6 to 9 months prey density     
    ## a) Unadjusted: methlyation by 6 to 9 months prim.prey density
      cub.6.9.prim.prey.unadj <- glm(methylation ~ prim.prey.6.9, 
                                     data = luma_data_cub)
      
    ## b) Parameter estimates
      summary(cub.6.9.prim.prey.unadj)  # model parameter estimates
      confint(cub.6.9.prim.prey.unadj)  # 95% CIs 
      
      plot(cub.6.9.prim.prey.unadj) # view the residual and QQ plots
      
    ## c) Adjusted: methlyation by 6 to 9 months prim.prey density
      cub.6.9.prim.prey.adj <- glm(methylation ~ prim.prey.6.9 + 
                                    prim.prey.3.6 + prim.prey.birth.3 + 
                                    prim.prey.gest + prim.prey.peri.concpt + 
                                    sex + age.mon, 
                                   data = luma_data_cub)
      
    ## d) Parameter estimates
      summary(cub.6.9.prim.prey.adj)  # model parameter estimates
      confint(cub.6.9.prim.prey.adj)  # 95% CIs 
      
    ## e) Check for heteroskedacity, normality, and outliers    
      plot(cub.6.9.prey.adj) # view the residual and QQ plots
      bptest(cub.6.9.prim.prey.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(cub.6.9.prim.prey.adj, vcov = vcovHC(cub.6.9.prim.prey.adj)) 

    ## f) Sensitivity: methlyation by 6 to 9 months prim.prey density
#      cub.6.9.prim.prey.sens <- glm(methylation ~ prim.prey.6.9 + 
#                                     prim.prey.3.6 + prim.prey.birth.3 + 
#                                     prim.prey.gest + prim.prey.peri.concpt + 
#                                     sex + age.mon + samp_year_cnt, 
#                                   data = luma_data_cub)
      
    ## g) Parameter estimates
#      summary(cub.6.9.prim.prey.sens)  # model parameter estimates
#      confint(cub.6.9.prim.prey.sens)  # 95% CIs 
        
      
  ### 8.11 Graph of cub %CCGG methylation by soc and ecologica factors 
    ## a) make tidy tables of glm model parameter estimates using broom and
      # dotwhisker pckgs for all adjusted cub models
      # terms, including intercept can be dropped from tidy table for graphing
      # purposes, and estimates can be re-labeled 
      rank_est_cub <- tidy(cub.mom.rank.adj) %>%
        filter(term != "(Intercept)" &
               term != "sexm" &
               term != "age.mon") %>%
        relabel_predictors(c(mom.strank.quartQ2 = "Q2 (maternal rank)",                       
                             mom.strank.quartQ3 = "Q3 (maternal rank)", 
                             mom.strank.quartQ4 = "Q4 (maternal rank)")) 
      
      lit_size_est_cub <- tidy(cub.lit.size.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon") %>%
        relabel_predictors(c(lit.sizetwin = "Twin litter")) 
      
      hum_est_cub <- tidy(cub.hum.pres.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon") %>%
        relabel_predictors(c(hum.presmed = "Medium human disturbance",
                             hum.preshi = "High human disturbance")) 
                             
      peri_pre_est_cub <- tidy(cub.peri.prim.prey.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon") %>%
        relabel_predictors(c(prim.prey.peri.concpt = 
                               "Periconceptional prey density")) 
      
      gest_prey_est_cub <- tidy(cub.gest.prim.prey.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon" &
                 term != "prim.prey.peri.concpt") %>%
        relabel_predictors(c(prim.prey.gest = 
                               "Gestational prey density"))
      
      zero_prey_est_cub <- tidy(cub.birth.3.prim.prey.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon" &
                 term != "prim.prey.peri.concpt" &
                 term != "prim.prey.gest") %>%
        relabel_predictors(c(prim.prey.birth.3 = 
                               "Birth - 3 month prey density"))
      
      three_prey_est_cub <- tidy(cub.3.6.prim.prey.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon" &
                 term != "prim.prey.peri.concpt" &
                 term != "prim.prey.gest" &
                 term != "prim.prey.birth.3") %>%
        relabel_predictors(c(prim.prey.3.6 = 
                               "3-6 month prey density"))
      
      six_prey_est_cub <- tidy(cub.6.9.prim.prey.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon" &
                 term != "prim.prey.peri.concpt" &
                 term != "prim.prey.gest" &
                 term != "prim.prey.birth.3" &
                 term != "prim.prey.3.6") %>%
        relabel_predictors(c(prim.prey.6.9 = 
                               "6-9 month prey density"))
      
    ## b) Combine tidy tables of glm estimates into a single tidy table
      soc_ecolog_meth_ests <- bind_rows(rank_est_cub, lit_size_est_cub, 
                                        hum_est_cub, peri_pre_est_cub, 
                                        gest_prey_est_cub, zero_prey_est_cub,
                                        three_prey_est_cub, six_prey_est_cub)
      
    ## c) Graph of the beta estimates from all models 1 which models each
      # predictor with %CCGG methylation in separate models 
      # (control for cub age and sex, as well as previous prey density)
      # uses dotwhisker, broom, dplyr, and ggplot2 packages
      dwplot(soc_ecolog_meth_ests,
             vline = geom_vline(xintercept = 0, colour = "red", 
                                linetype = 2)) + # line at zero behind coefs
        #geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
        ggtitle("Adjusted cub models: Associations of each 
explanatory variable with %CCGG methylation",
                subtitle = '(Each model is controlled for cub sex and age, as well as previous prey density in prey models)')+
        theme(plot.title = element_text( hjust = 0.5)) + # center title
        theme(plot.subtitle = element_text(hjust = 0.5, size = 12)) + 
        # bold and size title and axes labels
        theme(text = element_text(size=18, face = "bold")) +
        theme(legend.position = "none") + # remove legend
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = "white")) +
        # add major axes
        #theme(axis.line = element_line(colour = "lightgrey", 
        #                               size = 1, linetype = "solid")) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face="bold", color="black", 
                                         size=16, angle=0,
                                         margin = margin(t = 10, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face="bold", color="black", 
                                         size=16, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10))) +
        #scale_color_grey (start = 0, end = 0) + # make color estimates black
        scale_color_manual(values=c("black")) + # make color estimates black
        xlab(expression(atop(bold("Coefficient Estimate"), 
                             paste(symbol('\254'),
                                   italic("Difference in %CCGG methylation"), 
                                   symbol('\256'))))) + 
        ylab("") 
      
      ## h) Save Plot
      # use ggsave to save the linearization plot
      ggsave("cub_models1_plot.pdf", plot = last_plot(), device = NULL,
             path = paste0(here(),"/output/output_luma_ecolog"), 
             scale = 1, width = 11.5,
             height = 9,
             units = c("in"), dpi = 300, limitsize = TRUE)

   
  ### 8.12 Cub model: mutual adjustment by significatn predictors 
    ## a) Adjusted: methlyation by mom rank and birth to 3 months prey density
      cub.mutual.adj <- glm(methylation ~ mom.strank.quart + hum.pres + 
                              prim.prey.peri.concpt + prim.prey.gest + 
                              prim.prey.birth.3 + sex + age.mon, 
                                 data = luma_data_cub)
      
    ## b) Parameter estimates
      summary(cub.mutual.adj)  # model parameter estimates
      confint(cub.mutual.adj)  # 95% CIs 
 
    ## c) make a tidy table of glm model parameter estimates using broom pckg
      # terms, including intercept can be dropped from tidy table for graphing
      # purposes, and estimates can be re-labeled
      sig_pred_mut_adg <- tidy(cub.mutual.adj) %>%
        filter(term != "(Intercept)" &
                 term != "prim.prey.gest" &
                 term != "sexm" &
                 term != "age.mon") %>%
        relabel_predictors(c(mom.strank.quartQ2 = "Q2 (maternal rank)",                       
                             mom.strank.quartQ3 = "Q3 (maternal rank)", 
                             mom.strank.quartQ4 = "Q4 (maternal rank)", 
                             hum.presmed = "Medium (anthro. distubance)", 
                             hum.preshi = "High (anthro. distubance)", 
                             prim.prey.peri.concpt = "Periconception Prey",
                             prim.prey.birth.3 = "Birth to 3 month Prey"))
      
    ## d) Graph of the beta estimates from model 2 which includes significant
      # predictors from model now mutally adjustd for each other
      # uses dotwhisker, broom, dplyr, and ggplot2 packages
      dwplot(sig_pred_mut_adg,
             vline = geom_vline(xintercept = 0, colour = "red", 
                                linetype = 2)) + # line at zero behind coefs
        #geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
        ggtitle("Model 2: Mutually adjusted predictors 
of %CCGG methylation in cubs") +
        # bold and size title and axes labels
        theme(text = element_text(size=18, face = "bold")) +
        theme(plot.title = element_text( hjust = 0.5)) + # center title
        theme(legend.position = "none") + # remove legend
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = "white")) +
        # add major axes
        #theme(axis.line = element_line(colour = "lightgrey", 
        #                               size = 1, linetype = "solid")) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=16, angle=0,
                                   margin = margin(t = 10, r = 0, 
                                                   b = 10, l = 0)),
              axis.text.y = element_text(face="bold", color="black", 
                                    size=16, angle=0, 
                                    margin = margin(t = 0, r = 0, 
                                                    b = 0, l = 10))) +
        #scale_color_grey (start = 0, end = 0) + # make color estimates black
        scale_color_manual(values=c("black")) + # make color estimates black
        xlab(expression(atop(bold("Coefficient Estimate"), 
                             paste(symbol('\254'),
                                   italic("Difference in %CCGG methylation"), 
                                   symbol('\256'))))) + 
        ylab("") 

    ## e) Save Plot
      # use ggsave to save the linearization plot
      ggsave("cub_mut_adjust_plot.pdf", plot = last_plot(), device = NULL,
       path = paste0(here(),"/output/output_luma_ecolog"), 
       scale = 1, width = 11,
       height = 8,
       units = c("in"), dpi = 300, limitsize = TRUE)


      
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
      # generate p-value from type II Wald test
      car::Anova(sub.age.mon.adj,Type ="II", test = "Wald") 
      
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
      
      plot(sub.mom.rank.unadj) # view the residual and QQ plots
      
    ## d) Adjusted: methlyation by mom.strank.quart
      sub.mom.rank.adj <- glm(methylation ~ mom.strank.quart + sex + 
                                    age.mon,
                                  data = luma_data_sub)
      
    ## e) Parameter estimates
      summary(sub.mom.rank.adj)  # print model summary, effects and SE
      confint(sub.mom.rank.adj)  # print 95% CIs for parameter estimates
      
    ## f) Check for heteroskedacity, normality, and outliers    
      plot(sub.mom.rank.adj) # view the residual and QQ plots
      bptest(sub.mom.rank.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(sub.mom.rank.adj, vcov = vcovHC(sub.mom.rank.adj))
      
    ## g) Sesnitivity: methlyation by mom.strank.quart
#      sub.mom.rank.sens <- glm(methylation ~ mom.strank.quart + sex + 
#                                age.mon + samp_year_cnt,
#                              data = luma_data_sub)
      
    ## h) Parameter estimates
#      summary(sub.mom.rank.sens)  # print model summary, effects and SE
#      confint(sub.mom.rank.sens)  # print 95% CIs for parameter estimates


    ## i) Extract mom.strank.quart estimates and 
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
      
    ## j) Graph sub.mom.rank effects
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
      
    ## k) Save Plot
      # use ggsave to save the linearization plot
      ggsave("mat_rank_sub_mod_beta.pdf", plot = last_plot(), device = NULL,
             path = paste0(here(),("/output/output_luma_ecolog")),
                           scale = 1, width = 7,
             height = 5.5,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## l) Do a post-hoc test to determine if Q2, Q3, and Q4 differ
      pairwise.t.test(luma_data_sub$methylation,
                      luma_data_sub$mom.strank.quart,
                      p.adj = "none")
      TukeyHSD(aov(methylation ~ mom.strank.quart + sex,
                   data = luma_data_sub))
      
    ## m) Combine quartiles 2-4 into a single category
      luma_data_sub$mom.strank.quart.comb <- as.factor(
        ifelse(luma_data_sub$mom.strank.quart ==  "Q1 (lowest)",
               "Q1 (lowest)", "Q2-Q4 (highest)"))
      
    ## n)  Adjusted: methlyation by mom.strank.quart binned
      sub.rank.adj2 <- glm(methylation ~ mom.strank.quart.comb + sex + 
                             age.mon + samp_year_cnt,
                           data = luma_data_sub)
      
      summary(sub.rank.adj2)  # print model summary, effects and SE
      confint(sub.rank.adj2)  # print 95% CIs for parameter estimates

      
      
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
      
      plot(sub.lit.size.unadj) # view the residual and QQ plots
      
    ## d) Adjusted: methlyation by age.mon
      sub.lit.size.adj <- glm(methylation ~  lit.size + age.mon + sex,
                              data = luma_data_sub)
      
    ## e) Parameter estimates
      summary(sub.lit.size.adj)  # model parameter estimates
      confint(sub.lit.size.adj)  # 95% CIs 

      
    ## f) Check for heteroskedacity, normality, and outliers    
      plot(sub.lit.size.adj) # view the residual and QQ plots
      bptest(sub.lit.size.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(sub.lit.size.adj, vcov = vcovHC(sub.lit.size.adj))
    
    ## g) Sensitivity: methlyation by age.mon
#      sub.lit.size.sens <- glm(methylation ~  lit.size + age.mon + sex +
#                                samp_year_cnt,
#                              data = luma_data_sub)
      
    ## h) Parameter estimates
#      summary(sub.lit.size.sens)  # model parameter estimates
#      confint(sub.lit.size.sens)  # 95% CIs 

        
      
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
      
      plot(sub.hum.pres.unadj) # view the residual and QQ plots
      
    ## d) Adjusted: methlyation by age.mon
      sub.hum.pres.adj <- glm(methylation ~  hum.pres + age.mon + sex,
                              data = luma_data_sub)
      
    ## e) Parameter estimates
      summary(sub.hum.pres.adj)  # model parameter estimates
      confint(sub.hum.pres.adj)  # 95% CIs 
      
    ## f) Check for heteroskedacity, normality, and outliers    
      plot(sub.hum.pres.adj) # view the residual and QQ plots
      bptest(sub.hum.pres.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(sub.hum.pres.adj, vcov = vcovHC(sub.hum.pres.adj))
    
    ## g) Sensitivity: methlyation by age.mon
      sub.hum.pres.sens <- glm(methylation ~  hum.pres + age.mon + sex +
                                samp_year_cnt,
                              data = luma_data_sub)
      
    ## h) Calculate VIF
      vif(sub.hum.pres.sens)
      
    ## i) Parameter estimates
      summary(sub.hum.pres.sens)  # model parameter estimates
      confint(sub.hum.pres.sens)  # 95% CIs 
      
    ## j) Sensitivity: methlyation by samp_year_cnt
      sub.samp.age.sens <- glm(methylation ~ age.mon + sex +
                                 samp_year_cnt,
                               data = luma_data_sub)
      
    ## k) Parameter estimates
      summary(sub.samp.age.sens)  # model parameter estimates
      confint(sub.samp.age.sens)  # 95% CIs 
      
    ## l) Sensitivity: meth by mat_rank, while controling for samp_year_cnt
      sub.rank.samp.age.sens <- glm(methylation ~ mom.strank.quart + age.mon +
                                      sex + samp_year_cnt,
                                    data = luma_data_sub)
      
    ## m) Parameter estimates
      summary(sub.rank.samp.age.sens)  # model parameter estimates
      confint(sub.rank.samp.age.sens)  # 95% CIs 
      
      
  ### 9.6 sub model: methylation by periconceptional prey density     
    ## a) Unadjusted: methlyation by periconceptional prim.prey density
      sub.peri.prim.unadj <- glm(methylation ~ prim.prey.peri.concpt, 
                                     data = luma_data_sub)
      
    ## b) Parameter estimates
      summary(sub.peri.prim.unadj)  # model parameter estimates
      confint(sub.peri.prim.unadj)  # 95% CIs 
      
    ## c) Adjusted: methlyation by periconceptional prim.prey density
      sub.peri.prim.adj <- glm(methylation ~ prim.prey.peri.concpt + sex + 
                                     age.mon, 
                                   data = luma_data_sub)
      
    ## d) Parameter estimates
      summary(sub.peri.prim.adj)  # model parameter estimates
      confint(sub.peri.prim.adj)  # 95% CIs
      
      plot(sub.peri.prim.unadj) # view the residual and QQ plots
      
    ## e) Check for heteroskedacity, normality, and outliers    
      plot(sub.peri.prim.adj) # view the residual and QQ plots
      bptest(sub.peri.prim.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(sub.peri.prim.adj, vcov = vcovHC(sub.peri.prim.adj))
      
    ## f) Sensitivity: methlyation by periconceptional prim.prey density
    #  sub.peri.prim.sens <- glm(methylation ~ prim.prey.peri.concpt + sex + 
    #                             age.mon + samp_year_cnt, 
    #                           data = luma_data_sub)
      
    ## g) Parameter estimates
    #  summary(sub.peri.prim.sens)  # model parameter estimates
    #  confint(sub.peri.prim.sens)  # 95% CIs   
      
      
  ### 9.7 sub model: methylation by gestational prey density     
    ## a) Unadjusted: methlyation by gestational prim.prey density
      sub.gest.prim.unadj <- glm(methylation ~ prim.prey.gest, 
                                     data = luma_data_sub)
      
    ## b) Parameter estimates
      summary(sub.gest.prim.unadj)  # model parameter estimates
      confint(sub.gest.prim.unadj)  # 95% CIs 
      
      plot(sub.gest.prim.unadj) # view the residual and QQ plots
      
    ## c) Adjusted: methlyation by gestational prim.prey density
      sub.gest.prim.adj <- glm(methylation ~ prim.prey.gest + 
                                     prim.prey.peri.concpt + sex + 
                                     age.mon, 
                                   data = luma_data_sub)
      
    ## d) Parameter estimates
      summary(sub.gest.prim.adj)  # model parameter estimates
      confint(sub.gest.prim.adj)  # 95% CIs 
      
    ## e) Check for heteroskedacity, normality, and outliers    
      plot(sub.gest.prim.adj) # view the residual and QQ plots
      bptest(sub.gest.prim.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(sub.gest.prim.adj, vcov = vcovHC(sub.gest.prim.adj))
    
    ## f) Sensitivity: methlyation by gestational prim.prey density
#      sub.gest.prim.sens <- glm(methylation ~ prim.prey.gest + 
#                                 prim.prey.peri.concpt + sex + 
#                                 age.mon + samp_year_cnt, 
#                               data = luma_data_sub)
      
     ## g) Parameter estimates
#      summary(sub.gest.prim.sens)  # model parameter estimates
#      confint(sub.gest.prim.sens)  # 95% CIs  
       
      
  ### 9.8 sub model: methylation by birth to 3 months prey density     
    ## a) Unadjusted: methlyation by birth to 3 months prim.prey density
      sub.birth.3.prim.unadj <- glm(methylation ~ prim.prey.birth.3, 
                                        data = luma_data_sub)
      
    ## b) Parameter estimates
      summary(sub.birth.3.prim.unadj)  # model parameter estimates
      confint(sub.birth.3.prim.unadj)  # 95% CIs 
      
      plot(sub.birth.3.prim.unadj) # view the residual and QQ plots
      
    ## c) Adjusted: methlyation by birth to 3 months prim.prey density
      sub.birth.3.prim.adj <- glm(methylation ~ prim.prey.birth.3 + 
                                        prim.prey.gest + prim.prey.peri.concpt + 
                                        sex + age.mon, 
                                      data = luma_data_sub)
      
    ## d) Parameter estimates
      summary(sub.birth.3.prim.adj)  # model parameter estimates
      confint(sub.birth.3.prim.adj)  # 95% CIs 
      
    ## e) Check for heteroskedacity, normality, and outliers    
      plot(sub.birth.3.prim.adj) # view the residual and QQ plots
      bptest(sub.birth.3.prim.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(sub.birth.3.prim.adj, vcov = vcovHC(sub.birth.3.prim.adj))
      
    ## f) Sensitivity: methlyation by birth to 3 months prim.prey density
#      sub.birth.3.prim.sens <- glm(methylation ~ prim.prey.birth.3 + 
#                                    prim.prey.gest + prim.prey.peri.concpt + 
#                                    sex + age.mon + samp_year_cnt, 
#                                  data = luma_data_sub)
      
    ## g) Parameter estimates
#      summary(sub.birth.3.prim.sens)  # model parameter estimates
#      confint(sub.birth.3.prim.sens)  # 95% CIs 
      
      
  ### 9.9 sub model: methylation by 3 to 6 months prey density     
    ## a) Unadjusted: methlyation by 3 to 6 months prim.prey density
      sub.3.6.prim.unadj <- glm(methylation ~ prim.prey.3.6, 
                                    data = luma_data_sub)
      
    ## b) Parameter estimates
      summary(sub.3.6.prim.unadj)  # model parameter estimates
      confint(sub.3.6.prim.unadj)  # 95% CIs 
      
      plot(sub.3.6.prim.unadj) # view the residual and QQ plots
      
    ## c) Adjusted: methlyation by 3 to 6 months prim.prey density
      sub.3.6.prim.adj <- glm(methylation ~ prim.prey.3.6 + 
                                    prim.prey.birth.3 + prim.prey.gest + 
                                    prim.prey.peri.concpt + sex + 
                                    age.mon, 
                                  data = luma_data_sub)
      
    ## d) Parameter estimates
      summary(sub.3.6.prim.adj)  # model parameter estimates
      confint(sub.3.6.prim.adj)  # 95% CIs 
      
    ## e) Check for heteroskedacity, normality, and outliers    
      plot(sub.3.6.prim.adj) # view the residual and QQ plots
      bptest(sub.3.6.prim.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(sub.3.6.prim.adj, vcov = vcovHC(sub.3.6.prim.adj))
      
    ## f) Sensitivity: methlyation by 3 to 6 months prim.prey density
#      sub.3.6.prim.sens <- glm(methylation ~ prim.prey.3.6 + 
#                                prim.prey.birth.3 + prim.prey.gest + 
#                                prim.prey.peri.concpt + sex + 
#                                age.mon + samp_year_cnt, 
#                              data = luma_data_sub)
      
    ## g) Parameter estimates
#      summary(sub.3.6.prim.sens)  # model parameter estimates
#      confint(sub.3.6.prim.sens)  # 95% CIs   
      
   
  ### 9.10 sub model: methylation by 6 to 9 months prey density     
    ## a) Unadjusted: methlyation by 6 to 9 months prim.prey density
      sub.6.9.prim.unadj <- glm(methylation ~ prim.prey.6.9, 
                                    data = luma_data_sub)
      
    ## b) Parameter estimates
      summary(sub.6.9.prim.unadj)  # model parameter estimates
      confint(sub.6.9.prim.unadj)  # 95% CIs 
      
      plot(sub.6.9.prim.unadj) # view the residual and QQ plots
      
    ## c) Adjusted: methlyation by 6 to 9 months prim.prey density
      sub.6.9.prim.adj <- glm(methylation ~ prim.prey.6.9 + 
                                    prim.prey.3.6 + prim.prey.birth.3 + 
                                    prim.prey.gest + prim.prey.peri.concpt + 
                                    sex + age.mon, 
                                  data = luma_data_sub)
      
    ## d) Parameter estimates
      summary(sub.6.9.prim.adj)  # model parameter estimates
      confint(sub.6.9.prim.adj)  # 95% CIs 
      
    ## e) Check for heteroskedacity, normality, and outliers    
      plot(sub.6.9.prim.adj) # view the residual and QQ plots
      bptest(sub.6.9.prim.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(sub.6.9.prim.adj, vcov = vcovHC(sub.6.9.prim.adj))  
      
    ## f) Sensitivity: methlyation by 6 to 9 months prim.prey density
#      sub.6.9.prim.sens <- glm(methylation ~ prim.prey.6.9 + 
#                                prim.prey.3.6 + prim.prey.birth.3 + 
#                                prim.prey.gest + prim.prey.peri.concpt + 
#                                sex + age.mon + samp_year_cnt, 
#                              data = luma_data_sub)
      
    ## g) Parameter estimates
#      summary(sub.6.9.prim.sens)  # model parameter estimates
#      confint(sub.6.9.prim.sens)  # 95% CIs 
      
      
  ### 9.11 Graph of sub %CCGG methylation by soc and ecologica factors 
    ## a) make tidy tables of glm model parameter estimates using broom and
      # dotwhisker pckgs for all adjusted sub models
      # terms, including intercept can be dropped from tidy table for graphing
      # purposes, and estimates can be re-labeled 
      rank_est_sub <- tidy(sub.mom.rank.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon") %>%
        relabel_predictors(c(mom.strank.quartQ2 = "Q2 (maternal rank)",                       
                             mom.strank.quartQ3 = "Q3 (maternal rank)", 
                             mom.strank.quartQ4 = "Q4 (maternal rank)")) 
      
      lit_size_est_sub <- tidy(sub.lit.size.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon") %>%
        relabel_predictors(c(lit.sizetwin = "Twin litter")) 
      
      hum_est_sub <- tidy(sub.hum.pres.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon") %>%
        relabel_predictors(c(hum.presmed = "Medium human disturbance",
                             hum.preshi = "High human disturbance")) 
      
      peri_pre_est_sub <- tidy(sub.peri.prim.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon") %>%
        relabel_predictors(c(prim.prey.peri.concpt = 
                               "Periconceptional prey density")) 
      
      gest_prey_est_sub <- tidy(sub.gest.prim.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon" &
                 term != "prim.prey.peri.concpt") %>%
        relabel_predictors(c(prim.prey.gest = 
                               "Gestational prey density"))
      
      zero_prey_est_sub <- tidy(sub.birth.3.prim.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon" &
                 term != "prim.prey.peri.concpt" &
                 term != "prim.prey.gest") %>%
        relabel_predictors(c(prim.prey.birth.3 = 
                               "Birth - 3 month prey density"))
      
      three_prey_est_sub <- tidy(sub.3.6.prim.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon" &
                 term != "prim.prey.peri.concpt" &
                 term != "prim.prey.gest" &
                 term != "prim.prey.birth.3") %>%
        relabel_predictors(c(prim.prey.3.6 = 
                               "3-6 month prey density"))
      
      six_prey_est_sub <- tidy(sub.6.9.prim.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon" &
                 term != "prim.prey.peri.concpt" &
                 term != "prim.prey.gest" &
                 term != "prim.prey.birth.3" &
                 term != "prim.prey.3.6") %>%
        relabel_predictors(c(prim.prey.6.9 = 
                               "6-9 month prey density"))
      
    ## b) Combine tidy tables of glm estimates into a single tidy table
      soc_ecolog_meth_ests_sub <- bind_rows(rank_est_sub, lit_size_est_sub, 
                                        hum_est_sub, peri_pre_est_sub, 
                                        gest_prey_est_sub, zero_prey_est_sub,
                                        three_prey_est_sub, six_prey_est_sub)
      
    ## c) Graph of the beta estimates from all models 1 which models each
      # predictor with %CCGG methylation in separate models 
      # (control for sub age and sex, as well as previous prey density)
      # uses dotwhisker, broom, dplyr, and ggplot2 packages
      dwplot(soc_ecolog_meth_ests_sub,
             vline = geom_vline(xintercept = 0, colour = "red", 
                                linetype = 2)) + # line at zero behind coefs
        #geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
        ggtitle("Adjusted subadult models: Associations of each \n explanatory variable with %CCGG methylation",
                subtitle = '(Each model is controlled for subadult sex and age, as well as previous prey density in prey models)')+
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(plot.subtitle = element_text(hjust = 0.5, size = 12)) + 
        # bold and size title and axes labels
        theme(text = element_text(size=18, face = "bold")) +
        theme(legend.position = "none") + # remove legend
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = "white")) +
        # add major axes
        #theme(axis.line = element_line(colour = "lightgrey", 
        #                               size = 1, linetype = "solid")) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face="bold", color="black", 
                                         size=16, angle=0,
                                         margin = margin(t = 10, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face="bold", color="black", 
                                         size=16, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10))) +
        #scale_color_grey (start = 0, end = 0) + # make color estimates black
        scale_color_manual(values=c("black")) + # make color estimates black
        xlab(expression(atop(bold("Coefficient Estimate"), 
                             paste(symbol('\254'),
                                   italic("Difference in %CCGG methylation"), 
                                   symbol('\256'))))) + 
        ylab("") 
      
    ## e) Save Plot
      # use ggsave to save the linearization plot
      ggsave("sub_models1_plot.pdf", plot = last_plot(), device = NULL,
             path = paste0(here(),"/output/output_luma_ecolog"), 
             scale = 1, width = 11.6,
             height = 9,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
      
      
      
###############################################################################
##############                  10. Adult models                 ##############
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
      
      plot(adult.mom.rank.unadj) # view the residual and QQ plots
      
    ## d) Adjusted: methlyation by mom.strank.quart
      adult.mom.rank.adj <- glm(methylation ~ mom.strank.quart + sex + 
                                age.mon,
                              data = luma_data_adult)
      
    ## e) Parameter estimates
      summary(adult.mom.rank.adj)  # print model summary, effects and SE
      confint(adult.mom.rank.adj)  # print 95% CIs for parameter estimates
      
    ## f) Check for heteroskedacity, normality, and outliers    
      plot(adult.mom.rank.adj) # view the residual and QQ plots
      bptest(adult.mom.rank.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(adult.mom.rank.adj, vcov = vcovHC(adult.mom.rank.adj))  
    
    ## g) Adjusted own rank: methlyation by strank.quart
#      adult.own.rank.adj <- glm(methylation ~ strank.quart + 
#                                  age.mon,
#                                data = luma_data_adult)
   
    ## h) Parameter estimates
#      summary(adult.own.rank.adj)  # print model summary, effects and SE
#      confint(adult.own.rank.adj)  # print 95% CIs for parameter estimates
      
    ## i) Adjusted mom and own rank: methlyation by mom.strank.quart own strank
      adult.mom.own.rank.adj <- glm(methylation ~ mom.strank.quart + 
                                      strank.quart + age.mon,
                                data = luma_data_adult)
      
    ## j) Parameter estimates
      summary(adult.mom.own.rank.adj)  # print model summary, effects and SE
      confint(adult.mom.own.rank.adj)  # print 95% CIs for parameter estimates
      Anova(adult.mom.own.rank.adj, Type ="II", test = "Wald") # Wald test p  
      
      vif(adult.mom.own.rank.adj)
      
    ## k) Sensitivity: methlyation by mom.strank.quart
#      adult.mom.rank.sens <- glm(methylation ~ mom.strank.quart + sex + 
#                                  age.mon + samp_year_cnt,
#                                data = luma_data_adult)
      
    ## l) Parameter estimates
#      summary(adult.mom.rank.sens)  # print model summary, effects and SE
#      confint(adult.mom.rank.sens)  # print 95% CIs for parameter estimates
      
    ## m) Extract mom.strank.quart estimates and 
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
      
    ## n) Graph adult.mom.rank effects
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
      
    ## o) Save Plot
      # use ggsave to save the linearization plot
      ggsave("mat_rank_adult_mod_beta.pdf", plot = last_plot(), device = NULL,
             path = paste0(here(),"/output/output_luma_ecolog"), 
             scale = 1, width = 7,
             height = 5.5,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## p) Do a post-hoc test to determine if Q2, Q3, and Q4 differ
      pairwise.t.test(luma_data_adult$methylation,
                      luma_data_adult$mom.strank.quart,
                      p.adj = "none")
      TukeyHSD(aov(methylation ~ mom.strank.quart + sex,
                   data = luma_data_adult))
      
    ## q) Combine quartiles 2-4 into a single category
      luma_data_adult$mom.strank.quart.comb <- as.factor(
        ifelse(luma_data_adult$mom.strank.quart ==  "Q1 (lowest)",
               "Q1 (lowest)", "Q2-Q4 (highest)"))
      
    ## r)  Adjusted: methlyation by mom.strank.quart binned
      adult.rank.adj2 <- glm(methylation ~ mom.strank.quart.comb + sex + 
                             age.mon + samp_year_cnt,
                           data = luma_data_adult)
      
      summary(adult.rank.adj2)  # print model summary, effects and SE
      confint(adult.rank.mod2)  # print 95% CIs for parameter estimates
      
      
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
      
      plot(adult.lit.size.unadj) # view the residual and QQ plots

    ## d) Adjusted: methlyation by age.mon
      adult.lit.size.adj <- glm(methylation ~  lit.size + age.mon + sex,
                              data = luma_data_adult)
      
    ## e) Parameter estimates
      summary(adult.lit.size.adj)  # model parameter estimates
      confint(adult.lit.size.adj)  # 95% CIs 
      
    ## f) Check for heteroskedacity, normality, and outliers    
      plot(adult.lit.size.adj) # view the residual and QQ plots
      bptest(adult.lit.size.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(adult.lit.size.adj, vcov = vcovHC(adult.lit.size.adj))  
      
    ## g) Sensitivity: methlyation by age.mon
#      adult.lit.size.sens <- glm(methylation ~  lit.size + age.mon + sex +
#                                  samp_year_cnt,
#                                data = luma_data_adult)
      
    ## h) Parameter estimates
#      summary(adult.lit.size.sens)  # model parameter estimates
#      confint(adult.lit.size.sens)  # 95% CIs
   
         
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
      
      plot(adult.hum.pres.unadj) # view the residual and QQ plots
      
    ## d) Adjusted: methlyation by age.mon
      adult.hum.pres.adj <- glm(methylation ~  hum.pres + age.mon + sex,
                              data = luma_data_adult)
      
    ## e) Parameter estimates
      summary(adult.hum.pres.adj)  # model parameter estimates
      confint(adult.hum.pres.adj)  # 95% CIs 
      
    ## f) Check for heteroskedacity, normality, and outliers    
      plot(adult.hum.pres.adj) # view the residual and QQ plots
      bptest(adult.hum.pres.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(adult.hum.pres.adj, vcov = vcovHC(adult.hum.pres.adj))  
      
      
    ## g) Sensitivity: methlyation by age.mon
      adult.hum.pres.sens <- glm(methylation ~   age.mon + sex +
                                  samp_year_cnt,
                                data = luma_data_adult)
    
    ## h) Calculate VIF
        vif(adult.hum.pres.sens)
      
    ## i) Parameter estimates
      summary(adult.hum.pres.sens)  # model parameter estimates
      confint(adult.hum.pres.sens)  # 95% CIs 
      
    ## j) Sensitivity: methlyation by samp_year_cnt
      adult.samp.age.sens <- glm(methylation ~ age.mon + sex +
                                 samp_year_cnt,
                               data = luma_data_adult)
      
    ## k) Parameter estimates
      summary(adult.samp.age.sens)  # model parameter estimates
      confint(adult.samp.age.sens)  # 95% CIs 
      
    ## l) Sensitivity: meth by mat_rank, while controling for samp_year_cnt
      adult.rank.samp.age.sens <- glm(methylation ~ mom.strank.quart + age.mon +
                                      sex + samp_year_cnt,
                                    data = luma_data_adult)
      
    ## m) Parameter estimates
      summary(adult.rank.samp.age.sens)  # model parameter estimates
      confint(adult.rank.samp.age.sens)  # 95% CIs 
      
      
  ### 10.6 adult model: methylation by periconceptional prey density     
    ## a) Unadjusted: methlyation by periconceptional prim.prey density
      adult.peri.prim.unadj <- glm(methylation ~ prim.prey.peri.concpt, 
                                     data = luma_data_adult)
      
    ## b) Parameter estimates
      summary(adult.peri.prim.unadj)  # model parameter estimates
      confint(adult.peri.prim.unadj)  # 95% CIs 
      
      plot(adult.peri.prim.unadj) # view the residual and QQ plots

    ## c) Adjusted: methlyation by periconceptional prim.prey density
      adult.peri.prim.adj <- glm(methylation ~ prim.prey.peri.concpt + sex + 
                                     age.mon, 
                                   data = luma_data_adult)

    ## d) Parameter estimates
      summary(adult.peri.prim.adj)  # model parameter estimates
      confint(adult.peri.prim.adj)  # 95% CIs 
   
    ## e) Check for heteroskedacity, normality, and outliers    
      plot(adult.peri.prim.adj) # view the residual and QQ plots
      bptest(adult.peri.prim.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(adult.peri.prim.adj, vcov = vcovHC(adult.peri.prim.adj))  
      
    ## f) Sensitivity: methlyation by periconceptional prim.prey density
#      adult.peri.prim.sens <- glm(methylation ~ prim.prey.peri.concpt + sex + 
#                                   age.mon + samp_year_cnt, 
#                                 data = luma_data_adult)
      
    ## g) Parameter estimates
#      summary(adult.peri.prim.sens)  # model parameter estimates
#      confint(adult.peri.prim.sens)  # 95% CIs 
      
        
  ### 10.7 adult model: methylation by gestational prey density     
    ## a) Unadjusted: methlyation by gestational prim.prey density
      adult.gest.prim.unadj <- glm(methylation ~ prim.prey.gest, 
                                     data = luma_data_adult)
      
    ## b) Parameter estimates
      summary(adult.gest.prim.unadj)  # model parameter estimates
      confint(adult.gest.prim.unadj)  # 95% CIs 
      
      plot(adult.gest.prim.unadj) # view the residual and QQ plots
      
    ## c) Adjusted: methlyation by gestational prim.prey density
      adult.gest.prim.adj <- glm(methylation ~ prim.prey.gest + 
                                     prim.prey.peri.concpt + sex + 
                                     age.mon, 
                                   data = luma_data_adult)
      
    ## d) Parameter estimates
      summary(adult.gest.prim.adj)  # model parameter estimates
      confint(adult.gest.prim.adj)  # 95% CIs 
      
    ## e) Check for heteroskedacity, normality, and outliers    
      plot(adult.gest.prim.adj) # view the residual and QQ plots
      bptest(adult.gest.prim.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(adult.gest.prim.adj, vcov = vcovHC(adult.gest.prim.adj))
      
    ## f) Sensitivity: methlyation by gestational prim.prey density
#      adult.gest.prim.sens <- glm(methylation ~ prim.prey.gest + 
#                                   prim.prey.peri.concpt + sex + 
#                                   age.mon + samp_year_cnt, 
#                                 data = luma_data_adult)
      
    ## g) Parameter estimates
#      summary(adult.gest.prim.sens)  # model parameter estimates
#      confint(adult.gest.prim.sens)  # 95% CIs   
      
      
  ### 10.8 adult model: methylation by birth to 3 months prey density     
    ## a) Unadjusted: methlyation by birth to 3 months prim.prey density
      adult.birth.3.prim.unadj <- glm(methylation ~ prim.prey.birth.3, 
                                        data = luma_data_adult)
      
    ## b) Parameter estimates
      summary(adult.birth.3.prim.unadj)  # model parameter estimates
      confint(adult.birth.3.prim.unadj)  # 95% CIs 
      
      plot(adult.birth.3.prim.unadj) # view the residual and QQ plots
      
    ## c) Adjusted: methlyation by birth to 3 months prim.prey density
      adult.birth.3.prim.adj <- glm(methylation ~ prim.prey.birth.3 + 
                                        prim.prey.gest + prim.prey.peri.concpt + 
                                        sex + age.mon, 
                                      data = luma_data_adult)
      
    ## d) Parameter estimates
      summary(adult.birth.3.prim.adj)  # model parameter estimates
      confint(adult.birth.3.prim.adj)  # 95% CIs 
      
    ## e) Check for heteroskedacity, normality, and outliers    
      plot(adult.birth.3.prim.adj) # view the residual and QQ plots
      bptest(adult.birth.3.prim.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(adult.birth.3.prim.adj, vcov = vcovHC(adult.birth.3.prim.adj)) 
      
    ## f) Sensitivity: methlyation by birth to 3 months prim.prey density
#      adult.birth.3.prim.sens <- glm(methylation ~ prim.prey.birth.3 + 
#                                      prim.prey.gest + prim.prey.peri.concpt + 
#                                      sex + age.mon + samp_year_cnt, 
#                                    data = luma_data_adult)
      
    ## g) Parameter estimates
#      summary(adult.birth.3.prim.sens)  # model parameter estimates
#      confint(adult.birth.3.prim.sens)  # 95% CIs   
      
      
  ### 10.9 adult model: methylation by 3 to 6 months prey density     
    ## a) Unadjusted: methlyation by 3 to 6 months prim.prey density
      adult.3.6.prim.unadj <- glm(methylation ~ prim.prey.3.6, 
                                    data = luma_data_adult)
      
    ## b) Parameter estimates
      summary(adult.3.6.prim.unadj)  # model parameter estimates
      confint(adult.3.6.prim.unadj)  # 95% CIs 
      
      plot(adult.3.6.prim.unadj) # view the residual and QQ plots
      
    ## c) Adjusted: methlyation by 3 to 6 months prim.prey density
      adult.3.6.prim.adj <- glm(methylation ~ prim.prey.3.6 + 
                                    prim.prey.birth.3 + prim.prey.gest + 
                                    prim.prey.peri.concpt + sex + 
                                    age.mon, 
                                  data = luma_data_adult)
      
    ## d) Parameter estimates
      summary(adult.3.6.prim.adj)  # model parameter estimates
      confint(adult.3.6.prim.adj)  # 95% CIs 
      
    ## e) Check for heteroskedacity, normality, and outliers    
      plot(adult.3.6.prim.adj) # view the residual and QQ plots
      bptest(adult.3.6.prim.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(adult.3.6.prim.adj, vcov = vcovHC(adult.3.6.prim.adj))
      
    ## f) Sensitivity: methlyation by 3 to 6 months prim.prey density
#      adult.3.6.prim.sens <- glm(methylation ~ prim.prey.3.6 + 
#                                  prim.prey.birth.3 + prim.prey.gest + 
#                                  prim.prey.peri.concpt + sex + 
#                                  age.mon + samp_year_cnt, 
#                                data = luma_data_adult)
      
    ## g) Parameter estimates
#      summary(adult.3.6.prim.sens)  # model parameter estimates
#      confint(adult.3.6.prim.sens)  # 95% CIs   
      
      
  ### 10.10 adult model: methylation by 6 to 9 months prey density     
    ## a) Unadjusted: methlyation by 6 to 9 months prim.prey density
      adult.6.9.prim.unadj <- glm(methylation ~ prim.prey.6.9, 
                                    data = luma_data_adult)
      
    ## b) Parameter estimates
      summary(adult.6.9.prim.unadj)  # model parameter estimates
      confint(adult.6.9.prim.unadj)  # 95% CIs 
      
      plot(adult.6.9.prim.unadj) # view the residual and QQ plots
      
    ## c) Adjusted: methlyation by 6 to 9 months prim.prey density
      adult.6.9.prim.adj <- glm(methylation ~ prim.prey.6.9 + 
                                    prim.prey.3.6 + prim.prey.birth.3 + 
                                    prim.prey.gest + prim.prey.peri.concpt + 
                                    sex + age.mon, 
                                  data = luma_data_adult)
      
    ## d) Parameter estimates
      summary(adult.6.9.prim.adj)  # model parameter estimates
      confint(adult.6.9.prim.adj)  # 95% CIs 
      
    ## e) Check for heteroskedacity, normality, and outliers    
      plot(adult.6.9.prim.adj) # view the residual and QQ plots
      bptest(adult.6.9.prim.adj) # heteroskedasticity using the 
      # Breusch-Pagan test; if significant, then generate test
      # Robust Standard Errors
      # Robust Standard Errors (HC3 method)
      coeftest(adult.6.9.prim.adj, vcov = vcovHC(adult.6.9.prim.adj))  
      
    ## f) Sensitivity: methlyation by 6 to 9 months prim.prey density
#      adult.6.9.prim.sens <- glm(methylation ~ prim.prey.6.9 + 
#                                  prim.prey.3.6 + prim.prey.birth.3 + 
#                                  prim.prey.gest + prim.prey.peri.concpt + 
#                                  sex + age.mon + samp_year_cnt, 
#                                data = luma_data_adult)
      
    ## g) Parameter estimates
#      summary(adult.6.9.prim.sens)  # model parameter estimates
#      confint(adult.6.9.prim.sens)  # 95% CIs 
      
      
      
  ### 10.11 Graph of adult %CCGG methylation by soc and ecologica factors 
    ## a) make tidy tables of glm model parameter estimates using broom and
      # dotwhisker pckgs for all adjusted adult models
      # terms, including intercept can be dropped from tidy table for graphing
      # purposes, and estimates can be re-labeled 
      rank_est_adult <- tidy(adult.mom.rank.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon") %>%
        relabel_predictors(c(mom.strank.quartQ2 = "Q2 (maternal rank)",                       
                             mom.strank.quartQ3 = "Q3 (maternal rank)", 
                             mom.strank.quartQ4 = "Q4 (maternal rank)")) 
      
      lit_size_est_adult <- tidy(adult.lit.size.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon") %>%
        relabel_predictors(c(lit.sizetwin = "Twin litter")) 
      
      hum_est_adult <- tidy(adult.hum.pres.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon") %>%
        relabel_predictors(c(hum.presmed = "Medium human disturbance",
                             hum.preshi = "High human disturbance")) 
      
      peri_pre_est_adult <- tidy(adult.peri.prim.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon") %>%
        relabel_predictors(c(prim.prey.peri.concpt = 
                               "Periconceptional prey density")) 
      
      gest_prey_est_adult <- tidy(adult.gest.prim.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon" &
                 term != "prim.prey.peri.concpt") %>%
        relabel_predictors(c(prim.prey.gest = 
                               "Gestational prey density"))
      
      zero_prey_est_adult <- tidy(adult.birth.3.prim.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon" &
                 term != "prim.prey.peri.concpt" &
                 term != "prim.prey.gest") %>%
        relabel_predictors(c(prim.prey.birth.3 = 
                               "Birth - 3 month prey density"))
      
      three_prey_est_adult <- tidy(adult.3.6.prim.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon" &
                 term != "prim.prey.peri.concpt" &
                 term != "prim.prey.gest" &
                 term != "prim.prey.birth.3") %>%
        relabel_predictors(c(prim.prey.3.6 = 
                               "3-6 month prey density"))
      
      six_prey_est_adult <- tidy(adult.6.9.prim.adj) %>%
        filter(term != "(Intercept)" &
                 term != "sexm" &
                 term != "age.mon" &
                 term != "prim.prey.peri.concpt" &
                 term != "prim.prey.gest" &
                 term != "prim.prey.birth.3" &
                 term != "prim.prey.3.6") %>%
        relabel_predictors(c(prim.prey.6.9 = 
                               "6-9 month prey density"))
      
    ## b) Combine tidy tables of glm estimates into a single tidy table
      soc_ecolog_meth_ests_adult <- bind_rows(rank_est_adult, 
                                            lit_size_est_adult, 
                                            hum_est_adult, peri_pre_est_adult, 
                                            gest_prey_est_adult, 
                                            zero_prey_est_adult,
                                            three_prey_est_adult, 
                                            six_prey_est_adult)
      
    ## c) Graph of the beta estimates from all models 1 which models each
      # predictor with %CCGG methylation in separate models 
      # (control for adult age and sex, as well as previous prey density)
      # uses dotwhisker, broom, dplyr, and ggplot2 packages
      dwplot(soc_ecolog_meth_ests_adult,
             vline = geom_vline(xintercept = 0, colour = "red", 
                                linetype = 2)) + # line at zero behind coefs
        #geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
        ggtitle("Adjusted adult models: Associations of each \n explanatory variable with %CCGG methylation",
                subtitle = '(Each model is controlled for adult sex and age, as well as previous prey density in prey models)')+
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(plot.subtitle = element_text(hjust = 0.5, size = 12)) + 
        # bold and size title and axes labels
        theme(text = element_text(size=18, face = "bold")) +
        theme(legend.position = "none") + # remove legend
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = "white")) +
        # add major axes
        #theme(axis.line = element_line(colour = "lightgrey", 
        #                               size = 1, linetype = "solid")) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face="bold", color="black", 
                                         size=16, angle=0,
                                         margin = margin(t = 10, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face="bold", color="black", 
                                         size=16, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10))) +
        #scale_color_grey (start = 0, end = 0) + # make color estimates black
        scale_color_manual(values=c("black")) + # make color estimates black
        xlab(expression(atop(bold("Coefficient Estimate"), 
                             paste(symbol('\254'),
                                   italic("Difference in %CCGG methylation"), 
                                   symbol('\256'))))) + 
        ylab("") 
      
    ## e) Save Plot
      # use ggsave to save the linearization plot
      ggsave("adult_models1_plot.pdf", plot = last_plot(), device = NULL,
             path = paste0(here(),"/output/output_luma_ecolog"), 
             scale = 1, width = 11.6,
             height = 9,
             units = c("in"), dpi = 300, limitsize = TRUE)
      

      
###############################################################################
##############            11. Save data tables as .csv           ##############
###############################################################################         

  # Save intermediate tables as spreadsheets with a .cvs extension and today's
  # date. Files are saved in the 'data' folder or the 'output' folder
  # in the working directory.
  
  
  ### 11.1 Set up date parameters
    # print today's date
    today <- Sys.Date()
    date <- format(today, format="%d%b%Y")
    
  
  ### 11.2 Generate File Names
  # For each table that will be saved as a .csv file, first generate a file 
  # name to save each table
    
    ## a) File name for luma_data table used in analysis of manuscript
      csv.file.name.luma <- paste("~/R/R_wd/fisi/project/", 
                                  "3_hy_GR_global_DNA_meth/",
                                  "LUMA/soc_eco_detrmnts_ms/",
                                  "luma_data",".csv", sep= "")   
    
    ## b) File name for luma_data_group table used in analysis of manuscript
      csv.file.name.luma_data_group <- paste("~/R/R_wd/fisi/project/", 
                                             "3_hy_GR_global_DNA_meth/",
                                              "LUMA/soc_eco_detrmnts_ms/",
                                              "luma_data_group",".csv", 
                                             sep= "")   
    
    
  ### 11.3 Save Tables 
    # Save each data frame as a .csv file (a spreadsheet/table) into the 
    # data folder in the working directory.
    
    ## a) Save luma_data_no_out table
    write.csv (luma_data, file = csv.file.name.luma)
    
    ## b) Save re_runs table
    write.csv (luma_data_group, file = csv.file.name.luma_data_group)
       