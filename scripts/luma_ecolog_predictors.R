###########################################################
########   Spotted Hyena Global DNA Methylation    ########
########         LUMA Ecological Predictors        ########
########             By: Zach Laubach              ########
########        last updated: 12 Sept 2017         ########
###########################################################

### PURPOSE: This code is desingned to analyze LUMA data.

### Special thanks to all contributors and commentors on stackoverflow 
### and other open source communities for creative and informative coding
### solutions.


  # Code Blocks
    # 1: Configure Workspace
    # 2: Set Working Directory
    # 3: Import Data
    # 4: Tidy Data
    # 5: Variable Formatting 
    # 6: Univariate Analysis (Predictors)
    # 7: Bi-Variate Analysis
    # 8: 


###########################################################
####            1  Configure Workspace                 ####
###########################################################

  ### 1.1 clear global environment
    rm(list = ls())


  ### 1.2 Require - load a packages  into a current R session
    ## a) Data Manipulation and Descriptive Stats Packages
      
      require(plyr)
      require(dplyr)
      #require(purrr)
      require(readr)
      #require(reshape)
      require(reshape2)
      require(tidyr)
      require(lubridate)
      options(gsubfn.engine = "R") #fixes tcltk bug; run before require sqldf
      require(sqldf)
    
    ## b) Graph Plotting and Visualization Packages
      require(ggplot2)
      #require(sjPlot)
      #require(sciplot)
      #require(effects)
      #require(plotrix)
      #require(scatterplot3d)
      #require(texi2dvi)
      #require(stargazer)
      #require(xtable)
      require(gridExtra)
  
      
    ## c) Modeling Packages
      # Modeling Packages
      require(lme4)
      require(arm)
      require(car)
      #require(lsmeans)
      #require(boot)
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
    access_data_path <- paste("~/R/R_wd/fisi/access_fisi_data/",
                              sep = '')
    
  
  ### 2.4 Create path to sample selection data folder (from repository query)
    prey_data_path <- paste("~/Git/fisi_lab/hy_prey_density/",
                              sep = '')
    
    
  ### 2.5 Create path to sample selection data folder (from repository query)
    fas_data_path <- paste("/Volumes/Holekamp/CurrentCollaborations/",
                           "maternal_care/R_code/fas/data/",
                           "final_data_with_variables/",
                              sep = '')
    
    
       
###########################################################
####                  3  Import Data                   ####
###########################################################      
  
  ### 3.1 Import LUMA data file
    ## a) Import LUMA data, which has undergone QAQC (see R script 
      # luma_prep_analysis.R)
      luma_data <- read_csv(paste(luma_data_path,
                                    "luma_data_no_out.csv", sep = ''))
  
    
  ### 3.2 Import prey density file
    ## a) Import prey density data which was created with the R script
      # 'calc_prey_density'
      prey_density <- read_csv(paste(prey_data_path, "output/",
                                  "talek_prey_density31Aug2017.csv", sep = ''))  
         
    
  ### 3.3 Import FAS behavior file
    ## a) Import FAS data from L drive
      fas_data <- read_csv(paste(fas_data_path,
                                     "1988-1996FAS_25June2017.csv", sep = ''))  
     
       
  ### 3.4 Import tblFemalerank file
    ## a) Import tblFemalerank from copy of Access backend
      female_rank_data <- read_csv(paste(access_data_path,
                                 "tblFemalerank_jun2017.csv", sep = ''))
  
    
  ### 3.5 Import tblHyenas behavior file
    ## a) Import tblHyenas from copy of Access backend
      tblHyena_data <- read_csv(paste(access_data_path,
                                         "tblHyenas_jun2017.csv", sep = ''))     
      
      
                
###########################################################
####                  4  Tidy Data                     ####
###########################################################   
 
  ### 4.1 Tidy prey_density
    ## a) Select a subset of the prey_density data by column names
      prey_density <- prey_density %>% 
        dplyr:: select(grep("tot", names(prey_density)), # contains 'tot'
                       grep("num", names(prey_density)), # contains 'num'
                       grep("^ID$", names(prey_density))) # exact match 'ID'
      
      
  ### 4.2 Tidy luma_data
    ## a) Filter all rows from LUMA data by hyena's categorical age
     #luma_data <- luma_data %>%
    #   dplyr::filter(grepl('^adult$', Age))
     
    ## b) Select Variables
     # Column names to retain in a list
     luma_data_var_names <- c("meth_adjust", "sample_ID", "Hyena", "KayCode", 
                              "DartingDate", "Sex", "Age", "AgeMonths", 
                              "EstimatedAgeMo", "Clan")
     # Select columns based names from a list  
     luma_data <- (luma_data[,c(luma_data_var_names)]) 
     
    ## c) Rename 'Hyena' variable as 'ID
     luma_data <- rename_(luma_data, "ID" = "Hyena")
  
     
  ### 4.3 Tidy fas_data
    ## a) Filter all rows from FAS data where data are 'ready'
       fas_data <- fas_data %>%
         dplyr::filter(grepl('^ready$', notes)) 
    
    ## b) Calculate Behavioral rates based on time together 
      # (behavior/time together)   
      # nursing rate
       fas_data$n.rate <- (fas_data$n / fas_data$together)
      # close proximity rate
       fas_data$c.rate <- (fas_data$c / fas_data$together)
    
    ## c) Create a categorical age variable
       fas_data$cat.age <- ifelse(fas_data$cubAgeMonths <= 3 , "0-3", "old")
       fas_data$cat.age <- ifelse(fas_data$cubAgeMonths > 3 & 
                                    fas_data$cubAgeMonths <= 6, "3-6", 
                                  fas_data$cat.age)
       fas_data$cat.age <- ifelse(fas_data$cubAgeMonths > 6 & 
                                    fas_data$cubAgeMonths <= 9, "6-9", 
                                  fas_data$cat.age)
       
    ## d) Select Variables
      # Column names to retain in a list
       fas_data_var_names <- c("focal", "date", "cub", "mom", "sex", "cubDOB",
                               "cubAgeMonths", "momrank", "numlitmates",
                               "n.rate", "c.rate", "cat.age")
      # Select columns based names from a list  
       fas_data <- (fas_data[,c(fas_data_var_names)])                      
    
    ## e) Filter rows that are in the appropriate age range
       fas_data <- fas_data %>%
         dplyr::filter(!grepl('^old$', cat.age)) %>%
         rename_("ID" = "cub")  # rename cub variable as ID 
       
    ## f) Summarize the FAS data
       # for repeated measurements within a categorical age, average the 
       # behavior rates for each cub
       fas_data <- plyr::ddply (fas_data, .(ID, cat.age), summarize, 
                          n.rate.mean = mean(n.rate, na.rm = T),
                          c.rate.mean = mean(c.rate, na.rm = T),
                          mom = mom[1], # retain the first element for mom ID
                          cubDOB = cubDOB[1],
                          momrank = momrank [1], # need a single value for
                                                     # mom rank for when FAS 
                                                     # conducted in two years
                          numlitmates = numlitmates[1])
    
    ## g) Reshape FAS data
       fas_data <- fas_data %>%
         gather(key = "var.name", value = "value", c(3:4, 7)) %>% 
            # gather (long format) behavior rate by behavior type and  
            # standardized mom rank by mom rank as two new columns
            # 'var.name' and 'value'
         unite(var.name.age, var.name, cat.age, sep = ".") %>%
            # update the var.name to include the categorical age  
         spread(var.name.age,value)
            # spread to wide format so there is 1 uniqe row for each hyena 
            # with new variables corresponding to var.name.age 
       
  
  ### 4.4 Tidy female_rank_data
    ## a) Make Hyena ID lowercase
       female_rank_data$ID <- sapply(female_rank_data$ID, tolower)
       
    ## b) Pattern recognize numbers from Year variable and copy; gets rid of
      # unwanted text characters
       female_rank_data$year <- as.numeric(regmatches(female_rank_data$Year, 
                                           gregexpr("[[:digit:]]+",
                                                    female_rank_data$Year)))
    ## c) rename 'ID' variable as 'mom' 
      female_rank_data <- rename_(female_rank_data, "mom" = "ID")  
      
    ## d) rename 'Year' variable as 'orig.year' 
      female_rank_data <- rename_(female_rank_data, "orig.year" = "Year")
       
    ## c) convert cub first seen date
      # female_rank_data$Year <- as.POSIXct(as.character 
       #                                    (female_rank_data$Year),
       #                                    format='%Y')
  
            
  ### 4.5 Tidy tbhHyena_data
    ## a) convert cub first seen date
      tblHyena_data$FirstSeen <- as.POSIXct(as.character 
                                        (tblHyena_data$FirstSeen),
                                        format='%m/%d/%y')
    
    ## b) convert cub Disappeared date
      tblHyena_data$Disappeared <- as.POSIXct(as.character 
                                            (tblHyena_data$Disappeared),
                                            format='%m/%d/%y')    
       
    ## c) convert cub Birth date
      tblHyena_data$Birthdate <- as.POSIXct(as.character 
                                            (tblHyena_data$Birthdate),
                                            format='%m/%d/%y')
       
    ## d) convert cub Death date
      tblHyena_data$DeathDate <- as.POSIXct(as.character 
                                              (tblHyena_data$DeathDate),
                                              format='%m/%d/%y')   
      
    ## e) convert cub Weaned date
      tblHyena_data$Weaned <- as.POSIXct(as.character 
                                            (tblHyena_data$Weaned),
                                            format='%m/%d/%y')  
     
                 
  ### 4.6 Join Data Sets: Append covariates to LUMA data
    ## a) Left join prey_density to luma_data   
      luma_data <- sqldf("SELECT
                        luma_data.*           
                        , tblHyena_data. Status, FirstSeen, DenGrad, Disappeared,
                          Mom, Birthdate, NumberLittermates, Litrank, 
                          Fate, DeathDate, Weaned, park
                        FROM luma_data      
                        LEFT JOIN tblHyena_data       
                        ON tblHyena_data.ID = luma_data.ID")     
    
    ## b) extract the year for date of interest (here Birthdate) using lubridate
      # and make a new variable
      luma_data$year <- year(luma_data$Birthdate)
      
    ## c) Left join female_rank_data to luma_data   
      luma_data <- sqldf("SELECT
                        luma_data.*           
                        , female_rank_data. absrank, strank
                        FROM luma_data      
                        LEFT JOIN female_rank_data       
                        ON female_rank_data.mom = luma_data.Mom
                            AND female_rank_data.year = luma_data.year") 
       
    ## d) Left join prey_density to luma_data   
      luma_data<- left_join(luma_data,
                            prey_density, by = "ID") 
      
    ## e) Left join FAS behavior data to luma_data 
#      luma_data_fas <- left_join(luma_data, 
#                             fas_data, by = "ID") 
      
      
      
      
###########################################################
####              5  Variable Formatting               ####
###########################################################  
  
  ### 5.1 Maternal Rank Variable
    ## a) Create Quartiles of Maternal Rank
      # convert strank to numeric
      luma_data$strank <- as.numeric(luma_data$strank)
      # use quantile function to cut strank into 4 levels (approximately same
      # number in each level)
      luma_data <- within(luma_data, strank.quart.order
                          <- as.integer(cut(strank, 
                                         quantile(strank, probs=0:4/4,
                                                  na.rm = T), 
                                         include.lowest = T)))
    
    ## b) Rename Quartile Levels
      #rename maternal rank variable from 1-4 to high, mid, low, bottom as a new 
      # variable. The strank.quart.order variable is retained to define 
      # reference and set group order for graphing
      luma_data$strank.quart <- as.factor (cut(luma_data$strank.quart.order, 
                                                   breaks=c(0, 1, 2, 3, 4), 
                                                   labels=c("high", "mid", 
                                                            "low", "bottom")))
  
  
  ### 5.2 Age Variable
    ## a) Re-order Age Variable
      # change Age variable to a factor from a character
      luma_data$Age <- as.factor(luma_data$Age)
      # re-order the age variable based so it is not base on not alphabetic 
      # order this sets the reference level for ANOVA to adult
      luma_data <- transform(luma_data, 
                             Age = factor(Age,
                                          levels = c("cub", 
                                                     "subadult", "adult"), 
                                          ordered = TRUE))

                                
  ### 5.3  Reduce Data Set
    ## a) Select variables to drop
      luma_data_group <- luma_data %>%
        dplyr::select(- c(sample_ID))    
      
    ## b) Cobmine Age columns
      # replace NA in AgeMonths column where there is a value in 
      # EstimatedAgeMo column
        luma_data_group$AgeMonths <- ifelse(is.na(luma_data_group$AgeMonths),
                                            luma_data_group$EstimatedAgeMo,
                                            luma_data_group$AgeMonths) 
      # Drop extra EstimatedAgeMo column
        luma_data_group <- luma_data_group %>%
          select (- c(EstimatedAgeMo)) 
      
    ## c) Group rows with same ID and within an age group
        luma_data_group <- luma_data_group %>% 
          filter (!is.na(meth_adjust))%>% # remove rows where meth_adjust NA
          group_by (ID, Age) %>%    # set grouping same ID within same cat age
          group_by (meth_adjust = mean(meth_adjust)) %>%  # take avg meth_adjust
                                                  # value for repeat ID within
                                                  # same categorical age
          group_by (AgeMonths = mean(AgeMonths)) %>% # Avg age in months
          mutate (n = n()) %>%       # make var of n for uniqur rows            
          distinct (meth_adjust, .keep_all=TRUE) # keep all rows where var value
                                                # is unique while retaining all
                                                # all other variables
        
      # OR long form listing specific variables to retain
      #  luma_data_group <- luma_data_group %>% 
      #    group_by (ID, Age) %>%
      #    summarise ( meth_adjust = mean(meth_adjust),
      #                KayCode = first(KayCode),
      #                DartingDate = first(DartingDate),
      #                Sex = first(Sex),
      #                AgeMonths = mean(AgeMonths),
      #                Clan = first(Clan),
      #                Status = first(Status),
      #                FirstSeen = first(FirstSeen),
      #                DenGrad = first(DenGrad),
      #                Disappeared = first(Disappeared), 
      #                Mom = first(Mom),
      #                Birthdate = first(Birthdate))
        
      
      
###########################################################
####       6  UniVariate Analysis (Predictors)         ####
###########################################################          
      
  ### 6.1 Univariate Stats of Predictors
    ## a) Predictive variable summary
      pred_var_summry <- ddply (luma_data_group, .(), summarise,
                        mean.age = round (mean(AgeMonths, na.rm = T),
                                          2),
                        fem.sex.raito = (sum(grepl("f",
                                               luma_data_group$Sex))) /
                                    (sum(grepl("m",
                                               luma_data_group$Sex) |
                                           grepl("f",
                                                 luma_data_group$Sex))),
                        mean.peri.prey = round (mean(total.peri_concpt, 
                                                     na.rm = T), 2),
                        mean.gest.prey = round (mean(total.gest, na.rm = T),
                                                2),
                        mean.birth.3.prey = round (mean(total.birth.3, 
                                                        na.rm = T), 2),
                        mean.3.6.prey = round (mean(total.3.6, na.rm = T),
                                               2),
                        mean.6.9.prey = round (mean(total.6.9, na.rm = T),
                                               2))
                        
    
  ### 6.3 Center Predictive Variables
    ## a) center function based on column means    
      center_fxn <- function(x) {
        xcenter = colMeans(x, na.rm = T)
        x - rep(xcenter, rep.int(nrow(x), ncol(x)))
      }
      
    ## b) make a list of variable names to center, value = T is necessary
      # or column positions will be returned
      vars_to_center <- grep("total", names(luma_data_group), value = T)
      
      ## c) Center data in the select columns and replace non  
      luma_data_group[ ,c(vars_to_center)] <- center_fxn(luma_data_group
                                                   [ ,c(vars_to_center)])
      

  

###########################################################
####              7  Bi-Variate Analysis               ####
###########################################################          
      
  ### 7.1 cross tabulation of predicitive variables to get n per group
      table(luma_data_group$Sex, luma_data_group$strank.quart)
      table(luma_data_group$Age, luma_data_group$strank.quart)
      table(luma_data_group$Age, luma_data_group$Sex)
      table(luma_data_group$Age, luma_data_group$strank.quart, 
            luma_data_group$Sex) 
      
      
  ### 7.2 Bivariate Statistics Methylation by Sex
    ## a) Summary Stats Methylation by Sex
      meth_by_sex = ddply(luma_data_group, .(Sex), summarise,
                             N = sum(!is.na(meth_adjust)),
                             mean = round (mean(meth_adjust, na.rm = T), 2),
                             median =  round (quantile(meth_adjust, c(.5), 
                                                       na.rm = T), 2),
                             sd = round (sd(meth_adjust, na.rm = T), 2))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf("output/meth_by_sex.pdf", height = 6, width = 7)
      grid.table(meth_by_sex)
      dev.off()
      
    ## c) Plot Methylation by Sex
      ggplot(data = subset(luma_data_group, !is.na(x = Sex)), 
             aes(x = Sex, y = meth_adjust, color = Sex)) + 
        geom_point(shape = 1) +
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        labs(title = "Percent Global DNA 
Mehtylation by Sex") +
        ylab("% Global DNA Methylation") +
        xlab("Sex")
      
    ## d) boxplot of DNA methylation by Sex
      boxplot(meth_adjust ~ Sex,data = luma_data_group, 
              main = "Boxplots of Hyena Global DNA 
              Methylation by Sex",
              xlab = "Sex", 
              ylab = "% Global DNA Methylation", horizontal=F, 
              col=grey.colors(2))
      
    ## e) Bivariate Regression Methylatin by Sex
      sex.mod <- lmer(meth_adjust ~ Sex + (1|ID), 
                      data = luma_data_group, na.action = exclude)
 
      summary(sex.mod)  # print model summary, effects and SE
      confint(sex.mod)  # print 95% CIs for parameter estimates 
     
      
  ### 7.3 Bivariate Statistics Methylation by Age 
    ## a) Summary Stats Methylation by Age
      meth_by_age = ddply(luma_data_group, .(Age), summarise,
                          N = sum(!is.na(meth_adjust)),
                          mean = round (mean(meth_adjust, na.rm = T), 2),
                          median =  round (quantile(meth_adjust, c(.5), 
                                                    na.rm = T), 2),
                          sd = round (sd(meth_adjust, na.rm = T), 2))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf("output/meth_by_age.pdf", height = 6, width = 7)
      grid.table(meth_by_age)
      dev.off()
      
    ## c) Plot Methylation by Age
      # graph of the raw data for percent global DNA methylaiton by maternal rank
      # stratified by age
      ggplot(data = subset(luma_data_group, !is.na(x = Age)), 
             aes(x = Age, y = meth_adjust, color = Age)) + 
        geom_point(shape = 1) +
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        labs(title = "Percent Global DNA 
Mehtylation by Age") +
        ylab("% Global DNA Methylation") +
        xlab("Categorical Age")
      
    ## d) Boxplot of DNA methylation across factors of age
      boxplot(meth_adjust ~ Age,data = luma_data_group, 
              main = "Boxplots of Hyena Global DNA 
              Methylation by Age",
              xlab = "Age", 
              ylab = "% Global DNA Methylation", horizontal=F, 
              col=grey.colors(3))
      
    ## e) Bivariate Regression Methylatino by Age  
      age.mod <- lmer(meth_adjust ~ Age + (1|ID), 
                      data = luma_data_group)
    
      summary(age.mod)  # print model summary, effects and SE
      confint(age.mod)  # print 95% CIs for parameter estimates
      
      
  ### 7.4 Bivariate Statistics Methylation by Rank
    ## a) Summary Stats Methylation by Rank
      meth_by_rank = ddply(luma_data_group, .(strank.quart), summarise,
                          N = sum(!is.na(meth_adjust)),
                          mean = round (mean(meth_adjust, na.rm = T), 2),
                          median =  round (quantile(meth_adjust, c(.5), 
                                                    na.rm = T), 2),
                          sd = round (sd(meth_adjust, na.rm = T), 2))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf("output/meth_by_rank.pdf", height = 6, width = 7)
      grid.table(meth_by_rank)
      dev.off()
      
    ## c) Plot Methylation by Sex
      # graph of the raw data for percent global DNA methylaiton by maternal rank
      ggplot(data = subset(luma_data_group, !is.na(x = strank.quart)), 
             aes (x = strank.quart.order, y = meth_adjust, 
                  color = strank.quart)) +
        geom_point(shape = 1) +
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        labs (title = "Percent Global DNA 
Mehtylation by Maternal Rank") +
        ylab ("% Global DNA Methylation") +
        xlab ("Maternal Rank")
      
    ## d) boxplot of DNA methylation across quartiles of maternal rank
      boxplot(meth_adjust ~ strank.quart, data = luma_data_group, 
              main = "Boxplots of Hyena Global DNA Methylation 
              by Maternal Rank",
              xlab = "Maternal Rank Quartiles", 
              ylab = "% Global DNA Methylation", horizontal=F, 
              col=grey.colors(4))
      
    ## e)  MODEL Mehtylation by Maternal Rank Quartile 
      mat.rank.mod <- lmer(meth_adjust ~ strank.quart + (1|ID), 
                           data = luma_data_group)
      
      summary(mat.rank.mod)  # print model summary, effects and SE
      confint(mat.rank.mod)  # print 95% CIs for parameter estimates
      
      #head(model.matrix(mat.rank.mod)) #view first six rows of design matrix
      #lsmip(mat.rank.mod, ~ strank.quart) #graph estimated means
      #lsmeans(mat.rank.mod, pairwise ~ strank.quart)  # pairwise comparisons
      # plot of beta coeffiecints and SEs within each group using CAR package
      #arm::coefplot(mat.rank.mod, int=T, var.las=0, 
      #              h.axis=T, cex.pts=1, vertical=F, 
      #              main= "Treatment Constrast Esimates:
      #              Comparing Maternal Rank Quartiles", lwd=3)
      
      
      
################ Prey bivariatie meth
      
      ## 3. MODEL with PERICONCEPTION PREY
      peri_concept.mod <- lmer(meth_adjust ~ total.peri_concpt + (1|ID), 
                      data = luma_data_group)
      #display(fit.anova) # this is an alternative that shows less info than summary
      summary(peri_concept.mod)  # print model summary, effects and SE
      confint(peri_concept.mod)  # print 95% CIs for parameter estimates
      
      
      
      
      ## 3. MODEL with Gestational PREY
      gest.mod <- lmer(meth_adjust ~ total.gest + (1|ID), 
                               data = luma_data_group)
      #display(fit.anova) # this is an alternative that shows less info than summary
      summary( gest.mod)  # print model summary, effects and SE
      confint( gest.mod)  # print 95% CIs for parameter estimates
      
      
      
      
      ## 3. MODEL with Birth to 3 Months PREY
      birth_3.mod <- lmer(meth_adjust ~ total.birth.3 + (1|ID), 
                       data = luma_data_group)
      #display(fit.anova) # this is an alternative that shows less info than summary
      summary( birth_3.mod)  # print model summary, effects and SE
      confint( birth_3.mod)  # print 95% CIs for parameter estimates
      
      
      
      ## 3. MODEL with 3-6 Months PREY
      three_6.mod <- lmer(meth_adjust ~ total.3.6 + (1|ID), 
                          data = luma_data_group)
      #display(fit.anova) # this is an alternative that shows less info than summary
      summary( three_6.mod)  # print model summary, effects and SE
      confint( three_6.mod)  # print 95% CIs for parameter estimates
      
      
      
      ## 3. MODEL with 6-9 Months PREY
      six_9.mod <- lmer(meth_adjust ~ total.6.9 + (1|ID), 
                          data = luma_data_group)
      #display(fit.anova) # this is an alternative that shows less info than summary
      summary(six_9.mod)  # print model summary, effects and SE
      confint(six_9.mod)  # print 95% CIs for parameter estimates
    
      
      
      
      
      
      
      
      
      
      # graph of the raw data for percent global DNA methylaiton by maternal rank
      # stratified by sex
      ggplot(luma_data_group, aes(x = strank.quart.order, y = meth_adjust, color = Sex)) + 
        geom_point(shape = 1) +
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        geom_smooth(method = lm, se = F) + # Add linear regression best fit lines
        labs(title = "Percent Global DNA Mehtylation by Maternal 
             Rank for Males and Females", fill = "Sex" ) +
        ylab("% Global DNA Methylation") +
        xlab("Maternal Rank")
      
    
      
      # graph of the raw data for percent global DNA methylaiton by maternal rank
      # stratified by age
      ggplot(subset(luma_data_group, !is.na(x = Age)), 
             aes(x = strank.quart.order, y = meth_adjust, color = Age )) + 
        geom_point(shape = 1) +
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        geom_smooth(method = lm, se = F) + # Add linear regression best fit lines
        labs(title = "Percent Global DNA Mehtylation 
             by Maternal Rank by Age",
             fill = "age" ) +
        ylab("% Global DNA Methylation") +
        xlab("Maternal Rank")
      
   
      
      
      
      
      
      

###########################################################
####                     8 Models                      ####
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
        
        
        
        
        
        
        
        
            
