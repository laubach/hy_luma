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
    # 6: Univariate Analysis 
    # 7: Bi-Variate Analysis
    # 8: Models


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
      require(multcomp)
      require(coefplot2)
      #require(sjPlot)
      #require(sciplot)
      require(effects)
      #require(plotrix)
      #require(scatterplot3d)
      #require(texi2dvi)
      #require(stargazer)
      #require(xtable)
      require(gridExtra)
  
      
    ## c) Modeling Packages
      # Modeling Packages
      #require(nlme)
      require(lmerTest)
      #require(lme4)
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
                        , tblHyena_data. Status, FirstSeen, DenGrad, 
                          Disappeared, Mom, Birthdate, NumberLittermates, 
                          Litrank, Fate, DeathDate, Weaned, park
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
      #rename maternal rank variable from 1-4 to high, mid, low, bottom as a  
      # new variable. The strank.quart.order variable is retained to define 
      # reference and set group order for graphing
      luma_data$strank.quart <- as.factor (cut(luma_data$strank.quart.order, 
                                                   breaks = c(0, 1, 2, 3, 4), 
                                                   labels = c("Q1 (lowest)", 
                                                            "Q2", "Q3", 
                                                            "Q4 (highest)")))

      # re-order to sets the reference level for high
      luma_data <- transform (luma_data, 
                              strank.quart = factor(strank.quart,
                                                    levels = c("Q1 (lowest)",
                                                              "Q2","Q3",
                                                              "Q4 (highest)")))
      
      
  ### 5.2 Age and Sex Variables
    ## a) Re-order Age Variable
      # change Age variable to a factor from a character
      luma_data$Age <- as.factor(luma_data$Age)
      # re-order the age variable based so it is not base on not alphabetic 
      # order this sets the reference level for ANOVA to adult
      luma_data <- transform(luma_data, 
                             Age = factor(Age,
                                          levels = c("cub", 
                                                     "subadult", "adult")))

    ## b) Re-code Sex as a factor  
      luma_data$Sex <- as.factor(luma_data$Sex)
  
                                    
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
    ## c) Drop extra EstimatedAgeMo column
        luma_data_group <- luma_data_group %>%
          select (- c(EstimatedAgeMo)) 
      
    ## d) Group rows with same ID and within an age group
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
        
        
    ## e) Group rows with same ID (for calculation sex ratio etc.)
        luma_data_group_all <- luma_data_group %>% 
          filter (!is.na(meth_adjust))%>% # remove rows where meth_adjust NA
          group_by (ID) %>%    # set grouping same ID within same cat age
          group_by (meth_adjust = mean(meth_adjust)) %>%  # take avg meth_adjust
          # value for repeat ID within
          # same categorical age
          group_by (AgeMonths = mean(AgeMonths)) %>% # Avg age in months
          mutate (n = n()) %>%       # make var of n for uniqur rows            
          distinct (meth_adjust, .keep_all=TRUE) # keep all rows where var value
        # is unique while retaining all
        # all other variables    
        
      
      
###########################################################
####             6  UniVariate Analysis                ####
###########################################################          
  
  ### 6.1 Outcome Univariate 
    ## a) Descriptive Stats Outcome
      # calculate the mean, median and standard deviation of adjusted % 
      # methylation values (use ...group_all data frame)
      univar_meth = ddply(luma_data_group_all, .(), summarise,
                          N_adjust = sum(!is.na(meth_adjust)),
                          mean_adjust = mean(meth_adjust, 
                                             na.rm = T),
                          median_adjust =  quantile(meth_adjust, c(.5),
                                                    na.rm = T),
                          sd_adjust = sd(meth_adjust, na.rm = T))
        
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf("output/output_luma_ecolog/univar_meth.pdf", height = 4, width = 8)
      grid.table(univar_meth)
      dev.off()
      rm(univar_meth)
  
  ### 6.2 Univariate Stats of Predictors
    ## a) Sex Ratio summary (use ...group_all data frame)
      sex_ratio_summary <- ddply (luma_data_group_all, .(), summarise,
                        fem.sex.raito = (sum(grepl("f",
                                               luma_data_group$Sex))) /
                                    (sum(grepl("m",
                                               luma_data_group$Sex) |
                                           grepl("f",
                                                 luma_data_group$Sex))))
        
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf("output/output_luma_ecolog/sex_ratio_summary.pdf", height = 3, 
          width = 7)
      grid.table(sex_ratio_summary)
      dev.off() 
      rm(sex_ratio_summary)
      
    ## c) Age summary  
      age_var_summary <- ddply (luma_data_group, .(Age), summarise,
                                N_adjust = sum(!is.na(AgeMonths)), 
                                mean.age = round (mean(AgeMonths, 
                                                        na.rm = T),2),
                                sd.age = round (sd(AgeMonths, na.rm = T), 2))
    
    ## d) save the data frame of summary stats out as a pdf into output file
      pdf("output/output_luma_ecolog/age_var_summary.pdf", height = 3, 
          width = 7)
      grid.table(age_var_summary)
      dev.off() 
      rm(age_var_summary)
      
    ## e) Prey variable summary NEEDS TO BE FIXED
      prey_var_summary <- ddply (luma_data_group_all, .(), summarise,
                            mean.peri.prey = round (mean(total.peri_concpt,
                                                         na.rm = T), 2),
                            mean.gest.prey = round (mean(total.gest, 
                                                         na.rm = T), 2),
                            mean.birth.3.prey = round (mean(total.birth.3,
                                                            na.rm = T), 2),
                            mean.3.6.prey = round (mean(total.3.6, 
                                                        na.rm = T), 2),
                            mean.6.9.prey = round (mean(total.6.9,
                                                        na.rm = T),2))
      ## b) save the data frame of summary stats out as a pdf into output file
      pdf("output/output_luma_ecolog/prey_var_summary.pdf", height = 3, 
          width = 8)
      grid.table(prey_var_summary)
      dev.off()     
      rm(prey_var_summary)
      
  ### 6.2 Center Predictive Variables
    ## a) center function based on column means    
      center_fxn <- function(x) {
        xcenter = colMeans(x, na.rm = T)
        x - rep(xcenter, rep.int(nrow(x), ncol(x)))
      }
      
    ## b) make a list of variable names to center, value = T is necessary
      # or column positions will be returned
      vars_to_center <- grep("total", names(luma_data_group), value = T)
      
      ## c) Center data in the select columns and replace non-centered values
        # Here Prey densities have been centered
        luma_data_group[ ,c(vars_to_center)] <- center_fxn(luma_data_group
                                                     [ ,c(vars_to_center)])

  

###########################################################
####              7  Bi-Variate Analysis               ####
###########################################################          
      
  ### 7.1 cross tabulation of predicitive variables to get n per group
      table(luma_data_group_all$Sex, luma_data_group_all$strank.quart)
      table(luma_data_group$Age, luma_data_group$strank.quart)
      table(luma_data_group$Age, luma_data_group$Sex)
      table(luma_data_group$Age, luma_data_group$strank.quart, 
            luma_data_group$Sex) 
      
    
  ### 7.2 Bivariate Statistics Methylation by Sex
    ## a) Summary Stats Methylation by Sex (using ...group_all ???)
      meth_by_sex = ddply(luma_data_group_all, .(Sex), summarise,
                             N = sum(!is.na(meth_adjust)),
                             mean = round (mean(meth_adjust, na.rm = T), 2),
                             median =  round (quantile(meth_adjust, c(.5), 
                                                       na.rm = T), 2),
                             sd = round (sd(meth_adjust, na.rm = T), 2))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf("output/output_luma_ecolog/meth_by_sex.pdf", height = 6, width = 7)
      grid.table(meth_by_sex)
      dev.off()
      rm(meth_by_sex)
      
    ## c) Plot Methylation by Sex (using ...group_all ???)
      ggplot(data = subset(luma_data_group_all, !is.na(x = Sex)), 
             aes(x = Sex, y = meth_adjust, color = Sex)) + 
        geom_boxplot() +
        theme(text = element_text(size=20))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        labs(title = "Percent Global DNA 
Mehtylation by Sex") +
        ylab("% Global DNA Methylation") +
        xlab("Sex")
      
    ## d) Save Plot
      # use ggsave to save the linearization plot
      ggsave("meth_by_sex_plot.pdf", plot = last_plot(), device = NULL, 
             path = "./output/output_luma_ecolog", scale = 1, width = 7, 
             height = 5, 
             units = c("in"), dpi = 300, limitsize = TRUE)  
      
    ## e) Bivariate Regression Methylatin by Sex
      # Using lmerTest to obtain Satterthwaite approximation of DF
      sex.mod <- lmer(meth_adjust ~ Sex + (1|ID),  
                      data = subset(luma_data_group, !is.na(x = Sex)))
 
      summary(sex.mod)  # print model summary, effects and SE
      confint(sex.mod)  # print 95% CIs for parameter estimates 
     
      
  ### 7.3 Bivariate Statistics Methylation by Age 
    ## a) Summary Stats Methylation by Age (using ...group_all ???)
      meth_by_age = ddply(luma_data_group_all, .(Age), summarise,
                          N = sum(!is.na(meth_adjust)),
                          mean = round (mean(meth_adjust, na.rm = T), 2),
                          median =  round (quantile(meth_adjust, c(.5), 
                                                    na.rm = T), 2),
                          sd = round (sd(meth_adjust, na.rm = T), 2))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf("output/output_luma_ecolog/meth_by_age.pdf", height = 6, width = 7)
      grid.table(meth_by_age)
      dev.off()
      rm(meth_by_age)
      
    ## c) Plot Methylation by Age
      # graph of the raw data for percent global DNA methylaiton by maternal rank
      # stratified by age (using ...group_all ???)
      ggplot(data = subset(luma_data_group_all, !is.na(x = Age)), 
             aes(x = Age, y = meth_adjust, color = Age)) + 
        geom_boxplot() +
        theme(text = element_text(size=20))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        labs(title = "Percent Global DNA 
Mehtylation by Age") +
        ylab("% Global DNA Methylation") +
        xlab("Categorical Age")
      
    ## d) Save Plot
      # use ggsave to save the linearization plot
      ggsave("meth_by_age_plot.pdf", plot = last_plot(), device = NULL, 
             path = "./output/output_luma_ecolog", scale = 1, width = 7, 
             height = 5, 
             units = c("in"), dpi = 300, limitsize = TRUE)  
  
    ## e) Bivariate Regression Methylatino by Age 
      # Using lmerTest to obtain Satterthwaite approximation of DF
      age.mod <- lmer(meth_adjust ~ Age + (1|ID), 
                    data = subset(luma_data_group, !is.na(x = Age)))
   
      summary(age.mod)  # print model summary, effects and SE
      confint(age.mod)  # print 95% CIs for parameter estimates
      anova(age.mod) # type 3 ANOVA with Satterthwaite DF
      means.age.mod <- Effect("Age", age.mod) # use Effects to see group means

################################################################################

### PLOTTING OPTIONS #####
coefplot2(age.mod)

# First possibility
tmp <- as.data.frame(confint(glht(age.mod, mcp(Age = "Tukey")))$confint)
tmp$Comparison <- rownames(tmp)
ggplot(tmp, aes(x = Comparison, y = Estimate, ymin = lwr, ymax = upr)) +
  geom_errorbar() + geom_point()

# Second possibility
tmp <- as.data.frame(confint(glht(age.mod))$confint)
tmp$Comparison <- rownames(tmp)
ggplot(tmp, aes(x = Comparison, y = Estimate, ymin = lwr, ymax = upr)) +
  geom_errorbar() + geom_point()

# Third possibility
age.mod <- lmer(value ~ 0 + Age + (1|experiment), data = dataset)
tmp <- as.data.frame(confint(glht(age.mod))$confint)
tmp$Comparison <- rownames(tmp)
ggplot(tmp, aes(x = Comparison, y = Estimate, ymin = lwr, ymax = upr)) +
  geom_errorbar() + geom_point()

################################################################################







  ### 7.4 Bivariate Statistics Methylation by Rank
    ## a) Summary Stats Methylation by Rank
      meth_by_rank = ddply(luma_data_group, .(strank.quart), summarise,
      N = sum(!is.na(meth_adjust)),
      mean = round (mean(meth_adjust, na.rm = T), 2),
      median =  round (quantile(meth_adjust, c(.5),
      na.rm = T), 2),
      sd = round (sd(meth_adjust, na.rm = T), 2))
    
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf("output/output_luma_ecolog/meth_by_rank.pdf", height = 6, width = 7)
      grid.table(meth_by_rank)
      dev.off()

    ## c) Plot Methylation by Sex
    # graph of raw data for percent global DNA methylaiton by maternal rank
      ggplot(data = subset(luma_data_group, !is.na(x = strank.quart)),
      aes (x = strank.quart, y = meth_adjust,
      color = strank.quart)) +
      geom_boxplot() +
      theme(text = element_text(size=20))+
      scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
      labs (title = "Percent Global DNA
      Mehtylation by Maternal Rank") +
      ylab ("% Global DNA Methylation") +
      xlab ("Maternal Rank")
    
    ## d) Save Plot
    # use ggsave to save the linearization plot
      ggsave("meth_by_rank_plot.pdf", plot = last_plot(), device = NULL,
      path = "./output/output_luma_ecolog", scale = 1, width = 7,
      height = 5,
      units = c("in"), dpi = 300, limitsize = TRUE)
      
      
    
    ## e)  MODEL Mehtylation by Maternal Rank Quartile
      mat.rank.mod <- lmer(meth_adjust ~ strank.quart + (1|ID),
      data = luma_data_group)
      
      summary(mat.rank.mod)  # print model summary, effects and SE
      confint(mat.rank.mod)  # print 95% CIs for parameter estimates
      anova(mat.rank.mod) # type 3 ANOVA with Satterthwaite DF
      # use Effects for group means
      means.rank.mod <- Effect("strank.quart", mat.rank.mod)

#head(model.matrix(mat.rank.mod)) #view first six rows of design matrix
#lsmip(mat.rank.mod, ~ strank.quart) #graph estimated means
#lsmeans(mat.rank.mod, pairwise ~ strank.quart)  # pairwise comparisons
# plot of beta coeffiecints and SEs within each group using CAR package
#arm::coefplot(mat.rank.mod, int=T, var.las=0,
#              h.axis=T, cex.pts=1, vertical=F,
#              main= "Treatment Constrast Esimates:
#              Comparing Maternal Rank Quartiles", lwd=3)








################################################################################
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
  
################################################################################
  
  
  
  
  
  

###########################################################
####                     8 Models                      ####
###########################################################


    # graph of the raw data for percent global DNA methylaiton by maternal rank
    # stratified by age
    ggplot(subset(luma_data_group, !is.na(x = Age)),
    aes(x = strank.quart.order, y = meth_adjust, color = Age )) +
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
    
    ## d) Save Plot
    # use ggsave to save the linearization plot
    ggsave("meth_by_rank_by_age_lm_plot.pdf", plot = last_plot(), device = NULL,
           path = "./output/output_luma_ecolog", scale = 1, width = 7,
           height = 5,
           units = c("in"), dpi = 300, limitsize = TRUE)
    
    
    ggplot(subset(luma_data_group, !is.na(x = strank.quart)&!is.na(x=Age)),
      aes(y = meth_adjust, x = factor(strank.quart), fill = Age))+
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
    ggsave("meth_by_rank_by_age_boxplot.pdf", plot = last_plot(), device = NULL,
    path = "./output/output_luma_ecolog", scale = 1, width = 7,
    height = 5,
    units = c("in"), dpi = 300, limitsize = TRUE)
    
    
    
    
    #### Cub model
    
    ### Subset that does not include cubs or immigrant males
    ## ) Filter rows that are not cubs or immigrant males
    luma_data_cub <- luma_data_group %>%
      dplyr::filter(grepl('^cub$', Age))
    
    
    luma_data_cub$Sex <- as.factor(luma_data_cub$Sex)
    
    
    ## e)  MODEL Mehtylation by Maternal Rank Quartile
    cub.rank.mod <- lm(meth_adjust ~ strank.quart + Sex + AgeMonths,
                       data = luma_data_cub)
    
    summary(cub.rank.mod)  # print model summary, effects and SE
    confint(cub.rank.mod)  # print 95% CIs for parameter estimates
    
    levels(luma_data_cub$strank.quart)
    
    cub.rank.ef <- effect("strank.quart", cub.rank.mod)
    
    summary(cub.rank.ef)
    
    cub.rank.ef.table <- as.data.frame(cub.rank.ef)
    
    # order this sets the reference level for ANOVA to adult
    cub.rank.ef.table <- transform( cub.rank.ef.table,
                                    strank.quart = factor(strank.quart,
                                                      levels = c("Q1 (lowest)", 
                                                              "Q2", "Q3", 
                                                              "Q4 (highest)")))
    
    ggplot( cub.rank.ef.table, aes(x = strank.quart, y = fit)) +
      geom_point() +
      geom_errorbar(aes(ymin= fit-se, ymax= fit+se), width=0.4) +
      theme(text = element_text(size=20))+
      #theme_bw(base_size=12) +
      labs(title = "Cub Percent Global DNA Mehtylation 
by Maternal Rank") +
      theme(plot.title = element_text(hjust = 0.5))+
      ylab("% Global DNA Methylation ± SE") +
      xlab("Maternal Rank")
    
    ## d) Save Plot
    # use ggsave to save the linearization plot
    ggsave("mat_rank_cub_mod_beta.pdf", plot = last_plot(), device = NULL,
           path = "./output/output_luma_ecolog", scale = 1, width = 7,
           height = 5.5,
           units = c("in"), dpi = 300, limitsize = TRUE)
    
    
  
  ##### COMBINE Lowest rank quartiles
    luma_data_cub$strank.quart.comb <- as.factor(
      ifelse(luma_data_cub$strank.quart.order <= 1,
             "Q1 (lowest)", "Q2-Q4 (highest)"))
    
    
    ## e)  MODEL Mehtylation by Maternal Rank Quartile
    cub.rank.mod2 <- lm(meth_adjust ~ strank.quart.comb + Sex + AgeMonths,
                       data = luma_data_cub)
    
    summary(cub.rank.mod2)  # print model summary, effects and SE
    confint(cub.rank.mod2)  # print 95% CIs for parameter estimates
    
    
    
    
  
    
######## SUB and ADULT FEMALES
    
    ### Subset that does not include cubs or immigrant males
    ## ) Filter rows that are not cubs or immigrant males
    luma_data_no_cub_imm <- luma_data_group %>%
    dplyr::filter(!grepl('^cub$', Age)) %>%
    #dplyr::filter(grepl('^resident$', Status)) %>%
    dplyr::filter(grepl('^f$', Sex))
    
    
    
    ## e)  MODEL Mehtylation by Maternal Rank Quartile 
    # (should I control for categorical age, continuous age or an interaction)
    mat.rank.Age.mod <- lmer(meth_adjust ~ strank.quart + AgeMonths + (1|ID),
    data = luma_data_no_cub_imm)
    
    summary(mat.rank.Age.mod)  # print model summary, effects and SE
    confint(mat.rank.Age.mod)  # print 95% CIs for parameter estimates
    
    levels(luma_data_no_cub_imm$strank.quart)
    
    mat.rank.Age.ef <- effect("strank.quart", mat.rank.Age.mod)
    
    summary(mat.rank.Age.ef)
    
    mat.rank.Age.ef.table <- as.data.frame(mat.rank.Age.ef)
    
    # order this sets the reference level for ANOVA to adult
    mat.rank.Age.ef.table <- transform( mat.rank.Age.ef.table,
                                        strank.quart = factor(strank.quart,
                                                    levels = c("Q1 (lowest)",
                                                              "Q2", "Q3",
                                                              "Q4 (highest)")))
    
    
    ggplot( mat.rank.Age.ef.table, aes(x = strank.quart, y = fit)) +
    geom_point() +
    geom_errorbar(aes(ymin= fit-se, ymax= fit+se), width=0.4) +
    theme(text = element_text(size=20))+
    #theme_bw(base_size=12) +
    labs(title = "Subadult and Adult Female Percent Global
DNA Mehtylation by Maternal Rank") +
    theme(plot.title = element_text(hjust = 0.5))+
    ylab("% Global DNA Methylation ± SE") +
    xlab("Maternal Rank")
    
    
    ## d) Save Plot
    # use ggsave to save the linearization plot
    ggsave("mat_rank_sub_adult_mod_beta.pdf", plot = last_plot(), device = NULL,
    path = "./output/output_luma_ecolog", scale = 1, width = 7.25,
    height = 5.5,
    units = c("in"), dpi = 300, limitsize = TRUE)

    
    ##### COMBINE highest rank quartiles
    luma_data_no_cub_imm$strank.quart.comb <- as.factor(
    ifelse(luma_data_no_cub_imm$strank.quart.order <= 3,
    "Q1-Q3 (lowest)", "Q4 (highest"))
    
    
    ## e)  MODEL Mehtylation by Maternal Rank Quartile 
    mat.rank.Age.mod2 <- lmer(meth_adjust ~ strank.quart.comb + AgeMonths + 
                                (1|ID), data = luma_data_no_cub_imm)
    
    summary(mat.rank.Age.mod2)  # print model summary, effects and SE
    confint(mat.rank.Age.mod2)  # print 95% CIs for parameter estimates
    
   
    
    
    
    

######## SUB Females
    
    ### Subset that does not include cubs or immigrant males
    ## ) Filter rows that are not cubs or immigrant males
    luma_data_sub_fem <- luma_data_group %>%
      dplyr::filter(grepl('^subadult$', Age)) %>%
      #dplyr::filter(grepl('^resident$', Status)) %>%
      dplyr::filter(grepl('^f$', Sex))
    
    
    
    ## e)  MODEL Mehtylation by Maternal Rank Quartile 
    # (should I control for categorical age, continuous age or an interaction)
    sub.rank.mod <- lm(meth_adjust ~ strank.quart + AgeMonths,
                             data = luma_data_sub_fem)
    
    summary(sub.rank.mod)  # print model summary, effects and SE
    confint(sub.rank.mod)  # print 95% CIs for parameter estimates
    
    
    sub.rank.mod.ef <- effect("strank.quart", sub.rank.mod)
    
    summary(sub.rank.mod.ef)
    
    sub.rank.mod.ef.table <- as.data.frame(sub.rank.mod.ef)
    
    # order this sets the reference level for ANOVA to adult
    sub.rank.mod.ef.table <- transform( sub.rank.mod.ef.table,
                                        strank.quart = factor(strank.quart,
                                                      levels = c("Q1 (lowest)",
                                                              "Q2", "Q3",
                                                              "Q4 (highest)")))
    
    
    ggplot( sub.rank.mod.ef.table, aes(x = strank.quart, y = fit)) +
      geom_point() +
      geom_errorbar(aes(ymin= fit-se, ymax= fit+se), width=0.4) +
      theme(text = element_text(size=20))+
      #theme_bw(base_size=12) +
      labs(title = "Subadult Female Percent Global
DNA Mehtylation by Maternal Rank") +
      theme(plot.title = element_text(hjust = 0.5))+
      ylab("% Global DNA Methylation ± SE") +
      xlab("Maternal Rank")
    
    
    ## d) Save Plot
    # use ggsave to save the linearization plot
    ggsave("mat_rank_sub_mod_beta.pdf", plot = last_plot(), device = NULL,
           path = "./output/output_luma_ecolog", scale = 1, width = 7.25,
           height = 5.5,
           units = c("in"), dpi = 300, limitsize = TRUE)    
    
    
    
    

    
######## Adult Females
    
    ### Subset that does not include cubs or immigrant males
    ## ) Filter rows that are not cubs or immigrant males
    luma_data_adult_fem <- luma_data_group %>%
      dplyr::filter(grepl('^adult$', Age)) %>%
      #dplyr::filter(grepl('^resident$', Status)) %>%
      dplyr::filter(grepl('^f$', Sex))
    
    
    
    ## e)  MODEL Mehtylation by Maternal Rank Quartile 
    # (should I control for categorical age, continuous age or an interaction)
    adult.rank.mod <- lm(meth_adjust ~ strank.quart + AgeMonths,
                       data = luma_data_adult_fem)
    
    summary(adult.rank.mod)  # print model summary, effects and SE
    confint(adult.rank.mod)  # print 95% CIs for parameter estimates
    
    
    adult.rank.mod.ef <- effect("strank.quart", adult.rank.mod)
    
    summary(adult.rank.mod.ef)
    
    adult.rank.mod.ef.table <- as.data.frame(adult.rank.mod.ef)
    
    # order this sets the reference level for ANOVA to adult
    adult.rank.mod.ef.table <- transform(adult.rank.mod.ef.table,
                                        strank.quart = factor(strank.quart,
                                                      levels = c("Q1 (lowest)",
                                                              "Q2", "Q3",
                                                              "Q4 (highest)")))
    
    
    ggplot( adult.rank.mod.ef.table, aes(x = strank.quart, y = fit)) +
      geom_point() +
      geom_errorbar(aes(ymin= fit-se, ymax= fit+se), width=0.4) +
      theme(text = element_text(size=20))+
      #theme_bw(base_size=12) +
      labs(title = "Adult Female Percent Global
           DNA Mehtylation by Maternal Rank") +
      theme(plot.title = element_text(hjust = 0.5))+
      ylab("% Global DNA Methylation ± SE") +
      xlab("Maternal Rank")
    
    
    ## d) Save Plot
    # use ggsave to save the linearization plot
    ggsave("mat_rank_adult_mod_beta.pdf", plot = last_plot(), device = NULL,
           path = "./output/output_luma_ecolog", scale = 1, width = 7.25,
           height = 5.5,
           units = c("in"), dpi = 300, limitsize = TRUE)        
    
    
           

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









