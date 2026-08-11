# =============================================================================
# File: fctns_preprocess.R
# Purpose: Convert per-country HMD life table rds files into nested lists of
#          histograms, density functions, and quantile functions, then apply
#          local Fréchet regression (LFR) as a pre-smoothing step.
# Dependencies: plyr, dplyr, frechet, fdadensity
# =============================================================================

data_gen <-
  function(start_yr, end_yr, data_dir, max_age = 110) {
    ## ==============================================================
    ## Function: data_gen
    ## Purpose:  Load mortality data (male/female) by country and year,
    ##           generate lists of histograms, density functions and quantile functions, 
    ##           and compute local Fréchet regression estimated quantile functions.
    ## 
    ## Inputs:
    ##   start_yr  : Starting year (integer)
    ##   end_yr    : Ending year (integer)
    ##   data_dir  : Directory containing RDS files holding the life tables per country, 
    ##               each contains a data frame of four columns: year, sex, age, and dx.
    ##   max_age   : Maximum age to be considered, default: 110.
    ##
    ## Output:
    ##   A list containing:
    ##      subj_labs   : Labels of each subject, in this case abbreviations for each country
    ##      comp_labs   : Labels of each component, in this case labels of each sex
    ##      Lhist       : List of histograms per sex per country per year, 
    ##                   e.g., Lhist$female$aus$`1989` is a list of two fields, breaks and counts,
    ##                   holding the histogram for Australian females.
    ##      Ldens       : List of density functions per sex per country per year,
    ##                    e.g., Ldens$female$aus$`1989` is a list of three fields, bw, x, and y, 
    ##                    holding the bandwidth used, support points, and values of density functions 
    ##                    for Australian females.
    ##      qSup        : Vector of support points of quantile functions.
    ##      Lqf         : List of quantile functions per sex per country,
    ##                    e.g., Lqf$female$aus is a matrix with rows corresponding to years 
    ##                    and columns corresponding to quantile values, 
    ##                    holding the quantile functions evaluated on qSup 
    ##                    for each year from start_yr to end_yr for Australian females.
    ##      LqfLFR      : List of local Fréchet regression estimated quantile functions per sex per country,
    ##                    with the same structure as Lqf.
    ##      resLFR      : List of output of local Fréchet regression per sex per country,
    ##                    e.g., resLFR$female$aus is a list containing the LocDenReg() output 
    ##                    for Australian females.
    ##      obsGrid     : Vector holding the time grid of calendar years
    ##      workGrid    : Vector holding the time grid given \code{obsGrid} normalized to be within [0,1].
    ## ==============================================================
    
    ## ---- load RDS files ----
    List_RDS <- list.files(data_dir, pattern = "_dx\\.rds$", full.names = TRUE)
    List_using <- lapply(List_RDS, readRDS)
    
    ## ---- merge East and West Germany ----
    
    idx <- grep( 'deut', List_RDS )
    List_using[[idx[1]]] <-  
      List_using[[idx[1]]] %>%
      left_join( y = List_using[[idx[2]]], by = c('year','sex','age') ) %>% 
      mutate( dx = dx.x + dx.y ) %>% 
      select( -(dx.x:dx.y) )
    List_using <- List_using[-idx[2]]
    countries_abbr <- sub("^.*/(.*?)_dx\\.rds$", "\\1", List_RDS[-idx[2]])
    countries_abbr[idx[1]] <- "deu"
    countries_abbr <- substr( countries_abbr, 1, 3 )
    
    ## ---- take data from start_yr to end_yr with ages less than max_age, turning character-valued dx into numeric ----
    
    List_using <- llply( List_using, function (df) {
      df <- df %>% filter( year >= start_yr, year <= end_yr, age < max_age )
      if ( !is.numeric(df$dx) ) {
        df$dx <- as.numeric(df$dx)
      }
      df
    })
    
    ## ---- turn life tables into histograms ---- 
    ## Note dx = number of people who died at the age in [x,x+1), not (x,x+1] as in regular histograms.
    comp_labs <- as.vector(unique(List_using[[1]]$sex))
    Lhist <- comp_labs %>%
      alply( 1, function (sex_val) {
        List_using %>% llply( function (df_country) {
          df_country %>% filter( sex == sex_val ) %>%
            dlply( .(year), function (df_country_sex_yr) {
              df_country_sex_yr %>% with(list(
                breaks = c(age,max_age),
                counts = dx
              ))
            })
        }) %>% set_names(countries_abbr)
      }) %>% set_names(comp_labs)
    ## e.g., Lhist$female$aus$`1989` is a list of two fields, breaks and counts,
    ## holding the histogram for Australian females.
    
    ## ---- turn histograms into density functions ----
    dSup <- seq(0,max_age,0.2)
    Ldens <- Lhist %>% llply( function (Lhist_sex) {
      Lhist_sex %>% llply( function (Lhist_sex_country) {
        Lhist_sex_country %>% llply( function (hist_sex_country_yr) {
          frechet::CreateDensity(
            histogram = hist_sex_country_yr,
            optns = list( outputGrid = dSup, infSupport = FALSE )
          )
        })
      })
    })
    ## e.g., Ldens$female$aus$`1989` is a list of three fields, bw, x, and y, 
    ## holding the density function for Australian females.
    
    ## ---- turn density functions into quantile functions ----
    qSup <- seq(0,1,1/200)
    Lqf <- Ldens %>% llply( function (Ldens_sex) {
      Ldens_sex %>% llply( function (Ldens_sex_country) {
        Ldens_sex_country %>% laply( function (dens_sex_country_yr) {
          fdadensity::dens2quantile(
            dens = dens_sex_country_yr$y,
            dSup = dens_sex_country_yr$x,
            qSup = qSup
          )
        })
      })
    })
    ## e.g., Lqf$female$aus is a matrix with rows corresponding to years 
    ## and columns corresponding to quantile values, 
    ## holding the quantile functions for each year from start_yr to end_yr 
    ## for Australian females.
    
    ## ---- local Fréchet regression (LFR) ----
    xout <- xin <- start_yr:end_yr
    resLFR <- Lqf %>% llply( function (Lqf_sex) {
      Lqf_sex %>% llply( function (qf_sex_country) {
        frechet::LocDenReg(
          xin = xin, qin = qf_sex_country, xout = xout,
          optns = list(
            qSup = qSup,
            lower = 0, upper = max_age,
            kernelReg = "epan", bwReg = "CV"
          ) 
        )
      })
    })
    ## e.g., resLFR$female$aus is a list containing the LocDenReg() output 
    ## for Australian females.
    
    ## ---- extract LFR estimated quantile functions ----
    LqfLFR <- resLFR  %>% llply( function (res_sex) {
      res_sex %>% llply( function (res_sex_country) {
        res_sex_country$qout
      })
    })
    ## e.g., LqfLFR$female$aus is a matrix with rows corresponding to years 
    ## and columns corresponding to quantile values
    ## holding the LFR estimated quantile functions for each year 
    ## from start_yr to end_yr for Australian females.
    
    
    ## ---- return output ----
    return(list(
      subj_labs = countries_abbr,
      comp_labs = comp_labs,
      Lhist = Lhist,
      Ldens = Ldens,
      qSup = qSup,
      Lqf = Lqf,
      LqfLFR = LqfLFR,
      resLFR = resLFR,
      obsGrid = start_yr:end_yr,
      workGrid = ((start_yr:end_yr) - start_yr) / (end_yr - start_yr)
    ))
  }
