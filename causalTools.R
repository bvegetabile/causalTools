################################################################################
# 
# Tools for performing causal analysis
# --- Author: Brian Vegetabile (bvegetab [AT] uci [DOTEDU])
#
################################################################################

#-------------------------------------------------------------------------------
# Balance Tools for assessing covariate balance
# --- bal_stats calculates balance statistics (wtd or not) for single variable
# --- bal_table calculates balance statistics for dataset
#-------------------------------------------------------------------------------

bal_stats <- function(var_data,
                     treat_ind,
                     datatype = 'continuous',      
                     wts=rep(1, length(treat_ind))){
  #-----------------------------------------------------------------------------
  # bal_stats is a function which is used to assess covariate balance for a 
  # single variable.  
  #
  # Input variables
  # ---  var_data  : vector of data for a single variable
  # ---  treat_ind : vector of treatment indicators.  Binary or Logical
  # ---  datatype  : type of data provided.  allows 'continuous' or 'binary'
  # ---  wts       : weights after estimating propensity score or similar.
  #                  defaults to 1 for all observations unless otherwise 
  #                  specified.  
  #-----------------------------------------------------------------------------
  
  # Converting 0/1's to logicals to ensure subsetting works as desired
  treat_ind <- as.logical(treat_ind)
  
  #-----------------------------------------------------------------------------
  # Weighted mean in the treated and control groups
  
  # Weighted sample mean in treated group
  xbar_t <- sum(wts[treat_ind] * var_data[treat_ind]) / sum(wts[treat_ind])
  # Weighted sample mean in control group
  xbar_c <- sum(wts[!treat_ind] * var_data[!treat_ind]) / sum(wts[!treat_ind])
  
  #-----------------------------------------------------------------------------
  # Weighted sample variance calculations in both groups for binary and 
  # continuous data.  Number of samples in each group also calculated. 
  
  if(datatype=='binary'){
    # Binary Sample Variance - Treatment Group
    s2_t <- xbar_t*(1-xbar_t) 
    # Binary Sample Variance - Control Group
    s2_c <- xbar_c*(1-xbar_c)
    
    # Sample sizes for binary data
    NT <- sum(as.logical(var_data) & treat_ind)
    NC <- sum(as.logical(var_data) & !treat_ind)
    
  } else if(datatype=='continuous') {
    # Weighted Sample Variance in Treatment Group
    s2_t <- sum(wts[treat_ind] * (var_data[treat_ind]- xbar_t)**2) / sum(wts[treat_ind])
    # Weighted Sample Variance in Control Group
    s2_c <- sum(wts[!treat_ind] * (var_data[!treat_ind]- xbar_c)**2) / sum(wts[!treat_ind])
    
    # Sample sizes for continuous data   
    NT <- sum(treat_ind)
    NC <- sum(!treat_ind)
  }
  
  #-----------------------------------------------------------------------------
  # Calculating relevant statistics.  
  # -- Standardized Difference in Means for Balance.  Want |std_diff| < 0.1
  # -- Log of ratio of Standard Deviation.  Checks balance of second moment. 
  #    Should be close to zero.
  
  std_diff <- (xbar_t - xbar_c) / sqrt((s2_t + s2_c) / 2)
  log_var_ratio <- log(sqrt(s2_t)) - log(sqrt(s2_c))
  
  return(c(NT, round(xbar_t,4), round(s2_t,4), 
           NC, round(xbar_c,4), round(s2_c,4), 
           round(std_diff,4), round(log_var_ratio, 4)))
}


bal_table <- function(dataset, 
                      col_ind, 
                      treat_ind, 
                      wts = rep(1, length(treat_ind)), 
                      max_uniq=5){
  #-----------------------------------------------------------------------------
  # bal_table is a function which provides a table of covariate balance stats. 
  #
  # Input variables
  # --- dataset   : matrix or dataframe of data
  # --- col_ind   : indices of columns to check balance. vector of integers
  # --- treat_ind : indicator of treatment assignment. binary or logical
  # --- wts       : weights for use after estimating the propensity score
  # --- max_uniq  : defines threshold for what to consider a continuous or 
  #                 categorical variable.  If the number of unique members of a
  #                 category is less than this number it will be converted to a 
  #                 factor if not already a factor.  
  #-----------------------------------------------------------------------------
  var_names <- names(dataset[,col_ind])
  treat_ind <- as.logical(treat_ind)
  outtable <- c()
  counter <- 1
  
  # Iteration based on the order of column numbers provided
  for(i in 1:length(col_ind)){
    c <- col_ind[i]
    col_data <- dataset[, c]
    # Conversion of variables to categorical if it has less than the max_uniq
    # categories.  
    if(length(unique(col_data)) <= max_uniq){
      col_data <- as.factor(col_data)
    }
    if(is.factor(col_data)){
      obs_lvls <- levels(col_data)
      outtable <- rbind(outtable, rep(NA, 8))
      row.names(outtable)[counter] <- var_names[i]
      counter <- counter + 1
      for(lvl in obs_lvls){
        stddiff <- bal_stats(col_data==lvl, treat_ind, 'binary', wts)  
        outtable <- rbind(outtable, stddiff)
        row.names(outtable)[counter] <- lvl
        counter <- counter + 1
      }
    } else {
      stddiff <- bal_stats(col_data, treat_ind, 'continuous', wts)  
      outtable <- rbind(outtable, stddiff)
      row.names(outtable)[counter] <- var_names[i]
      counter <- counter + 1
    }
  }
  colnames(outtable) <- c('NT', 'MeanT', 'VarT', 
                          'NC', 'MeanC', 'VarC', 
                          'StdDiff', 'LogRatio')
  return(outtable)
}


wtd_ecdf <- function (var_data, wts) {
  #-----------------------------------------------------------------------------
  # wtd_ecdf is a modification of the ecdf() function in base R.  It modifies
  # the function to be able to incorporate weights.  This is to visualize 
  # balance using the empirical cumulative distribution function for continuous 
  # covariates after weighting by the inverse of the propensity score (IPTW)
  # 
  # Input variables
  # --- var_data : covariate values - vector of data
  # --- wts      : weights for assessing cov balance by IPTW - vector of data.
  #-----------------------------------------------------------------------------
  ord <- order(var_data)
  var_ordered <- var_data[ord]
  wts_ordered <- wts[ord]
  n <- length(var_data)
  if (n < 1)
    stop("'var_data' must have 1 or more non-missing values")
  vals <- unique(var_ordered)
  matched_vals <- match(var_ordered, vals)
  weight_list <- aggregate(wts_ordered, by=list(matched_vals), sum)
  rval <- approxfun(vals, cumsum(weight_list[,2])/sum(wts_ordered),
                    method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}
