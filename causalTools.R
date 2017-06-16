#-------------------------------------------------------------------------------
# 
# Tools for performing causal analysis
# --- Author: Brian Vegetabile (bvegetab [AT] uci [DOTEDU])
# 
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Balance Tools for assessing covariate balance
# --- bal_stats calculates balance statistics (wtd or not) for single variable
# --- bal_table calculates balance statistics for dataset

bal_stats <- function(var_data,
                     treat_ind,
                     datatype = 'continuous',      
                     wts=rep(1, length(treat_ind))){
  #-----------------------------------------------------------------------------
  # bal_stats is a function which is used to assess covariate balance for a 
  # single variable.  
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