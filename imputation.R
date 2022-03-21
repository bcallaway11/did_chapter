#-----------------------------------------------------------------------------
#
# This file contains functions for computing ATT(g,t) using
# imputations and building on the infrastructure of the `pte`
# packge
#
#-----------------------------------------------------------------------------


#' Function to keep all available data for computing ATT(g,t)
#'
#' @param data a data frame
#' @param g group
#' @param tp time period
#'
#' @return all data but in correct format for computing ATT(g,t)
keep_all_subset <- function(data, g, tp, ...) {
  # drop post-treatment observations that are not in group g
  data$.treat <- 1*((data$G <= data$period) & (data$G != 0))
  this.data <- subset(data, (.treat==0) | (G==g))
  # This line drops any periods after the "current" period too (even
  # if these are untreated); it seems ambiguous/confusing what
  # exactly to do in pre-treatment periods.
  this.data <- subset(this.data, !(G==g & period > tp)) 
  this.data$D <- 1*((this.data$G==g) & this.data$period >= tp) 
  n1 <- length(unique(this.data$id))
  disidx <- unique(data$id) %in% unique(this.data$id) 
  list(gt_data = this.data, n1 = n1, disidx = disidx)
}

#' Function to impute ATT(g,t)
#'
#' @param gt_data data for this group and time period
#' @param xformla formulat to use for the covariates
#'
#' @return estimate of ATT(g,t)
imputation_attgt <- function (gt_data, xformla, ...) {
  treat_data <- subset(gt_data, D==1)
  untreat_data <- subset(gt_data, D==0)
  formla <- toformula("Y", rhs.vars(xformla))
  untreated_reg <- feols(formla, fixef=c("period", "id"), data=untreat_data)
  treat_data$.y0 <- predict(untreated_reg, newdata=treat_data)
  attgt <- list()
  attgt$ATT <- mean(treat_data$Y - treat_data$.y0)
  attgt_noif(attgt = attgt)
}


