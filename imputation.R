library(pte)
library(BMisc)
library(tidyr)
library(dplyr)


data4_imp <- data4
data4_imp$G[data4_imp$G==2007] <- 0
data4_imp <- subset(data4_imp, year < 2007)
data4_imp <- droplevels(data4_imp)

imp2 <- pte(yname="lemp",
            tname="year",
            idname="countyreal",
            gname="G",
            data=data2,
            setup_pte_fun=setup_pte_basic,
            subset_fun=keep_all_subset,
            attgt_fun=imputation_attgt,
            xformla=~as.factor(region)*as.factor(year) + lpop + lavg_pay,
            biters=100,
            cl=1)

summary(imp2)


btt2 <- pte(yname="lemp",
            tname="year",
            idname="countyreal",
            gname="G",
            data=data2,
            setup_pte_fun=setup_pte_basic,
            subset_fun=keep_all_subset,
            attgt_fun=build_the_trend_attgt,
            xformla=~1,
            biters=10,
            cl=1)

summary(btt2)
ggpte(btt2)

keep_all_subset <- function(data, g, tp, ...) {
  # this would probably need to be updated for multiplier bootstrap

  # drop post-treatment observations that are not in group g
  data$.treat <- 1*((data$G <= data$period) & (data$G != 0))
  this.data <- subset(data, (.treat==0) | (G==g))
  # not clear if this is totally right in pre-treatment periods
  # as it drops any periods after the "current" period too (even
  # if these are untreated)
  this.data <- subset(this.data, !(G==g & period > tp)) 
  this.data$D <- 1*((this.data$G==g) & this.data$period >= tp) 
  n1 <- length(unique(this.data$id))
  disidx <- unique(data$id) %in% unique(this.data$id) 
  list(gt_data = this.data, n1 = n1, disidx = disidx)
}


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


