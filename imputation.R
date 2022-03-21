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

build_the_trend_attgt <- function(gt_data, ...) {
 
  treat_data <- subset(gt_data, D==1)
  untreat_data <- subset(gt_data, D==0)

  treat_id <- treat_data$id
  temp_treat_pre <- subset(gt_data, (id %in% treat_id) & (D==0) )
  Ybar_pre <- temp_treat_pre %>%
    group_by(id) %>%
    summarize(Ybar_pre=mean(Y)) %>%
    as.data.frame()
  treat_data <- merge(treat_data, Ybar_pre, by="id")
  att1 <- mean(treat_data$Y - treat_data$Ybar_pre)

  g <- unique(treat_data$G)
  this.period <- unique(treat_data$period)

  tlist <- unique(untreat_data$period)
  tlist <- tlist[tlist <= this.period]
  nt_pre <- sum(tlist < min(g,this.period)) # handle pre- and post-periods slightly differently; in pre-treatment periods, only use up to the current period as pre-periods
  pre_count <- 2
  att2 <- 0
  
  for (tt in 2:length(tlist)) {
    # weights
    btt_w <- ifelse(tlist[tt] < min(g,this.period), (pre_count-1)/nt_pre, 1)
    pre_count <- pre_count+1
      
    dY_ids <- subset(untreat_data, period %in% tlist[tt])$id
    dY <- subset(untreat_data, period==tlist[tt] & id %in% dY_ids)$Y - subset(untreat_data, period==tlist[tt-1] & id %in% dY_ids)$Y
    att2 <- att2 + mean(dY)*btt_w
  }

  #if (g==2006 & this.period == 2005) browser()
  
  ## # alternative version as average of marcus and sant'anna using all base periods
  ## base_period_idx <- which((tlist < g) & (tlist < this.period)) # use only prior periods too
  ## this_period_idx <- which(tlist==this.period)
  ## cf_path2_by_base_period <- sapply(base_period_idx, function(base_period) {
  ##   this.att3 <- 0
  ##   for (tt in (base_period+1):this_period_idx) {
  ##     dY_ids <- subset(untreat_data, period %in% tlist[tt])$id
  ##     dY <- subset(untreat_data, period==tlist[tt] & id %in% dY_ids)$Y - subset(untreat_data, period==tlist[tt-1] & id %in% dY_ids)$Y
  ##     this.att3 <- this.att3 + mean(dY)
  ##   }
  ##   this.att3
  ## })
  ## cf_path2 <- mean(cf_path2_by_base_period)

  ## if ( !(isTRUE(all.equal(cf_path2, att2)))) browser()
         
  attgt <- list()
  attgt$ATT <- att1 - att2
  attgt_noif(attgt = attgt)
}
