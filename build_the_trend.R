#-----------------------------------------------------------------------------
#
# This file contains functions for computing the build-the-trend
# version of Callaway and Sant'Anna discussed in the chapter
# using the infrastructure of the `pte` package
#
#-----------------------------------------------------------------------------

#' Function to compute the build the trend estimator
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
  # handle pre- and post-periods slightly differently;
  # in pre-treatment periods, only use up to the current period
  # as pre-periods
  nt_pre <- sum(tlist < min(g,this.period)) 
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
    
  attgt <- list()
  attgt$ATT <- att1 - att2
  attgt_noif(attgt = attgt)
}
