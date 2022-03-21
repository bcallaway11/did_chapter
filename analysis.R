#-----------------------------------------------------------------------------
#
# Main file for DID chapter
#
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# load external packages
#-----------------------------------------------------------------------------
library(did)
library(BMisc)
library(fixest)
library(bacondecomp)
library(ggplot2)
library(modelsummary)
library(HonestDiD)

#-----------------------------------------------------------------------------
# load/process data
#-----------------------------------------------------------------------------
load("mw_data_ch.RData")
data <- mw_data_ch
rm(mw_data_ch)
data <- subset(data, emp0A01_BS > 0)
data <- subset(data, annual_avg_pay > 0)
data <- data[complete.cases(data[,-c("state_mw", "fed_mw")]),]
data <- makeBalancedPanel(data, "countyreal", "year")
data$lemp <- log(data$emp0A01_BS)
data$pop <- as.numeric(data$pop)
data$lpop <- log(data$pop)
data$lavg_pay <- log(data$annual_avg_pay)
data$region <- 1*(data$censusdiv==1 | data$censusdiv==2) + 2*(data$censusdiv==3 | data$censusdiv==4) +
    3*(data$censusdiv==5 | data$censusdiv==6 | data$censusdiv==7) + 4*(data$censusdiv==8 | data$censusdiv==9)
data$region <- as.factor(data$region)
data$ever_treated <- 1*(data$G != 0)

#-----------------------------------------------------------------------------
# creat subsets of data used in chapter
#-----------------------------------------------------------------------------
# drops already and early-treated, similar to data used in CS
data2 <- subset(data, (G==0) | (G>2001))

# drop never-treated 
data3 <- subset(data, G!=0)

# additionally drop always-treated
data4 <- subset(data3, G > 2001)

# drop years 2004 or earlier
data5 <- subset(data4, year > 2004)


#-----------------------------------------------------------------------------
# summary statistics
#-----------------------------------------------------------------------------
data2$fever_treated <- as.factor(data2$ever_treated)
data2$emp100 <- data2$emp0A01_BS/100
data2$pop1000 <- data2$pop/1000
data2$pay1000 <- data2$annual_avg_pay/1000
datasummary_balance(~ fever_treated,
                    data=dplyr::select(subset(data2,year==2001),
                                       emp100,
                                       pop1000,
                                       pay1000,
                                       region,
                                       fever_treated),
                    output="latex")


#-----------------------------------------------------------------------------
# Callaway and Sant'Anna estimates
# * Different results come from changing the `data` argument among
#   `data`, `data2`, `data3`, `data4`, and `data5`; from changing
#   the `base_period` among "varying" and "universal"; and the
#   `control_group` among "nevertreated" and "notyettreated".
#-----------------------------------------------------------------------------
cs_res <- att_gt(yname="lemp",
                 tname="year",
                 idname="countyreal",
                 gname="G",
                 data=data2,
                 base_period ="varying",
                 control_group="nevertreated",
                 cband=FALSE)

# print att_gt results
cs_res

# event study aggregation
cs_dyn <- aggte(cs_res, type="dynamic")
cs_dyn

# overall treatment effect aggregation
cs_o <- aggte(cs_res, type="group")
cs_o

# generates plots of att_gt, event study, and overall treatment effect
ggdid(cs_res, ylim=c(-.18,.1))
ggdid(cs_dyn)
ggdid(cs_o)

# results under conditional parallel trends
# * can apply same options here as for unconditional parallel trends above 
csX_res <- att_gt(yname="lemp",
                 tname="year",
                 idname="countyreal",
                 gname="G",
                 data=subset(data2, region !=1),
                 xformla=~lpop + lavg_pay + region,
                 base_period ="varying",
                 control_group="nevertreated")
csX_res

csX_dyn <- aggte(csX_res, type="dynamic", na.rm=TRUE)
csX_o <- aggte(csX_res, type="group", na.rm=TRUE)
ggdid(csX_dyn)

# twfe regression using different data
twfe1 <- feols(lemp ~ treated | countyreal + year,
               data=data)
twfe2 <- feols(lemp ~ treated | countyreal + year,
               data=data2)
twfe3 <- feols(lemp ~ treated | countyreal + year,
               data=data3)
twfe4 <- feols(lemp ~ treated | countyreal + year,
               data=data4)
twfe5 <- feols(lemp ~ treated | countyreal + year,
               data=data5)

twfe_covs <- feols(lemp ~ treated + as.factor(year)*as.factor(region) + lpop + lavg_pay | countyreal, data=data4)

modelsummary(list(twfe1, twfe2, twfe3, twfe4, twfe5),
             output="latex",
             stars=TRUE)

# run bacon decomposition
bacon_res <- bacon(lemp ~ post, 
                   data=data5,
                   id_var="countyreal",
                   time_var="year")

# confirm same estimate
sum(bacon_res$estimate * bacon_res$weight)

# bacon decomp results
bacon_res

# get weights by type
aggregate(bacon_res$weight, by=list(bacon_res$type), sum)

# find fraction of "good" weight on particular group
this_group <- 2007
good_weight <- subset(bacon_res, type %in% c("Earlier vs Later Treated", "Treated vs Untreated"))

sum(good_weight[good_weight$treated==this_group,]$weight) / sum(good_weight$weight)

sum(bacon_res[bacon_res$treated==this_group,]$weight)


mean(data2$G==this_group)/mean(data2$G>0)

mean(data2$G==2007)/mean(data2$G>0)

# plot
ggplot(data=bacon_res, 
       mapping=aes(x=weight,
                   y=estimate,
                   color=as.factor(type))) + 
  geom_point(size=5) + 
  scale_color_discrete(name="") + 
  theme_classic() + 
  theme(legend.position="bottom")



# twfe event study
this_data <- data2

this_data$e <- ifelse(this_data$G==0, 0, this_data$year - this_data$G)

twfe_es <- feols(lemp ~ i(e, ever_treated, ref=-1) | countyreal + year,
                 data=this_data)
summary(twfe_es)
iplot(twfe_es)

# try to construct event study regression weigths
Dit <- i(this_data$e, this_data$ever_treated, ref=-1)

ddotDit <- demean(Dit, this_data[,c("countyreal", "year")])

Yit <- this_data$lemp

bete <- solve(t(ddotDit)%*%Dit) %*% t(ddotDit) %*% Yit
bete

# weight reg
e <- 0# pick which event time
g <- 2007
es_weight_outcome <- 1*( (this_data$G == g) & (g+e == this_data$year))

es_weight_e <- solve(t(ddotDit)%*%Dit) %*% t(ddotDit) %*% es_weight_outcome
es_weight_e

# get the att_gt's
attgt_idx <- cs_res$group == g & cs_res$t == g+e
attgt <- cs_res$att[attgt_idx]
cs_dyn$att.egt[cs_dyn$egt==e]

get_attgt <- function(g,e) {
  attgt_idx <- cs_res$group == g & cs_res$t == g+e
  if( !(any(attgt_idx)) ) return(0)
  attgt <- cs_res$att[attgt_idx]
  attgt
}

p_g <- function(g) {
  mean(this_data$G==g)
}

tlist <- unique(this_data$year)
nT <- length(tlist)
he <- function(e,g,t) {
  ( 1*(g+e==t) - 1*( (g+e) %in% tlist )/nT) * (g != 0)
}
e.seq <- sort(unique(this_data$e[this_data$e != -1]))
hinner <- function(g,t) {
  sapply(e.seq, he, g=g,t=t)
}
h <- function(g,t) {
  hinner(g,t) - apply(sapply(glist, function(gg) hinner(gg,t)*p_g(gg)), 1, sum)
}

n <- length(unique(this_data$countyreal))
inner <- function(g,e,weightonly=FALSE) {
  part1 <- solve(t(ddotDit)%*%Dit/n)
  part2 <- as.matrix(h(g,g+e))*( (g+e) %in% tlist )*(g!=0)*p_g(g)
  if (!weightonly) part2 <- part2*get_attgt(g,e)
  part1 %*% part2
  #browser()
  #1+1
  ## # sa weights
  ## # this appears to work...
  ## data2$weight_outcome <- data2$ever_treated*(data2$e==e)*(data2$G==g)
  ## if(get_attgt(g,e)==0) return(as.matrix(rep(0, length(e.seq))))
  ## as.matrix(coef(feols( weight_outcome ~ i(e, ever_treated, ref=-1) | countyreal + year,
  ##               data=data2)))*get_attgt(g,e)
}

glist <- sort(unique(this_data$G))[-1]
coff <- function() {
  out <- matrix(data=0, nrow=length(e.seq))
  for (e in e.seq) {
    for (g in glist) {
      out <- out + inner(g,e)
      #for (t in tlist) {
      #  out <- out + inner(e,g,t)
      #}
    }
  }
  out
}
coff()

eval <- 5
eval.idx <- which(e.seq==eval)
eg <- expand.grid(e.seq, glist)
e_weight <- sapply(1:nrow(eg), function(i) inner(eg[i,2], eg[i,1], weightonly=TRUE)[eval.idx])
sum(e_weight[eg[,1]==eval])
sum(e_weight[eg[,1]!=eval])
max(abs(e_weight[eg[,1] != eval]))



# weight reg from SA
coef(feols( I(ever_treated*(e==0)*(G==2007)) ~ i(e, ever_treated, ref=-1) | countyreal + year,
                 data=this_data))


# put on same plot as callaway and santanna
plot_df1 <- cbind.data.frame(att=cs_dyn$att.egt, ciL=cs_dyn$att.egt-1.96*cs_dyn$se.egt, ciU=cs_dyn$att.egt+1.96*cs_dyn$se.egt, e=cs_dyn$egt, type="CS")
twfe_plot_data <- iplot(twfe_es)
plot_df2 <- cbind.data.frame(att=twfe_plot_data$prms$estimate,
                             ciL=twfe_plot_data$prms$ci_low,
                             ciU=twfe_plot_data$prms$ci_high,
                             e=twfe_plot_data$prms$estimate_names,
                             type="ES Reg. Data 2")
# event study with full data to put on plot
data$e <- ifelse(data$G==0, 0, data$year - data$G)
twfe_es_all <- feols(lemp ~ i(e, ever_treated, ref=-1) | countyreal + year,
                 data=data)
twfe_all_plot_data <- iplot(twfe_es_all)
plot_df3 <- cbind.data.frame(att=twfe_all_plot_data$prms$estimate,
                             ciL=twfe_all_plot_data$prms$ci_low,
                             ciU=twfe_all_plot_data$prms$ci_high,
                             e=twfe_all_plot_data$prms$estimate_names,
                             type="ES Reg. Data 1")

plot_df <- rbind.data.frame(plot_df1, plot_df2, plot_df3)

ggplot(data=plot_df, mapping=aes(x=e,y=att,color=type,fill=type)) +
  geom_point(position=position_dodge(width=.2)) +
  geom_errorbar(aes(ymin=ciL, ymax=ciU), position=position_dodge(width=.2), width=.5) +
  theme_classic() +
  scale_x_continuous(breaks=seq(-6,6)) +
  theme(legend.position="bottom", legend.title=element_blank()) 
#ggsave("cs_twfe_es_comparison.pdf", height=5, width=8)

# twfe event study with covariates
twfeX <- feols(lemp ~ treated + lpop + lavg_pay | countyreal + year + region^year,
               data=data)
summary(twfeX)

twfeX_es <- feols(lemp ~ i(e, ref=-1) + lpop + lavg_pay  | countyreal + year + region^year,
                 data=data)
summary(twfeX_es)
iplot(twfeX_es)



#-----------------------------------------------------------------------------
#
# now, let's drop the never-treated group
#
#-----------------------------------------------------------------------------
data_no_untreated <- subset(data, G > 0)

nu_twfe <- feols(lemp ~ treated | countyreal + year,
              data=data_no_untreated)
summary(nu_twfe)

nu_twfe_es <- feols(lemp ~ i(e, ref=-1) | countyreal + year,
                 data=data_no_untreated)
summary(nu_twfe_es)
iplot(nu_twfe_es)


# twfe event study with covariates
nu_twfeX <- feols(lemp ~ treated + lpop + lavg_pay  | countyreal + year + region^year,
               data=data_no_untreated)
summary(nu_twfeX)

nu_twfeX_es <- feols(lemp ~ i(e, ref=-1) + lpop + lavg_pay  | countyreal + year + region^year,
                 data=data_no_untreated)
summary(nu_twfeX_es)
iplot(nu_twfeX_es)


#-----------------------------------------------------------------------------
#
# now, let's drop the already treated
#
#-----------------------------------------------------------------------------
data_no_already_treated <- subset(data, (G==0) | (G>2002))

nat_twfe <- feols(lemp ~ treated | countyreal + year,
              data=data_no_already_treated)
summary(nat_twfe)

nat_twfe_es <- feols(lemp ~ i(e, ref=-1) | countyreal + year,
                 data=data_no_already_treated)
summary(nat_twfe_es)
iplot(nat_twfe_es)


# twfe event study with covariates
nat_twfeX <- feols(lemp ~ treated + lpop + lavg_pay  | countyreal + year + region^year,
               data=data_no_already_treated)
summary(nat_twfeX)

nat_twfeX_es <- feols(lemp ~ i(e, ref=-1) + lpop + lavg_pay  | countyreal + year + region^year,
                 data=data_no_already_treated)
summary(nat_twfeX_es)
iplot(nat_twfeX_es)



#-----------------------------------------------------------------------------
#
# honest did
#
#-----------------------------------------------------------------------------

hd_rm <- honest_did(cs_dyn, type="relative_magnitude", Mbarvec=c(0.5,1,1.5,2), gridPoints=100, grid.lb=-.25, grid.ub=.25)

createSensitivityPlot_relativeMagnitudes(hd_rm$robust_ci,
                                         hd_rm$orig_ci) +
  theme_classic() +
  scale_color_discrete("", labels=c("Original", "RR Bounds")) +
  theme(legend.position="bottom")
#ggsave("rr_bounds.pdf", width=8, height=5)

hd_lt <- honest_did(cs_dyn, type="relative_magnitude", Mbarvec=c(0.5,1,1.5,2), bound="deviation from linear trend", gridPoints=100, grid.lb=-.25, grid.ub=.25)
createSensitivityPlot_relativeMagnitudes(hd_lt$robust_ci,
                                         hd_lt$orig_ci) + theme_classic()

hd_rm1 <- honest_did(es=cs_dyn, type="relative_magnitude", Mbarvec=c(0.5,1,1.5,2), e=1, gridPoints=100, grid.lb=-1, grid.ub=1)
createSensitivityPlot_relativeMagnitudes(hd_rm1$robust_ci,
                                         hd_rm1$orig_ci) + theme_classic()

hd_rm2 <- honest_did(es=cs_dyn, type="relative_magnitude", Mbarvec=c(0.5,1,1.5,2), e=2, gridPoints=100, grid.lb=-1, grid.ub=1)
createSensitivityPlot_relativeMagnitudes(hd_rm2$robust_ci,
                                         hd_rm2$orig_ci) + theme_classic()

hd_lt2 <- honest_did(es=cs_dyn, type="relative_magnitude", Mbarvec=c(0.5,1,1.5,2), bound="deviation from linear trend", e=2, gridPoints=100, grid.lb=-1, grid.ub=1)
createSensitivityPlot_relativeMagnitudes(hd_lt2$robust_ci,
                                         hd_lt2$orig_ci) + theme_classic()



hd_smooth <- honest_did(cs_dyn, type="smoothness")
createSensitivityPlot(hd_lt$robust_ci,
                      hd_lt$orig_ci)

hd_smooth2 <- honest_did(es=cs_dyn, e=2, type="smoothness")
createSensitivityPlot(hd_lt2$robust_ci,
                      hd_lt2$orig_ci)


#-----------------------------------------------------------------------------
#
# imputation estimator
#
#-----------------------------------------------------------------------------

compute.did_imputation <- function(upo_formula,
                                   data,
                                   treatedname,
                                   gname,
                                   event_times) {

  untreated_data <- data[data[,treatedname]==0,]

  untreated_twfe <- feols(upo_formula,
                          data=untreated_data)

  data$.imputation <- predict(untreated_twfe, newdata=data)

  imputation_event_study <- function(event_time) {
    this_data <- subset(data, e==event_time)
    mean(this_data$lemp - this_data$.imputation, na.rm=TRUE)
  }

  ies <- sapply(event_times, imputation_event_study)

  #glist <- sort(unique(data[,gname]))
  
  return(list(imputation_event_study=ies, imputation_overall_att=0))
}


data <- as.data.frame(data)
compute.did_imputation(lemp ~ lpop + lavg_pay | countyreal + year + region^year, data=subset(data,G!=2001), treatedname="treated", gname="G", event_times=-6:6)

did_imputation <- function(upo_formula,
                           data,
                           treatedname,
                           gname,
                           idname,
                           event_times,
                           biters=100,
                           cl=1) {

  est <- compute.did_imputation(upo_formula=upo_formula,
                                data=data,
                                treatedname=treatedname,
                                gname=gname,
                                event_times=event_times)

  
  # bootstrap
  cat("bootstrapping...\n")
  boot_out <- pblapply(1:biters, function(b) {
    boot_data <- blockBootSample(data, idname)
    boot_est <- compute.did_imputation(upo_formula=upo_formula,
                                       data=boot_data,
                                       treatedname=treatedname,
                                       gname=gname,
                                       event_times=event_times)
    boot_est
  }, cl=cl)

  boot_es <- t(simplify2array(BMisc::getListElement(boot_out, "imputation_event_study")))
  es_se <- apply(boot_es, 2, sd)
  o_se <- 0
  return(list(es=est$imputation_event_study, overall=est$imputation_overall_att, es_se=es_se, o_se=o_se))
}

did_imp <- did_imputation(lemp ~ lpop + lavg_pay | countyreal + year + region^year,
                          data=subset(data,G!=2001),
                          treatedname="treated",
                          gname="G",
                          idname="countyreal",
                          event_times=-6:6,
                          biters=100)

imp_plot_df <- cbind.data.frame(att=did_imp$es, ciL=did_imp$es-1.96*did_imp$es_se, ciU=did_imp$es + 1.96*did_imp$es_se, e=-6:6, type="Imputation")
plot_df1 <- cbind.data.frame(att=csX_dyn$att.egt, ciL=csX_dyn$att.egt-1.96*csX_dyn$se.egt, ciU=csX_dyn$att.egt+1.96*csX_dyn$se.egt, e=csX_dyn$egt, type="CS")
plot_df <- rbind.data.frame(plot_df1, imp_plot_df)

ggplot(data=plot_df, mapping=aes(x=e,y=att,color=type,fill=type)) +
  geom_point(position=position_dodge(width=.2)) +
  geom_errorbar(aes(ymin=ciL, ymax=ciU), position=position_dodge(width=.2), width=.5) +
  theme_classic() +
  scale_x_continuous(breaks=seq(-6,6)) +
  theme(legend.position="bottom")

ggsave("cs_imp_es_comparison_covs.pdf", height=5, width=8)
