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
library(pte)
library(tidyr)
library(dplyr)

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
data$id <- data$countyreal

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
                 base_period ="universal",
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

# print group-time average treatment effects
csX_res

# event study
csX_dyn <- aggte(csX_res, type="dynamic", na.rm=TRUE)
csX_dyn

# overall ATT
csX_o <- aggte(csX_res, type="group", na.rm=TRUE)
csX_o

# plot event study (other plots similar)
ggdid(csX_dyn)

#-----------------------------------------------------------------------------
# twfe regressions using different data
#-----------------------------------------------------------------------------
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

# twfe with covariates; in chapter this is reported
# with for `data2` and `data4`
twfe_covs <- feols(lemp ~ treated + as.factor(year)*as.factor(region) + lpop + lavg_pay | countyreal, data=data2)

# print results
modelsummary(list(twfe1, twfe2, twfe3, twfe4, twfe5),
             output="latex",
             stars=TRUE)

#-----------------------------------------------------------------------------
# bacon decomposition
#-----------------------------------------------------------------------------
bacon_res <- bacon(lemp ~ post, 
                   data=data2,
                   id_var="countyreal",
                   time_var="year")

# confirm same estimate
sum(bacon_res$estimate * bacon_res$weight)

# bacon decomp results
bacon_res

# get weights by type
aggregate(bacon_res$weight, by=list(bacon_res$type), sum)

# find fraction of weight on particular group
this_group <- 2006
sum(bacon_res[bacon_res$treated==this_group,]$weight)

# compare to relative size of group
mean(data2$G==this_group)/mean(data2$G>0)

# bacon decomposition plot; this is not reported in chapter
ggplot(data=bacon_res, 
       mapping=aes(x=weight,
                   y=estimate,
                   color=as.factor(type))) + 
  geom_point(size=5) + 
  scale_color_discrete(name="") + 
  theme_classic() + 
  theme(legend.position="bottom")

#-----------------------------------------------------------------------------
# imputation
#
# Note: this uses custom code (mainly to make it be more compatible
# with other code in the chapter), but `didimputation` and `did2s`
# are likely to be better options in actual applications.
#-----------------------------------------------------------------------------

# load code to run imputation
source("imputation.R")

# for imputation with `data4`, need to modify it slightly
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
ggpte(imp2)

#-----------------------------------------------------------------------------
# Callaway and Sant'Anna build-the-trend
#-----------------------------------------------------------------------------

source("build_the_trend.R")

btt2 <- pte(yname="lemp",
            tname="year",
            idname="countyreal",
            gname="G",
            data=data2,
            setup_pte_fun=setup_pte_basic,
            subset_fun=keep_all_subset,
            attgt_fun=build_the_trend_attgt,
            xformla=~1,
            biters=100,
            cl=1)

summary(btt2)
ggpte(btt2)


#-----------------------------------------------------------------------------
# event study regression
#-----------------------------------------------------------------------------
# which data to use in this section
this_data <- data2

# add event time variable
this_data$e <- ifelse(this_data$G==0, 0, this_data$year - this_data$G)

# run event study regression
twfe_es <- feols(lemp ~ i(e, ever_treated, ref=-1) | countyreal + year,
                 data=this_data)

# summary and plot
summary(twfe_es)
iplot(twfe_es)

# replicate event study using more primitive calculations
Dit <- i(this_data$e, this_data$ever_treated, ref=-1)

ddotDit <- demean(Dit, this_data[,c("countyreal", "year")])

Yit <- this_data$lemp

bete <- solve(t(ddotDit)%*%Dit) %*% t(ddotDit) %*% Yit
bete

# load code for computing regression weights
source("event_study_reg_weights.R")

# should replicate bete
es_weights_reg()

# get weights for particular groups across event times
e <- 5
g <- 2002
es_weights(g,e,weightonly=TRUE)

# compare to regression from Sun & Abraham
es_weight_outcome <- 1*( (this_data$G == g) & (g+e == this_data$year))
es_weight_e <- solve(t(ddotDit)%*%Dit) %*% t(ddotDit) %*% es_weight_outcome
es_weight_e


# compute underlying weights for a particular event time
eval.idx <- which(e.seq==e)
all_eg <- expand.grid(e.seq, glist)
colnames(all_eg) <- c("e", "g")
e_weight <- sapply(1:nrow(eg), function(i) es_weights(g=all_eg$g[i], e=all_eg$e[i], weightonly=TRUE)[eval.idx])
e_att <- sapply(1:nrow(eg), function(i) es_weights(g=all_eg$g[i], e=all_eg$e[i], weightonly=FALSE)[eval.idx])
all_eg <- cbind.data.frame(all_eg, weight=e_weight, att=e_att) 
e_weight_this_e <- all_eg[all_eg$e==e,]
e_weight_other_e <- all_eg[all_eg$e!=e,]

# check that weights on "correct" e sum to 1
sum(e_weight_this_e$weight)

# check that weights on "incorrect" e sum to 0
sum(e_weight_other_e$weight)

# find large weights on "incorrect" e
e_weight_other_e[abs(e_weight_other_e$weight) > .1,]




# event study plot reported in chapter that includes Callaway &
# Sant'Anna, event study regression using `data` and `data2`
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


# twfe event study with covariates
twfeX_es <- feols(lemp ~ i(e, ref=-1) + lpop + lavg_pay  | countyreal + year + region^year,
                 data=this_data)
summary(twfeX_es)
iplot(twfeX_es)


#-----------------------------------------------------------------------------
# honest did
#-----------------------------------------------------------------------------

# this is the relative magnitudes version of RR where Mbar
# is the fraction of deviation relative to max observed in
# pre-treatment periods.  this is what is reported in chapter.
#
# change `e` to get sensitivity analysis at different length
# of exposure
hd_rm <- honest_did(es=cs_dyn, type="relative_magnitude", Mbarvec=c(0.5,1,1.5,2), gridPoints=100, grid.lb=-.25, grid.ub=.25, e=0)

createSensitivityPlot_relativeMagnitudes(hd_rm$robust_ci,
                                         hd_rm$orig_ci) +
  theme_classic() +
  scale_color_discrete("", labels=c("Original", "RR Bounds")) +
  theme(legend.position="bottom")


# alternatively, can replace maximum violation of parallel trends
# with maximum violation of linear trends, which is here...
# this is not reported in chapter
hd_lt <- honest_did(es=cs_dyn, type="relative_magnitude", Mbarvec=c(0.5,1,1.5,2), bound="deviation from linear trend", gridPoints=100, grid.lb=-.25, grid.ub=.25, e=0)
createSensitivityPlot_relativeMagnitudes(hd_lt$robust_ci,
                                         hd_lt$orig_ci) + theme_classic()

# alternatively, RR also discuss smoothness restrictions; I don't
# think that this is related to violations in pre-treatment periods
# but rather how much linear trends can be violated, here is code...
# this is not reported in chapter
hd_smooth <- honest_did(es=cs_dyn, type="smoothness",e=0)
createSensitivityPlot(hd_smooth$robust_ci,
                      hd_smooth$orig_ci)

