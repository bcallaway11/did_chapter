## ---- message=FALSE-----------------------------------------------------------
library(did)
# for bacondecomp, dev version is much faster
devtools::install_github("evanjflack/bacondecomp")
library(bacondecomp)
library(fixest)
library(modelsummary)
library(ggplot2)
load("mw_data2.RData")


## -----------------------------------------------------------------------------
head(mw_data2)


## -----------------------------------------------------------------------------
# add post-treatment dummy variable
mw_data2$post <- 1*((mw_data2$year >= mw_data2$first.treat) & mw_data2$treat != 0) 

twfe_res <- feols(lemp ~ post | countyreal + year,
                  data=mw_data2,
                  cluster="countyreal")


## -----------------------------------------------------------------------------
modelsummary(twfe_res, gof_omit=".*")


## -----------------------------------------------------------------------------
# run bacon decomposition
bacon_res <- bacon(lemp ~ post, 
                   data=mw_data2,
                   id_var="countyreal",
                   time_var="year")

# confirm same estimate
sum(bacon_res$estimate * bacon_res$weight)

# bacon decomp results
head(bacon_res)


## ----eval=FALSE---------------------------------------------------------------
## # plot bacon decomposition
## ggplot(data=bacon_res,
##        mapping=aes(x=weight,
##                    y=estimate,
##                    color=as.factor(type))) +
##   geom_point(size=5) +
##   scale_color_discrete(name="") +
##   theme_bw() +
##   theme(legend.position="bottom")


## ----echo=FALSE, fig.align='center'-------------------------------------------
ggplot(data=bacon_res, 
       mapping=aes(x=weight,
                   y=estimate,
                   color=as.factor(type))) + 
  geom_point(size=5) + 
  scale_color_discrete(name="") + 
  theme_bw() + 
  theme(legend.position="bottom")


## ---- warning=FALSE-----------------------------------------------------------
# callaway and sant'anna
cs_res <- att_gt(yname="lemp",
                 tname="year",
                 idname="countyreal",
                 gname="first.treat",
                 data=mw_data2)


## ---- fig.align="center"------------------------------------------------------
ggdid(cs_res)


## -----------------------------------------------------------------------------
aggte(cs_res, type="group")


## ---- message=FALSE-----------------------------------------------------------
library(tidyr)
library(dplyr)
# simulation parameters
time.periods <- 20
groups <- c(5,15,time.periods+1)
pg <- c(0.5,0.5,0)
n <- 1000

# generate data (code omitted...)
# load file: sim_data.RDS


## ---- echo=FALSE--------------------------------------------------------------
G <- sample(groups, 
            size=n,
            prob=pg,
            replace=TRUE)
ind_fe <- rnorm(n, mean=G, sd=1)
time_fe <- seq(1,time.periods)
Y0 <- sapply(1:time.periods,
             function(tt) time_fe[tt] + ind_fe + rnorm(n))
data0 <- cbind.data.frame(id=1:n, group=G, Y0)
colnames(data0) <- c("id", "G", paste0("Y0_",1:time.periods))
data <- pivot_longer(data0, 
                     cols=starts_with("Y0"),
                     names_to = "time.period",
                     names_prefix="Y0_",
                     names_transform=list(time.period=as.numeric),
                     values_to = "Y0")
data$post <- 1*(data$time.period >= data$G)
data$e <- data$time.period - data$G

# setup actual outcomes
att0 <- 5
eseq <- 1:time.periods
data$Y <- data$Y0 + data$post*(att0 + (15-data$G)*(data$e>0)*eseq[pmax(data$e,1)])
saveRDS(data, file="sim_data.RDS")


## ---- eval=FALSE--------------------------------------------------------------
## # plot data
## plotdf <- data %>%
##   group_by(G, time.period) %>%
##   summarise(Yobs=mean(Y),
##             Y0=mean(Y0))
## plotdf_obs <- plotdf %>% select(-Y0)
## plotdf_obs$group <- paste0(plotdf$G,"-observed")
## plotdf0 <- plotdf %>% select(-Yobs)
## plotdf0$group <- paste0(plotdf0$G,"-untreated")
## ggplot(data=plotdf,
##        mapping=aes(x=time.period, y=Yobs, color=as.factor(G))) +
##   geom_point() +
##   geom_line() +
##   geom_point(aes(y=Y0)) +
##   geom_line(aes(y=Y0), linetype="dashed") +
##   scale_x_continuous(breaks=seq(2,time.periods,by=2)) +
##   ylab("Y") +
##   theme_bw()


## ---- echo=FALSE, fig.align="center", fig.width=10, fig.height=6, message=FALSE----
plotdf <- data %>%
  group_by(G, time.period) %>%
  summarise(Yobs=mean(Y),
            Y0=mean(Y0))
plotdf_obs <- plotdf %>% select(-Y0)
plotdf_obs$group <- paste0(plotdf$G,"-observed")
plotdf0 <- plotdf %>% select(-Yobs)
plotdf0$group <- paste0(plotdf0$G,"-untreated")
ggplot(data=plotdf,
       mapping=aes(x=time.period, y=Yobs, color=as.factor(G))) + 
  geom_point() +
  geom_line() + 
  geom_point(aes(y=Y0)) + 
  geom_line(aes(y=Y0), linetype="dashed") + 
  scale_x_continuous(breaks=seq(2,time.periods,by=2)) + 
  ylab("Y") + 
  theme_bw()


## -----------------------------------------------------------------------------
# TWFE
twfe_res <- feols(Y ~ post | id + time.period,
            data=data,
            cluster="id")
            
modelsummary(twfe_res, gof_omit=".*")


## -----------------------------------------------------------------------------
cs_res <- att_gt(yname="Y",
                 tname="time.period",
                 idname="id",
                 gname="G",
                 data=data,
                 control_group = "notyettreated")

round(aggte(cs_res, type="group")$overall.att, 3)













# create event time variable
mw_data2$e <- ifelse(mw_data2$treat==1, 
                  mw_data2$year - mw_data2$first.treat,
                  0)

# run event study regression
es_reg <- feols(lemp ~ i(e, ref=-1) | year + countyreal,
                data=mw_data2,
                cluster="countyreal")


D_es_reg <- feols(e ~ i(e, ref=-1) | year + countyreal,
                data=mw_data2,
                cluster="countyreal")


ddotD <- resid(D_es_reg)


es_demean_D <- function(e, data, ename, idname, tname) {
  fixest::demean(X=1*(data[,ename]==e),
                 f=data[,c(idname,tname)])
}

e_seq <- unique(mw_data2$e)
e_seq <- e_seq[e_seq != -1]

ddotD <- sapply(e_seq, function(e) es_demean_D(e, mw_data2, "e", "countyreal", "year"))
n <- length(unique(mw_data2$countyreal))
A <- t(ddotD)%*%ddotD / n
solve(A)

i1 <- rep(0, 9); i1[1] <- 1
t(i1)%*%solve(A)
