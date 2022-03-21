#-----------------------------------------------------------------------------
#
# event_study_reg_weights.R
#
# * This is code/functions for computing weights on underlying
#   an event study regression.
#
# * In order to run this code, one should first run a Callaway &
#   Sant'Anna att_gt function with a universal base period,
#   w/o covariates and using the "nevertreated" comparison group,
#   and have saved this in a variable called `cs_res`  The data
#   should also be saved in a variable called `this_data`.
#
#-----------------------------------------------------------------------------

# recover info about time periods
tlist <- unique(this_data$year)
nT <- length(tlist)

# recover unique event times
e.seq <- sort(unique(this_data$e[this_data$e != -1]))

# number of units
n <- length(unique(this_data$id))

# list of groups excluding never treated 
glist <- sort(unique(this_data$G))[-1]

#' Function to recover group-time average treatment effects
#' from a CS att_gt call
#'
#' @param g group
#' @param e event time
#'
#' @return the value of ATT(g,g+e)
get_attgt <- function(g,e) {
  attgt_idx <- cs_res$group == g & cs_res$t == g+e
  if( !(any(attgt_idx)) ) return(0)
  attgt <- cs_res$att[attgt_idx]
  attgt
}

#' Function to compute weights given by relative sizes of groups
#'
#' @param g group
#'
#' @return P(G=g)
p_g <- function(g) {
  mean(this_data$G==g)
}

#' Function to compute first part of particular elements of h_e
#'
#' @param e event time
#' @param g group
#' @param t time period
#'
#' @return first component of weights
he <- function(e,g,t) {
  ( 1*(g+e==t) - 1*( (g+e) %in% tlist )/nT) * (g != 0)
}

#' Function to compute vector of h_e
#'
#' @param g group
#' @param t time period
#'
#' @return vector of h_e
hinner <- function(g,t) {
  sapply(e.seq, he, g=g,t=t)
}

#' Function to compute h(g,t) from chapter; this applies `hinner`
#' and then de-means it
#'
#' @param g group
#' @param t time period
#'
#' @return vector of demeaned h_e
h <- function(g,t) {
  hinner(g,t) - apply(sapply(glist, function(gg) hinner(gg,t)*p_g(gg)), 1, sum)
}

#' Function to compute event study regression weights or event
#' study regression components
#'
#' @param g group
#' @param e event time
#' @param weightonly whether or not to return only the weights or
#'  to return the weights multiplied times the corresponding values
#'  of ATT(g,g+e)
#'
#' @return weights for group at particular event time or contribution
#'  of a particular ATT(g,g+e)
es_weights <- function(g,e,weightonly=FALSE) {
  part1 <- solve(t(ddotDit)%*%Dit/n)
  part2 <- as.matrix(h(g,g+e))*( (g+e) %in% tlist )*(g!=0)*p_g(g)
  if (!weightonly) part2 <- part2*get_attgt(g,e)
  part1 %*% part2
}

#' Function to recover event study regression parameters from
#' underlying components
#'
#' @return event study regression parameters
es_weights_reg <- function() {
  out <- matrix(data=0, nrow=length(e.seq))
  for (e in e.seq) {
    for (g in glist) {
      out <- out + es_weights(g,e)
    }
  }
  out
}

