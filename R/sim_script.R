#' simulation script
#'
#' This function was used to create the simulation results.
#' @export

# library(CFsurvival)
# library(dplyr)
# rm(list=ls())

# survival function generation #
simulWeib <- function(N, lambda, rho, beta, rateC, X,s, interactions = NULL, A = NULL, covs_int = NULL)
{
    set.seed(s)
    if(is.null(interactions))
        interactions = rep(0,ncol(X)-1)
    # Weibull latent event times
    v <- runif(n=N)
    # if(!is.null(dim(covs_int)))
    Tlat <- (- log(v) / (lambda * exp(X %*% beta + A*covs_int %*% interactions)))^(1 / rho)
    # else
    #     Tlat <- (- log(v) / (lambda * exp(X * beta + A*X * interactions )))^(1 / rho)

    # censoring times
    C <- rexp(n=N, rate=rateC)

    # follow-up times and event indicators
    time <- pmin(Tlat, C)
    status <- as.numeric(Tlat <= C)

    # data set
    data.frame(
        time=time,
        status=status,
        true_time = Tlat
    )
}

#####
simul_DML_causal <- function(N = 500, # n = 250,
                             X,
                             beta_R, # = runif(n=covs,-1,1),
                             beta_A, # = runif(n=covs,-1,1),
                             beta_T, # = c(runif(n=covs,-1,1),-1),
                             prop.RCT.SL.library = c("SL.mean", "SL.bayesglm"),
                             prop.SL.library = c("SL.mean", "SL.bayesglm"),
                             event.SL.library = c("survSL.km", "survSL.coxph", "survSL.weibreg", "survSL.expreg"),
                             cens.SL.library = c("survSL.km", "survSL.coxph", "survSL.weibreg", "survSL.expreg"),
                             lambda = 0.1, rho = 2, rateC = 0.05,s = 210793, quants = c(0.75,0.5,0.25),
                             # parameters for interactions #
                             beta_interactions = NULL
                             ){

    library(dplyr)
    set.seed(s)
    expit <- function(x,beta) 1/(1 + exp(-x %*% beta))
    covs <- ncol(X)
    # X <- matrix(runif(n = N*covs,-1,1), nrow = N, ncol = covs)

    # simulate RCT attribution #
    # beta_R <- beta_R
    gRs <- expit(X, beta_R)
    RCT <- rbinom(n = N, size = 1, prob = gRs)

    # treatment assignment #
    # beta_A <- beta_A
    g0s <- expit(X, beta_A)
    A <- rbinom(n=N, size=1, prob=g0s)
    sum(A)

    # time #
    # beta_T <- beta_T
    # beta_interactions <- c(1, rep(0, covs-1))
    lambdaT=0.1
    rhoT=2
    X_temp <- X
    dat <- simulWeib(N=N, lambda=lambdaT, rho=rhoT, beta=beta_T, rateC=rateC, X = cbind(X,A),s = s+1, interactions = beta_interactions, A = A, covs_int = X_temp)
    t0 <- simulWeib(N=N, lambda=lambdaT, rho=rhoT, beta=beta_T, rateC=rateC, X = cbind(X,0),s = s+1, interactions = beta_interactions, A = 0, covs_int = X_temp)
    t1 <- simulWeib(N=N, lambda=lambdaT, rho=rhoT, beta=beta_T, rateC=rateC, X = cbind(X,1),s = s+1, interactions = beta_interactions, A = 1, covs_int = X_temp)
    summary(dat$time)
    sum(dat$status)

    fit.times = seq(0, as.numeric(quantile(dat$time[RCT == 1], 0.995)), by=0.1)

    S0 <- function(t){
        # sum(t0$time > t)/(N - nrow(t0 %>% filter(time < t, status == 0)))
        sum(t0$true_time > t)/N

    }
    S1 <- function(t){
        sum(t1$true_time > t)/N
        # sum(t1$time > t)/(N - nrow(t1 %>% filter(time < t, status == 0)))
    }
    out_S0 <- sapply(X = fit.times, FUN = S0)
    out_S1 <- sapply(X = fit.times, FUN = S1)


    S0_RCT <- function(t){
        # sum(t0$time[RCT == 1] > t)/(sum(RCT) - nrow(t0[RCT == 1,] %>% filter(time < t, status == 0)))
        sum(t0$true_time[RCT == 1] > t)/sum(RCT)
    }
    S1_RCT <- function(t){
        # sum(t1$time[RCT == 1] > t)/(sum(RCT) - nrow(t1[RCT == 1,] %>% filter(time < t, status == 0)))
        sum(t1$true_time[RCT == 1] > t)/sum(RCT)
    }
    out_S0_RCT <- sapply(X = fit.times, FUN = S0_RCT)
    out_S1_RCT <- sapply(X = fit.times, FUN = S1_RCT)


    # Make final data #
    obs.time <- dat$time[RCT == 1]
    obs.event <- dat$status[RCT == 1]
    rx <- A[RCT == 1]
    confounders <- X[RCT == 1,]
    W_c <- X[RCT == 0,]

    # run model #
    fit <- CFsurvival(time = obs.time, event = obs.event, treat = rx,
                      confounders =  confounders, contrasts = NULL,
                      verbose = TRUE, fit.times = fit.times,
                      fit.treat = c(0,1),
                      nuisance.options = list(
                          prop.RCT.SL.library = prop.RCT.SL.library,
                          prop.SL.library = prop.SL.library,
                          event.SL.library = event.SL.library,
                          cens.SL.library = cens.SL.library),
                      W_c = NULL)

    ### adding in the cohort data ###
    fit_Cohort <- CFsurvival(time = obs.time, event = obs.event, treat = rx,
                             confounders =  confounders, contrasts = NULL,
                             verbose = TRUE, fit.times = fit.times,
                             fit.treat = c(0,1),
                             nuisance.options = list(
                                 prop.RCT.SL.library = prop.RCT.SL.library,
                                 prop.SL.library = prop.SL.library,
                                 event.SL.library = event.SL.library,
                                 cens.SL.library = cens.SL.library),
                             W_c = W_c)



    ###### metrics #####

    ### Overall ###

    ### bias survival ###
    Delta_0_bias <- function(t){
        temp <- fit_Cohort$surv.df %>% filter(trt == 0)
        sum(t0$time > t)/(N - nrow(t0 %>% filter(time < t, status == 0))) - temp$surv[temp$time == t]
    }
    Delta_0_bias_out <- sapply(X = fit.times, FUN = Delta_0_bias)

    Delta_1_bias <- function(t){
        temp <- fit_Cohort$surv.df %>% filter(trt == 1)
        sum(t1$time > t)/(N - nrow(t1 %>% filter(time < t, status == 0))) - temp$surv[temp$time == t]
    }
    Delta_1_bias_out <- sapply(X = fit.times, FUN = Delta_1_bias)

    ### MSE survival ###
    Delta_0_MSE <- function(t){
        temp <- fit_Cohort$surv.df %>% filter(trt == 0)
        (sum(t0$time > t)/(N - nrow(t0 %>% filter(time < t, status == 0))) - temp$surv[temp$time == t])^2
    }
    Delta_0_MSE_out <- sapply(X = fit.times, FUN = Delta_0_MSE)

    Delta_1_MSE <- function(t){
        temp <- fit_Cohort$surv.df %>% filter(trt == 1)
        (sum(t1$time > t)/(N - nrow(t1 %>% filter(time < t, status == 0))) - temp$surv[temp$time == t])^2
    }
    Delta_1_MSE_out <- sapply(X = fit.times, FUN = Delta_1_MSE)


    ### Quantile ###
    ### Quantile bias ###
    quantile_bias <- function(quants,time_surv,survival){
        time_quant <- c()
        for(i in 1:length(quants)){
            ind <- which.min(abs(quants[i] - survival))
            time_quant[i] <- time_surv[ind]
        }
        return(time_quant)
    }

    # no treatment #
    quant_bias_0 <- quantile_bias(quants = quants, time_surv = fit.times, survival = out_S0) -
        quantile_bias(quants = quants, time_surv = fit.times, survival = fit_Cohort$surv.df$surv[fit_Cohort$surv.df$trt == 0])
    names(quant_bias_0) <- quants

    # treatment #
    quant_bias_1 <- quantile_bias(quants = quants, time_surv = fit.times, survival = out_S1) -
        quantile_bias(quants = quants, time_surv = fit.times, survival = fit_Cohort$surv.df$surv[fit_Cohort$surv.df$trt == 1])
    names(quant_bias_1) <- quants

    ### Quantile MSE ###
    quant_MSE_0 <- quant_bias_0^2
    quant_MSE_1 <- quant_bias_1^2


    ####################
    # merge data #
    dat <- cbind(cbind(cbind(dat,X),A),RCT)

    return(list("fit" = fit,
                "fit_Cohort" = fit_Cohort,
                "out_S0" = out_S0,
                "out_S1" = out_S1,
                "out_S0_RCT" = out_S0_RCT,
                "out_S1_RCT" = out_S1_RCT,
                "dat" = dat,
                "Delta_0_bias" = Delta_0_bias_out,
                "Delta_1_bias" = Delta_1_bias_out,
                "Delta_0_MSE" = Delta_0_MSE_out,
                "Delta_1_MSE" = Delta_1_MSE_out,
                "quant_bias_0" = quant_bias_0,
                "quant_bias_1" = quant_bias_1,
                "quant_MSE_0" = quant_MSE_0,
                "quant_MSE_1" = quant_MSE_1,
                "gRs" = gRs
    ))

}



# # ######
# N = 1000
# covs <- 10
# set.seed(1)
# X <- matrix(runif(n = N*covs,-1,1), nrow = N, ncol = covs)
# beta_R = runif(n=covs,-1,1)
# beta_A = runif(n=covs,-1,1)
# beta_T = c(runif(n=covs,-1,1),-1)
#
# N = N
# covs = covs
# beta_R = beta_R
# beta_A = beta_A
# beta_T = beta_T
# lambda = 0.1
# rho = 2
# rateC = 0.05
#
# ex_simul <- simul_DML_causal(N = N, # n = n,
#                              X = X,
#                              beta_R = beta_R,
#                              beta_A = beta_A,
#                              beta_T = beta_T,
#                              prop.RCT.SL.library = c("SL.mean", "SL.glm"),
#                              prop.SL.library = c("SL.mean", "SL.glm"),
#                              event.SL.library = c( "survSL.coxph", "survSL.weibreg"),
#                              cens.SL.library = c("survSL.coxph", "survSL.expreg"),
#                              lambda = lambda, rho = rho, rateC = rateC#, beta_interactions = c(1,rep(0,9))
#                              )

# sum(ex_simul$dat$RCT)/N
# sum(ex_simul$dat$A)/N
# sum(ex_simul$dat$status)/N
#
# summary(ex_simul$fit_Cohort$nuisance$prop.pred.RCT)
# summary(ex_simul$fit_Cohort$nuisance$prop.pred)
#
# ex_simul$Delta_0_bias
# cumsum(ex_simul$Delta_0_bias)
#
# ex_simul$Delta_1_bias
# cumsum(ex_simul$Delta_1_bias)
#
# # quantile stuff #
# ex_simul$quant_bias_0
# ex_simul$quant_bias_1
#
# ex_simul$quant_MSE_0
# ex_simul$quant_MSE_1
#
# summary(ex_simul$fit_Cohort$nuisance$prop.pred.RCT)
#
#
# library(ggplot2)
# # First plot the survival curves + conf intervals + conf bands
# # fit$surv.df$true.surv <- c(S1(c(0, fit$fit.times)), S0(c(0, fit$fit.times)))
# ex_simul$fit$surv.df$true.surv <- c(ex_simul$out_S1, ex_simul$out_S0)
# ggplot(ex_simul$fit$surv.df) +
#     geom_line(aes(time, true.surv, group=trt), color='black') +
#     geom_step(aes(time, surv, color=as.factor(trt), group=trt)) +
#     geom_step(aes(time, ptwise.lower, color=as.factor(trt), group=trt), linetype=2) +
#     geom_step(aes(time, ptwise.upper, color=as.factor(trt), group=trt), linetype=2) +
#     geom_step(aes(time, unif.ew.lower, color=as.factor(trt), group=trt), linetype=3) +
#     geom_step(aes(time, unif.ew.upper, color=as.factor(trt), group=trt), linetype=3) +
#     scale_color_discrete("Treatment") +
#     xlab("Time") +
#     ylab("Survival") +
#     coord_cartesian(xlim=c(0,as.numeric(quantile(ex_simul$dat$time[ex_simul$dat$RCT == 1], 0.995))), ylim=c(0,1))
#
#
#
# ex_simul$fit_Cohort$surv.df$true.surv <- c(ex_simul$out_S1, ex_simul$out_S0)
# ggplot(ex_simul$fit_Cohort$surv.df) +
#     geom_line(aes(time, true.surv, group=trt), color='black') +
#     geom_step(aes(time, surv, color=as.factor(trt), group=trt)) +
#     geom_step(aes(time, ptwise.lower, color=as.factor(trt), group=trt), linetype=2) +
#     geom_step(aes(time, ptwise.upper, color=as.factor(trt), group=trt), linetype=2) +
#     geom_step(aes(time, unif.ew.lower, color=as.factor(trt), group=trt), linetype=3) +
#     geom_step(aes(time, unif.ew.upper, color=as.factor(trt), group=trt), linetype=3) +
#     scale_color_discrete("Treatment") +
#     xlab("Time") +
#     ylab("Survival") +
#     coord_cartesian(xlim=c(0,as.numeric(quantile(ex_simul$dat$time[ex_simul$dat$RCT == 1], 0.995))), ylim=c(0,1))
