#' simulation script
#'
#' This function was used to create the simulation results.
#' @export



# survival function generation #
simulWeib <- function(N, lambda, rho, beta, rateC, X,s, interactions = NULL, A = NULL, covs_int = NULL)
{
    set.seed(s)
    if(is.null(interactions))
        interactions = rep(0,ncol(X)-1)
    if(is.null(covs_int))
        covs_int <- 0
    # Weibull latent event times
    v <- runif(n=N)
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

#' simulation script
#'
#' This function was used to create the simulation results.
#' @export
#' @import
#' dplyr
#' survival

simul_DML_causal <- function(N = 500, # n = 250,
                             X,
                             beta_R, # = runif(n=covs,-1,1),
                             beta_A, # = runif(n=covs,-1,1),
                             beta_T, # = c(runif(n=covs,-1,1),-1),
                             beta_C,
                             prop.RCT.SL.library = c("SL.mean", "SL.bayesglm"),
                             prop.SL.library = c("SL.mean", "SL.bayesglm"),
                             event.SL.library = c("survSL.coxph", "survSL.weibreg", "survSL.expreg"),
                             cens.SL.library = c("survSL.coxph", "survSL.weibreg", "survSL.expreg"),
                             lambda = 0.1, rho = 2, rateC = 0.05,s = 210793, quants = c(0.75,0.5,0.25),
                             # parameters for interactions #
                             misspe = "",
                             beta_interactions = NULL, run_type = "correct", surv_type = "exp", cov_int = NULL,
                             betaT = 1, lambdaT = 5, betaC = 2, lambdaC = 10, fit.times = seq(0, 10, by=0.1)
){

    # library(dplyr)
    set.seed(s)
    expit <- function(x,beta,int = 0) 1/(1 + exp(-x %*% beta)+int)
    covs <- ncol(X)

    # simulate RCT attribution #

    # if("sampling" %in% misspe)
    # gRs <- expit(exp(-sqrt(X)/2)-1/2, beta_R)
    # else

    gRs <- expit(X, beta_R,1)
    RCT <- rbinom(n = N, size = 1, prob = gRs)
    # sum(RCT)/length(RCT)

    # treatment assignment #
    # if("treatment" %in% misspe)
    #     g0s <- ifelse(RCT == 1,1/2 ,expit(exp(-sqrt(X)/2)-1/2, beta_A))
    # else
    g0s <- ifelse(RCT == 1,1/2 ,expit(X, beta_A))
    A <- rbinom(n=N, size=1, prob=g0s)
    # sum(A)/length(A)
    # sum(A[RCT==1])/length(A[RCT==1])
    # sum(A[RCT==0])/length(A[RCT==0])

    # time #
    if(surv_type == "exp"){
        lambdaT=lambda
        rhoT=rho
        X_temp <- X
        dat <- simulWeib(N=N, lambda=lambdaT, rho=rhoT, beta=beta_T, rateC=rateC, X = cbind(X,A),s = s+1, interactions = beta_interactions, A = A, covs_int = X_temp)
        t0 <- simulWeib(N=N, lambda=lambdaT, rho=rhoT, beta=beta_T, rateC=rateC, X = cbind(X,0),s = s+1, interactions = beta_interactions, A = 0, covs_int = X_temp)
        t1 <- simulWeib(N=N, lambda=lambdaT, rho=rhoT, beta=beta_T, rateC=rateC, X = cbind(X,1),s = s+1, interactions = beta_interactions, A = 1, covs_int = X_temp)
        t0 <- t0$true_time
        t1 <- t1$true_time

        # Make final data #
        obs.time <- dat$time[RCT == 1]
        obs.event <- dat$status[RCT == 1]
        rx <- A[RCT == 1]
        confounders <- X[RCT == 1,]
        W_c <- X[RCT == 0,]
    }

    if(surv_type == "weibull"){
        # if(is.null(beta_interactions))
        #     beta_interactions <- rep(0,ncol(X) +1)
        A*covs_int %*% interactions
        cens.time <- rweibull(N, shape = betaC, scale = lambdaC * exp(-cbind(X,A) %*% beta_C))
        set.seed(s+1)
        event.time <- rweibull(N, shape = betaT, scale = lambdaT * exp(-cbind(X,A) %*% beta_T))
        set.seed(s+1)
        t0 <- rweibull(N, shape = betaT, scale = lambdaT * exp(-cbind(X,0) %*% beta_T))
        set.seed(s+1)
        t1 <- rweibull(N, shape = betaT, scale = lambdaT * exp(-cbind(X,1) %*% beta_T))

        cens.time[cens.time > 15] <- 15
        obs.time <- pmin(event.time, cens.time)
        obs.event <- as.numeric(event.time <= cens.time)

        # Make final data #
        obs.time <- obs.time[RCT == 1]
        obs.event <- obs.event[RCT == 1]
        rx <- A[RCT == 1]
        confounders <- X[RCT == 1,]
        W_c <- X[RCT == 0,]
    }

    # fit.times = seq(0, 10, by=0.1)

    S0 <- function(t){
        sum(t0 > t)/N

    }
    S1 <- function(t){
        sum(t1 > t)/N
    }
    out_S0 <- sapply(X = fit.times, FUN = S0)
    out_S1 <- sapply(X = fit.times, FUN = S1)


    S0_RCT <- function(t){
        sum(t0[RCT == 1] > t)/sum(RCT)
    }
    S1_RCT <- function(t){
        sum(t1[RCT == 1] > t)/sum(RCT)
    }
    out_S0_RCT <- sapply(X = fit.times, FUN = S0_RCT)
    out_S1_RCT <- sapply(X = fit.times, FUN = S1_RCT)


    ################
    # X_mod <- cbind(sqrt(abs(X[,1])/2)-1/2,exp(X[,2] - X[,1])^2, log((X[,3]+ X[,4])^2), sqrt(abs(X[,4]+X[,5])), (X[,3]+X[,4]+X[,5])^2)
    X_mod <- cbind(sqrt(abs(X[,1])/2),exp(X[,2] - X[,1]), log(abs(X[,3]+ X[,4] + 1)), sqrt(abs(X[,4]+X[,5])), (X[,5])^2)

    # X_mod <- cbind(sqrt(abs(X[,1])/2),exp(X[,2] - X[,1]), log(abs(X[,3]+ X[,4])+1), sqrt(abs(X[,4]+X[,5])), (X[,5])^2)

    if(covs > 5)
        X_mod <- cbind(X_mod, exp(X[,6:ncol(X)]))
    mod_conf <- X_mod[RCT==1,]
    mod_W_c <- X_mod[RCT ==0,]

    ### adding in the cohort data ###
    if(run_type == "correct")
        fit <- CFsurvival(time = obs.time, event = obs.event, treat = rx,
                          confounders =  confounders, contrasts = NULL,
                          verbose = FALSE, fit.times = fit.times,
                          fit.treat = c(0,1),
                          nuisance.options = list(
                              prop.RCT.SL.library = prop.RCT.SL.library,
                              prop.SL.library = prop.SL.library,
                              event.SL.library = event.SL.library,
                              cens.SL.library = cens.SL.library, verbose = F),
                          W_c = W_c
        )

    if(run_type == "marginal")
        fit <- CFsurvival(time = obs.time, event = obs.event, treat = rx,
                          confounders =  rep(1,sum(RCT == 1)), contrasts = NULL,
                          verbose = FALSE, fit.times = fit.times,
                          fit.treat = c(0,1),
                          nuisance.options = list(
                              prop.RCT.SL.library = prop.RCT.SL.library,
                              prop.SL.library = prop.SL.library,
                              event.SL.library = event.SL.library,
                              cens.SL.library = cens.SL.library, verbose = F),
                          W_c = rep(1,sum(RCT == 0))
        )

    if(run_type == "incorrect"){
        if(length(misspe) == 1 & misspe == "sampling")
            fit <- CFsurvival(time = obs.time, event = obs.event, treat = rx,
                              confounders =  confounders, contrasts = NULL,
                              verbose = FALSE, fit.times = fit.times,
                              fit.treat = c(0,1),
                              nuisance.options = list(
                                  prop.RCT.SL.library = prop.RCT.SL.library,
                                  prop.SL.library = prop.SL.library,
                                  event.SL.library = event.SL.library,
                                  cens.SL.library = cens.SL.library, verbose = F),
                              W_c = W_c,misspe = misspe, mod_conf = mod_conf,mod_W_c = mod_W_c
            )

        if(length(misspe) == 1 & misspe == "survival")
            fit <- CFsurvival(time = obs.time, event = obs.event, treat = rx,
                              confounders =  confounders, contrasts = NULL,
                              verbose = FALSE, fit.times = fit.times,
                              fit.treat = c(0,1),
                              nuisance.options = list(
                                  prop.RCT.SL.library = prop.RCT.SL.library,
                                  prop.SL.library = prop.SL.library,
                                  event.SL.library = event.SL.library,
                                  cens.SL.library = cens.SL.library, verbose = F),
                              W_c = W_c,misspe = misspe, mod_conf = mod_conf,mod_W_c = mod_W_c
            )

        if(sum(c("sampling", "survival") %in% misspe) == 2)
            fit <- CFsurvival(time = obs.time, event = obs.event, treat = rx,
                              confounders =  mod_conf, contrasts = NULL,
                              verbose = FALSE, fit.times = fit.times,
                              fit.treat = c(0,1),
                              nuisance.options = list(
                                  prop.RCT.SL.library = prop.RCT.SL.library,
                                  prop.SL.library = prop.SL.library,
                                  event.SL.library = event.SL.library,
                                  cens.SL.library = cens.SL.library, verbose = F),
                              W_c = mod_W_c
            )

    }
    fit$surv.df$true.surv <- c(out_S0[1:(nrow(fit$surv.df)/2)], out_S1[1:(nrow(fit$surv.df)/2)])
    ###### metrics #####

    ### Overall ###

    ### bias survival ###

    Delta_0_bias <- function(t){
        temp <- fit$surv.df %>% filter(trt == 0, time <= t)
        sum(abs(temp$surv - temp$true.surv))
    }
    Delta_0_bias_out <- sapply(X = fit.times, FUN = Delta_0_bias)
    Delta_1_bias <- function(t){
        temp <- fit$surv.df %>% filter(trt == 1, time <= t)
        sum(abs(temp$surv - temp$true.surv))
    }
    Delta_1_bias_out <- sapply(X = fit.times, FUN = Delta_1_bias)

    ### MSE survival ###
    Delta_0_MSE <- function(t){
        temp <- fit$surv.df %>% filter(trt == 0, time <= t)
        sum(abs(temp$surv - temp$true.surv)^2)
    }
    Delta_0_MSE_out <- sapply(X = fit.times, FUN = Delta_0_MSE)

    Delta_1_MSE <- function(t){
        temp <- fit$surv.df %>% filter(trt == 1, time <= t)
        sum(abs(temp$surv - temp$true.surv)^2)
    }
    Delta_1_MSE_out <- sapply(X = fit.times, FUN = Delta_1_MSE)


    # ### Quantile ###
    # ### Quantile bias ###
    quantile_bias <- function(quants,time_surv,survival){
        time_quant <- c()
        for(i in 1:length(quants)){
            ind <- which.min(abs(quants[i] - survival))
            time_quant[i] <- time_surv[ind]
        }
        return(time_quant)
    }

    # # no treatment #
    quant_bias_0 <- quantile_bias(quants = quants, time_surv = fit.times, survival = out_S0) -
        quantile_bias(quants = quants, time_surv = fit.times, survival = fit$surv.df$surv[fit$surv.df$trt == 0])
    names(quant_bias_0) <- quants

    # # treatment #
    quant_bias_1 <- quantile_bias(quants = quants, time_surv = fit.times, survival = out_S1) -
        quantile_bias(quants = quants, time_surv = fit.times, survival = fit$surv.df$surv[fit$surv.df$trt == 1])
    names(quant_bias_1) <- quants
    #
    # ### Quantile MSE ###
    quant_MSE_0 <- quant_bias_0^2
    quant_MSE_1 <- quant_bias_1^2


    ####################
    # merge data #
    dat <- cbind(cbind(cbind(dat,X),A),RCT,X_mod)

    return(list(
        "fit" = fit,
        "out_S0" = out_S0,
        "out_S1" = out_S1,
        # "out_S0_RCT" = out_S0_RCT,
        # "out_S1_RCT" = out_S1_RCT,
        "dat" = dat,
        "Delta_0_bias" = Delta_0_bias_out,
        "Delta_1_bias" = Delta_1_bias_out,
        "Delta_0_MSE" = Delta_0_MSE_out,
        "Delta_1_MSE" = Delta_1_MSE_out,
        "quant_bias_0" = quant_bias_0,
        "quant_bias_1" = quant_bias_1,
        "quant_MSE_0" = quant_MSE_0,
        "quant_MSE_1" = quant_MSE_1#,
        # "gRs" = gRs
    ))

}


# rm(list=ls())
# # N = 1000
# # covs <- c(5)
# # lambda = 0.1
# # rho = 2
# # rateC = 0.1
# # run <- 23
# # set.seed(run*covs)
# #
# N = 1000
# vars <- 1
# lambda = 0.1
# rho = 2
# rateC = 0.05
# runs <- 100
# covs <- 5
#
# X <- matrix(rnorm(n = N*covs,0,1), nrow = N, ncol = covs)
# beta_R = rep(0.5, covs)
# beta_A = rep(0.5, covs)
# beta_T = c(rep(0.5, covs),-1)
# beta_C = c(rep(0.5, covs),-1/5)
#
#
# N = N
# X = X
# beta_R = beta_R
# beta_A = beta_A
# beta_T = beta_T
# lambda = lambda
# rho = rho
# rateC = rateC
# s = 21071993 + 1
# beta_interactions = NULL
# surv_type = "exp"
# run_type = "incorrect"
# betaT = 1
# lambdaT = 5
# betaC = 2
# lambdaC = 10
# misspe = c("sampling")
#
# prop.RCT.SL.library = c("SL.glm", "SL.mean","SL.rpart")
# prop.SL.library = c("SL.mean")
# event.SL.library = c( "survSL.coxph", "survSL.weibreg", "survSL.expreg")
# cens.SL.library = c("survSL.coxph", "survSL.weibreg", "survSL.expreg")
#
# ex_simul <- simul_DML_causal(N = N, # n = n,
#                              X = X,
#                              beta_R = beta_R,
#                              beta_A = beta_A,
#                              beta_T = beta_T,
#                              beta_C = beta_C,
#                              prop.RCT.SL.library = prop.RCT.SL.library,
#                              prop.SL.library = prop.SL.library,
#                              event.SL.library = event.SL.library,
#                              cens.SL.library =cens.SL.library,
#                              lambda = lambda, rho = rho, rateC = rateC, s = s,
#                              run_type = "incorrect",
#                              surv_type = "exp",
#                              misspe = misspe
# )
#
# ex_simul$Delta_0_bias[80]
# ex_simul$Delta_1_bias[80]
# summary(ex_simul$fit$nuisance$prop.pred.RCT)


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
#     coord_cartesian(xlim=c(0,10), ylim=c(0,1))



######
# N = 500
# covs <- 1
# set.seed(133)
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
#                              beta_C = beta_C,
#                              prop.RCT.SL.library = c("SL.mean", "SL.glm"),
#                              prop.SL.library = c("SL.mean", "SL.glm"),
#                              event.SL.library = c( "survSL.coxph", "survSL.weibreg"),
#                              cens.SL.library = c("survSL.coxph", "survSL.expreg"),
#                              lambda = lambda, rho = rho, rateC = rateC
#                              # , beta_interactions = c(1,rep(0,9))
# )
#
# sum(ex_simul$dat$RCT)/N
# sum(ex_simul$dat$A)/N
# sum(ex_simul$dat$status)/N


# #
# # save(ex_simul, file = "~/Causal/CFsurvival/results/example.Rdata")
#
# # summary(ex_simul$fit_Cohort$nuisance$prop.pred.RCT)
# # summary(ex_simul$fit_Cohort$nuisance$prop.pred)
# #
# # ex_simul$Delta_0_bias
# # cumsum(ex_simul$Delta_0_bias)
# #
# # ex_simul$Delta_1_bias
# # cumsum(ex_simul$Delta_1_bias)
# #
# # # quantile stuff #
# # ex_simul$quant_bias_0
# # ex_simul$quant_bias_1
# #
# # ex_simul$quant_MSE_0
# # ex_simul$quant_MSE_1
# #
# # summary(ex_simul$fit_Cohort$nuisance$prop.pred.RCT)
#
# #
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
#     coord_cartesian(xlim=c(0,10), ylim=c(0,1))
#
# #
# #
# # fit_Cohort$surv.df$true.surv <- c(out_S1, out_S0)
# # ggplot(fit_Cohort$surv.df) +
# #     geom_line(aes(time, true.surv, group=trt), color='black') +
# #     geom_step(aes(time, surv, color=as.factor(trt), group=trt)) +
# #     geom_step(aes(time, ptwise.lower, color=as.factor(trt), group=trt), linetype=2) +
# #     geom_step(aes(time, ptwise.upper, color=as.factor(trt), group=trt), linetype=2) +
# #     geom_step(aes(time, unif.ew.lower, color=as.factor(trt), group=trt), linetype=3) +
# #     geom_step(aes(time, unif.ew.upper, color=as.factor(trt), group=trt), linetype=3) +
# #     scale_color_discrete("Treatment") +
# #     xlab("Time") +
# #     ylab("Survival") +
# #     coord_cartesian(xlim=c(0,as.numeric(quantile(dat$time[dat$RCT == 1], 0.995))), ylim=c(0,1))
