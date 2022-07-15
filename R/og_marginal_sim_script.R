#' simulation script
#'
#' This function was used to create the simulation results.
#' @export

# library(CFsurvival)
# library(dplyr)
# rm(list=ls())

#####
simul_DML_causal_misspec_marginal_og <- function(N = 500, # n = 250,
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
                                              beta_interactions = NULL,
                                              misspe = ""
){

    library(dplyr)
    set.seed(s)
    expit <- function(x,beta) 1/(1 + exp(-x %*% beta))
    covs <- ncol(X)

    # simulate RCT attribution #

    # if("sampling" %in% misspe)
    # gRs <- expit(exp(-sqrt(X)/2)-1/2, beta_R)
    # else

    gRs <- expit(X, beta_R)
    RCT <- rbinom(n = N, size = 1, prob = gRs)
    # sum(RCT)/length(RCT)

    # treatment assignment #
    # if("treatment" %in% misspe)
    #     g0s <- ifelse(RCT == 1,1/2 ,expit(exp(-sqrt(X)/2)-1/2, beta_A))
    # else
    g0s <- expit(X, beta_A)
    A <- rbinom(n=N, size=1, prob=g0s)
    # sum(A)/length(A)
    # sum(A[RCT==1])/length(A[RCT==1])
    # sum(A[RCT==0])/length(A[RCT==0])

    # time #
    lambdaT=lambda
    rhoT=rho
    X_temp <- X
    # if("survival" %in% misspe){
    #     dat <- simulWeib(N=N, lambda=lambdaT, rho=rhoT, beta=beta_T, rateC=rateC, X = cbind(exp(-sqrt(X)/2)-1/2,A),s = s+1, interactions = beta_interactions, A = A, covs_int = exp(-X_temp/2)-1/2)
    #     t0 <- simulWeib(N=N, lambda=lambdaT, rho=rhoT, beta=beta_T, rateC=rateC, X = cbind(exp(-sqrt(X)/2)-1/2,0),s = s+1, interactions = beta_interactions, A = 0, covs_int = exp(-X_temp/2)-1/2)
    #     t1 <- simulWeib(N=N, lambda=lambdaT, rho=rhoT, beta=beta_T, rateC=rateC, X = cbind(exp(-sqrt(X)/2)-1/2,1),s = s+1, interactions = beta_interactions, A = 1, covs_int = exp(-X_temp/2)-1/2)
    # }
    # else{
    dat <- simulWeib(N=N, lambda=lambdaT, rho=rhoT, beta=beta_T, rateC=rateC, X = cbind(X,A),s = s+1, interactions = beta_interactions, A = A, covs_int = X_temp)
    t0 <- simulWeib(N=N, lambda=lambdaT, rho=rhoT, beta=beta_T, rateC=rateC, X = cbind(X,0),s = s+1, interactions = beta_interactions, A = 0, covs_int = X_temp)
    t1 <- simulWeib(N=N, lambda=lambdaT, rho=rhoT, beta=beta_T, rateC=rateC, X = cbind(X,1),s = s+1, interactions = beta_interactions, A = 1, covs_int = X_temp)
    # }
    # summary(dat$time)
    # sum(dat$status)
    # sum(dat$status[RCT == 1])

    fit.times = seq(0, 10, by=0.1) # seq(0, as.numeric(quantile(dat$time[RCT == 1], 0.995)), by=0.1)

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
    confounders <- X[RCT == 1,] # X_mod[RCT == 1,]
    W_c <- X[RCT == 0,] # X_mod[RCT == 0,]

    ################
    # X_mod <- cbind(log(sqrt(X[,1])/2)-1/2, exp(X[,2] - X[,1]))
    # # X_mod <- matrix(rnorm(n = length(RCT)*ncol(X), mean = 4, sd = 4), nrow = length(RCT), ncol = ncol(X)) #cbind(rnorm(n = length(RCT),mean = 2, sd = 3),rnorm(n = length(RCT),mean = 5, sd = 1))
    # mod_conf <- X_mod[RCT==1,]
    # mod_W_c <- X_mod[RCT ==0,]


    fit_Cohort <- CFsurvival(time = obs.time, event = obs.event, treat = rx,
                             confounders =  rep(1,length(obs.time)),# confounders,
                             contrasts = NULL,
                             verbose = FALSE, fit.times = fit.times,
                             fit.treat = c(0,1),
                             nuisance.options = list(
                                 prop.RCT.SL.library = prop.RCT.SL.library,
                                 prop.SL.library = prop.SL.library,
                                 event.SL.library = event.SL.library,
                                 cens.SL.library = cens.SL.library, verbose = F)
    )



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

    return(list(
        # "fit" = fit,
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

