#' simulation script coxph
#'
#' This function was used to create the simulation results.
#' @export
simul_SL_causal <- function(N = 500, # n = 250,
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
                            beta_interactions = NULL, run_type = NULL, surv_type = "exp",
                            betaT = 1, lambdaT = 5, betaC = 2, lambdaC = 10, fit.times = seq(0, 10, by=1),
                            V = 10
){

    library(dplyr)
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
        # A*covs_int %*% interactions
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
    # X_mod <- cbind(sqrt(abs(X[,1])/2),exp(X[,2] - X[,1]), log(abs(X[,3]+ X[,4])+1), sqrt(abs(X[,4]+X[,5])), (X[,5])^2)
    # X_mod <- cbind(sqrt(abs(X[,1])/2),exp(X[,2] - X[,1]), abs(X[,3]+ X[,4]), sqrt(abs(X[,4]+X[,5])), (X[,5])^2)

    # X_mod <- cbind(sqrt(abs(X[,1])/2),exp(X[,2] - X[,1]), log(abs(X[,3]+ X[,4] + 1)), exp(sqrt(abs(X[,4]+X[,5]))), (X[,5])^2)
    # X_mod <- cbind(sqrt(abs(X[,1])/2),exp(X[,2] - X[,1]), abs(X[,3]+ X[,4]), sqrt(abs(X[,4]+X[,5])), (X[,5])^2)
    X_mod <- cbind(sqrt(abs(X[,1])/2),exp(X[,2] - X[,1]), log(abs(X[,3]+ X[,4] + 1)), sqrt(abs(X[,4]+X[,5])), (X[,5])^2)

    mod_conf <- X_mod[RCT==1,]
    mod_W_c <- X_mod[RCT ==0,]

    ######################################################
    ########## MAKE FOLDS SIMILAR TO THE CF ##############
    ######################################################

    nuis <- do.call("CFsurvival.nuisance.options", list(
        prop.RCT.SL.library = prop.RCT.SL.library,
        prop.SL.library = prop.SL.library,
        event.SL.library = event.SL.library,
        cens.SL.library = cens.SL.library, verbose = F))
    fit.treat <- c(0,1)
    ### marginal --> just a regular cox model ###
    if(run_type == "correct"){

        fit <- .estimate.conditional.survival(Y=obs.time,
                                              Delta=obs.event,
                                              A=rx,
                                              W=as.data.frame(confounders),
                                              newW=as.data.frame(X),
                                              event.SL.library=nuis$event.SL.library,
                                              fit.times=fit.times,
                                              fit.treat=fit.treat,
                                              cens.SL.library=nuis$cens.SL.library,
                                              survSL.control=nuis$survSL.control,
                                              survSL.cvControl = nuis$survSL.cvControl,
                                              cens.trunc=nuis$cens.trunc,
                                              save.fit = nuis$save.nuis.fits,
                                              verbose = nuis$verbose)

        # estimates for events under treat = 0 #
        S.hats.0 <- fit$event.pred.0
        # estimates for events under treat = 1 #
        S.hats.1 <- fit$event.pred.1

        final_pred_0 <- as_tibble(cbind(fit.times, apply(S.hats.0,2, mean),out_S0))
        colnames(final_pred_0) <- c('time','surv','true.surv')

        final_pred_1 <- as_tibble(cbind(fit.times, apply(S.hats.1,2, mean),out_S1))
        colnames(final_pred_1) <- c('time','surv','true.surv')

    }
    if(run_type == "marginal"){

        W <- as.data.frame(rep(1, length(rx)))
        newW <- as.data.frame(rep(1, length(A)))
        colnames(W) <- colnames(newW) <- "X"
        fit <- .estimate.conditional.survival(Y=obs.time,
                                              Delta=obs.event,
                                              A=rx,
                                              W=W,
                                              newW=newW,
                                              event.SL.library=nuis$event.SL.library,
                                              fit.times=fit.times,
                                              fit.treat=fit.treat,
                                              cens.SL.library=nuis$cens.SL.library,
                                              survSL.control=nuis$survSL.control,
                                              survSL.cvControl = nuis$survSL.cvControl,
                                              cens.trunc=nuis$cens.trunc,
                                              save.fit = nuis$save.nuis.fits,
                                              verbose = nuis$verbose)

        # estimates for events under treat = 0 #
        S.hats.0 <- fit$event.pred.0
        # estimates for events under treat = 1 #
        S.hats.1 <- fit$event.pred.1

        final_pred_0 <- as_tibble(cbind(fit.times, apply(S.hats.0,2, mean),out_S0))
        colnames(final_pred_0) <- c('time','surv','true.surv')

        final_pred_1 <- as_tibble(cbind(fit.times, apply(S.hats.1,2, mean),out_S1))
        colnames(final_pred_1) <- c('time','surv','true.surv')

    }

    if(run_type == "incorrect"){

        W <- as.data.frame(mod_conf)
        newW <- as.data.frame(X_mod)
        # colnames(W) <- colnames(newW) <- "X"
        fit <- .estimate.conditional.survival(Y=obs.time,
                                              Delta=obs.event,
                                              A=rx,
                                              W=W,
                                              newW=newW,
                                              event.SL.library=nuis$event.SL.library,
                                              fit.times=fit.times,
                                              fit.treat=fit.treat,
                                              cens.SL.library=nuis$cens.SL.library,
                                              survSL.control=nuis$survSL.control,
                                              survSL.cvControl = nuis$survSL.cvControl,
                                              cens.trunc=nuis$cens.trunc,
                                              save.fit = nuis$save.nuis.fits,
                                              verbose = nuis$verbose)

        # estimates for events under treat = 0 #
        S.hats.0 <- fit$event.pred.0
        # estimates for events under treat = 1 #
        S.hats.1 <- fit$event.pred.1

        final_pred_0 <- as_tibble(cbind(fit.times, apply(S.hats.0,2, mean),out_S0))
        colnames(final_pred_0) <- c('time','surv','true.surv')

        final_pred_1 <- as_tibble(cbind(fit.times, apply(S.hats.1,2, mean),out_S1))
        colnames(final_pred_1) <- c('time','surv','true.surv')


    }

    ###### metrics #####

    ### Overall ###

    ### bias survival ###

    Delta_0_bias <- function(t){
        temp <- final_pred_0 %>% filter(time <= t)
        sum(abs(temp$surv - temp$true.surv))
    }
    Delta_0_bias_out <- sapply(X = fit.times, FUN = Delta_0_bias)
    Delta_1_bias <- function(t){
        temp <- final_pred_1 %>% filter(time <= t)
        sum(abs(temp$surv - temp$true.surv))
    }
    Delta_1_bias_out <- sapply(X = fit.times, FUN = Delta_1_bias)

    ### MSE survival ###
    Delta_0_MSE <- function(t){
        temp <- final_pred_0 %>% filter(time <= t)
        sum(abs(temp$surv - temp$true.surv)^2)
    }
    Delta_0_MSE_out <- sapply(X = fit.times, FUN = Delta_0_MSE)

    Delta_1_MSE <- function(t){
        temp <- final_pred_1 %>% filter(time <= t)
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
        quantile_bias(quants = quants, time_surv = fit.times, survival = final_pred_0$surv) #fit$surv.df$surv[fit$surv.df$trt == 0]
    names(quant_bias_0) <- quants

    # # treatment #
    quant_bias_1 <- quantile_bias(quants = quants, time_surv = fit.times, survival = out_S1) -
        quantile_bias(quants = quants, time_surv = fit.times, survival = final_pred_1$surv)
    names(quant_bias_1) <- quants
    #
    # ### Quantile MSE ###
    quant_MSE_0 <- quant_bias_0^2
    quant_MSE_1 <- quant_bias_1^2


    ####################
    # merge data #
    if(surv_type == "exp")
        dat <- cbind(cbind(cbind(dat,X),A),RCT,X_mod)
    if(surv_type =="weibull")
        dat <- cbind(cbind(cbind(event.time, cens.time,X),A),RCT,X_mod)




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


# library(CFsurvival)
# rm(list=ls())
# library(survival)
# N = 1000
# vars <- c(5)
# lambda = 0.1
# rho = 2
# rateC = 0.1
# runs <- 250
# covs <- 5
# run = 7
# set.seed(run*covs)
#
# X <- matrix(rnorm(n = N*covs,0,1), nrow = N, ncol = covs)
# beta_R = rep(0.5, covs)
# beta_A = rep(0.5, covs)
# beta_T = c(rep(0.5, covs),-1)
# beta_C = c(rep(0.5, covs),-1/5)
# # s = 21072028
# N = N
# X = X
# beta_R = beta_R
# beta_A = beta_A
# beta_T = beta_T
# beta_C = beta_C
# prop.RCT.SL.library = c("SL.glm", "SL.glmnet", "SL.rpart")
# prop.SL.library = c("SL.mean")
# event.SL.library = c( "survSL.coxph", "survSL.weibreg", "survSL.expreg","survSL.rfsrc","survSL.gam")
# cens.SL.library = c("survSL.coxph", "survSL.weibreg", "survSL.expreg","survSL.rfsrc","survSL.gam")
# # [1] "survSL.coxph"     "survSL.expreg"    "survSL.gam"       "survSL.km"        "survSL.loglogreg"
# # [6] "survSL.pchreg"    "survSL.pchSL"     "survSL.rfsrc"     "survSL.template"  "survSL.weibreg"
# lambda = lambda
# rho = rho
# rateC = rateC
# s = 21071993 + run*covs
# run_type = "correct"
# surv_type = "exp"
# fit.times = seq(0,10, by = 1)
# V = 10
# beta_interactions = NULL
#
#
# source("./R/estimate_nuisances.R")
#
# CFsurvival.nuisance.options <- function(cross.fit = TRUE, V = ifelse(cross.fit, 10, 1), folds = NULL, eval.times = NULL,
#                                         event.SL.library = lapply(c("survSL.km", "survSL.coxph", "survSL.expreg", "survSL.weibreg", "survSL.loglogreg", "survSL.gam", "survSL.rfsrc"), function(alg) c(alg, "survscreen.glmnet", "survscreen.marg", "All") ),  event.pred.0 = NULL, event.pred.1 = NULL, event.pred.0_c = NULL, event.pred.1_c = NULL,
#                                         cens.SL.library = lapply(c("survSL.km", "survSL.coxph", "survSL.expreg", "survSL.weibreg", "survSL.loglogreg", "survSL.gam", "survSL.rfsrc"), function(alg) c(alg, "survscreen.glmnet", "survscreen.marg", "All") ),  cens.trunc=0, cens.pred.0 = NULL, cens.pred.1 = NULL, cens.pred.0_c = NULL, cens.pred.1_c = NULL,
#                                         survSL.control = list(initWeightAlg = "survSL.rfsrc", verbose=FALSE), survSL.cvControl = list(V = 10), save.nuis.fits = FALSE,
#                                         prop.SL.library = lapply(c("SL.mean", "SL.glm", "SL.gam", "SL.earth", "SL.ranger", "SL.xgboost"), function(alg) c(alg, "screen.glmnet", "screen.corRank", "All") ),
#                                         prop.RCT.SL.library = lapply(c("SL.mean", "SL.glm", "SL.gam", "SL.earth", "SL.ranger", "SL.xgboost"), function(alg) c(alg, "screen.glmnet", "screen.corRank", "All") ),
#                                         prop.trunc=0, prop.pred = NULL,
#                                         verbose=NULL) {
#     list(cross.fit = cross.fit, V = V, folds = folds, eval.times = eval.times,
#          event.SL.library = event.SL.library, event.pred.0 = event.pred.0, event.pred.1 = event.pred.1, event.pred.0_c = event.pred.0_c, event.pred.1_c = event.pred.1_c,
#          cens.SL.library = cens.SL.library, cens.trunc=cens.trunc, cens.pred.0 = cens.pred.0, cens.pred.1 = cens.pred.1, cens.pred.0_c = cens.pred.0_c, cens.pred.1_c = cens.pred.1_c,
#          survSL.control = survSL.control,
#          survSL.cvControl = survSL.cvControl,
#          save.nuis.fits = save.nuis.fits,
#          prop.SL.library = prop.SL.library,  prop.RCT.SL.library = prop.RCT.SL.library, prop.trunc=prop.trunc, prop.pred = prop.pred, verbose=verbose)
# }
#
# ex_simul <- simul_coxph_causal(N = N, # n = n,
#                                X = X,
#                                beta_R = beta_R,
#                                beta_A = beta_A,
#                                beta_T = beta_T,
#                                beta_C = beta_C,
#                                prop.RCT.SL.library = c("SL.glm", "SL.glmnet", "SL.rpart"),
#                                prop.SL.library = c("SL.mean"),
#                                event.SL.library = event.SL.library,
#                                cens.SL.library = cens.SL.library,
#                                lambda = lambda, rho = rho, rateC = rateC, s = 21071993 + run*covs,
#                                run_type = "correct",
#                                surv_type = "exp",
#                                fit.times = seq(0,10, by = 1)
# )
# ex_simul$Delta_0_bias[5]
# ex_simul$Delta_1_bias[5]

# rm(list=ls())
# library(survival)
# N = 1000
# vars <- c(5)
# lambda = 0.1
# rho = 2
# rateC = 0.1
# runs <- 250
# covs <- 5
#
# X <- matrix(rnorm(n = N*covs,0,1), nrow = N, ncol = covs)
# beta_R = rep(0.5, covs)
# beta_A = rep(0.5, covs)
# beta_T = c(rep(0.5, covs),-1)
# beta_C = c(rep(0.5, covs),-1/5)
#
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
# s = 7*5
# # s = sample(1:1000,1) # 21071993 + 12
# beta_interactions = NULL
# surv_type = "exp"
# run_type = "incorrect"
# betaT = 1
# lambdaT = 5
# betaC = 2
# lambdaC = 10
# misspe = c("survival")
# V = 10
#
# prop.RCT.SL.library = c("SL.glm","SL.bayesglm", "SL.glmnet", "SL.rpart")
# prop.SL.library = c("SL.mean")
# event.SL.library = c( "survSL.coxph", "survSL.weibreg", "survSL.expreg")
# cens.SL.library = c("survSL.coxph", "survSL.weibreg", "survSL.expreg")
# fit.times = seq(0,10, 1)
#
#
# ex_simul <- simul_coxph_causal(N = N, # n = n,
#                                X = X,
#                                beta_R = beta_R,
#                                beta_A = beta_A,
#                                beta_T = beta_T,
#                                beta_C = beta_C,
#                                prop.RCT.SL.library = prop.RCT.SL.library,
#                                prop.SL.library = prop.SL.library,
#                                event.SL.library = event.SL.library,
#                                cens.SL.library =cens.SL.library,
#                                lambda = lambda, rho = rho, rateC = rateC, s = s,
#                                run_type = "correct",
#                                surv_type = "exp",
#                                misspe = misspe,
#                                fit.times = fit.times
# )
#
# ex_simul$Delta_0_bias[5]
# ex_simul$Delta_1_bias[5]
#
# run = 7
# ex_simul <- simul_coxph_causal(N = N, # n = n,
#                                X = X,
#                                beta_R = beta_R,
#                                beta_A = beta_A,
#                                beta_T = beta_T,
#                                beta_C = beta_C,
#                                prop.RCT.SL.library = c("SL.glm", "SL.glmnet", "SL.rpart"),
#                                prop.SL.library = c("SL.mean"),
#                                event.SL.library = c( "survSL.coxph", "survSL.weibreg", "survSL.expreg"),
#                                cens.SL.library = c("survSL.coxph", "survSL.weibreg", "survSL.expreg"),
#                                lambda = lambda, rho = rho, rateC = rateC, s = 21071993 + run*covs,
#                                run_type = "correct",
#                                surv_type = "exp",
#                                fit.times = seq(0,10, by = 1)
# )


#
# ex_simul <- simul_coxph_causal(N = N, # n = n,
#                                X = X,
#                                beta_R = beta_R,
#                                beta_A = beta_A,
#                                beta_T = beta_T,
#                                beta_C = beta_C,
#                                prop.RCT.SL.library = prop.RCT.SL.library,
#                                prop.SL.library = prop.SL.library,
#                                event.SL.library = event.SL.library,
#                                cens.SL.library =cens.SL.library,
#                                lambda = lambda, rho = rho, rateC = rateC, s = s,
#                                run_type = "incorrect",
#                                surv_type = "exp",
#                                misspe = misspe,
#                                fit.times = fit.times
# )
#
# ex_simul$Delta_0_bias[5]
# ex_simul$Delta_1_bias[5]
#
#
# ex_simul <- simul_coxph_causal(N = N, # n = n,
#                                X = X,
#                                beta_R = beta_R,
#                                beta_A = beta_A,
#                                beta_T = beta_T,
#                                beta_C = beta_C,
#                                prop.RCT.SL.library = prop.RCT.SL.library,
#                                prop.SL.library = prop.SL.library,
#                                event.SL.library = event.SL.library,
#                                cens.SL.library =cens.SL.library,
#                                lambda = lambda, rho = rho, rateC = rateC, s = s,
#                                run_type = "marginal",
#                                surv_type = "exp",
#                                misspe = misspe,
#                                fit.times = fit.times
# )
#
# ex_simul$Delta_0_bias[5]
# ex_simul$Delta_1_bias[5]

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
