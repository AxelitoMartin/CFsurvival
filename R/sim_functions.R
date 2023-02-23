#### TO DO ####
# 1. remove the cross-validation for the other methods --> should simplify the coverage
# 2. switch everything to a difference between S1 and S0 --> S1-S0 (bias, MSE,etc)


### Script for data creation in simulations ###


#' Make data for simulations
#'
#' This function was used to create the simulation datasets.
#' @export
#'


simulWeibDat <- function(N, covs, beta_R, beta_A,
                         lambdaT, lambdaC, beta_T, beta_C,shapeT, shapeC, rho,
                         s, interactions = NULL, covs_int = NULL)
{
    set.seed(s)
    X <- matrix(rnorm(n = N*covs,0,1), nrow = N, ncol = covs)
    colnames(X) <- paste0("X",1:covs)

    if(is.null(interactions))
        interactions = rep(0,ncol(X)-1)
    if(is.null(covs_int))
        covs_int = rep(0, N)

    # generate covariates, sampling and treatment #
    expit <- function(x,beta,int = 0) 1/(1 + exp(-x %*% beta)+int)

    gRs <- expit(X, beta_R,1)
    RCT <- rbinom(n = N, size = 1, prob = gRs)
    # sum(RCT)/length(RCT)

    g0s <- ifelse(RCT == 1,1/2 ,expit(X, beta_A))
    A <- rbinom(n=N, size=1, prob=g0s)


    # generate time data #
    # set.seed(s)
    # cens.time <- rweibull(N, shape = shapeC, scale = lambdaC * exp(-cbind(X,A) %*% beta_C))
    # set.seed(s)
    # event.time <- rweibull(N, shape = shapeT, scale = lambdaT * exp(-cbind(X,A) %*% beta_T))
    # set.seed(s)
    # t0 <- rweibull(N, shape = shapeT, scale = lambdaT * exp(-cbind(X,0) %*% beta_T))
    # set.seed(s)
    # t1 <- rweibull(N, shape = shapeT, scale = lambdaT * exp(-cbind(X,1) %*% beta_T))


    # Weibull latent event times

    v <- runif(n=N)
    event.time <- (- log(v) / (lambdaT * exp(cbind(X,A) %*% beta_T)))^(1 / rho) # + A*covs_int %*% interactions
    k <- runif(n=N)
    cens.time <- (- log(k) / (lambdaC * exp(cbind(X,A) %*% beta_C)))^(1 / rho)
    t0 <- (- log(v) / (lambdaT * exp(cbind(X,0) %*% beta_T)))^(1 / rho)
    t1 <- (- log(v) / (lambdaT * exp(cbind(X,1) %*% beta_T)))^(1 / rho)
    # follow-up times and event indicators

    time <- pmin(event.time, cens.time)
    event <- as.numeric(event.time <= cens.time)
    out <- as.data.frame(cbind(X,RCT,A,time, event.time, cens.time, event, t0, t1))
    colnames(out) <- c(colnames(X),"RCT","A","time", "event.time", "cens.time", "event", "t0", "t1")
    return(out)

}


#' Make data for simulations
#'
#' This function was used to create the simulation datasets.
#' @export
#'


simulWeibDat_nonlinear <- function(N, covs, beta_R, beta_A,
                                   lambdaT, lambdaC, beta_T, beta_C,shapeT, shapeC, rho,
                                   s, interactions = NULL, covs_int = NULL)
{
    set.seed(s)
    X <- matrix(rnorm(n = N*covs,0,1), nrow = N, ncol = covs)
    colnames(X) <- paste0("X",1:covs)

    if(is.null(interactions))
        interactions = rep(0,ncol(X)-1)
    if(is.null(covs_int))
        covs_int = rep(0, N)

    # generate covariates, sampling and treatment #
    expit <- function(x,beta,int = 0) 1/(1 + exp(-x %*% beta)+int)

    gRs <- expit(X, beta_R,1)
    RCT <- rbinom(n = N, size = 1, prob = gRs)
    # sum(RCT)/length(RCT)

    g0s <- ifelse(RCT == 1,1/2 ,expit(X, beta_A))
    A <- rbinom(n=N, size=1, prob=g0s)


    # generate time data #
    # set.seed(s)
    # cens.time <- rweibull(N, shape = shapeC, scale = lambdaC * exp(-cbind(X,A) %*% beta_C))
    # set.seed(s)
    # event.time <- rweibull(N, shape = shapeT, scale = lambdaT * exp(-cbind(X,A) %*% beta_T))
    # set.seed(s)
    # t0 <- rweibull(N, shape = shapeT, scale = lambdaT * exp(-cbind(X,0) %*% beta_T))
    # set.seed(s)
    # t1 <- rweibull(N, shape = shapeT, scale = lambdaT * exp(-cbind(X,1) %*% beta_T))


    # Weibull latent event times

    v <- runif(n=N)
    event.time <- (- log(v) / (lambdaT * exp(cbind(X^3,A) %*% beta_T)))^(1 / rho) # + A*covs_int %*% interactions
    k <- runif(n=N)
    cens.time <- (- log(k) / (lambdaC * exp(cbind(X^3,A) %*% beta_C)))^(1 / rho)
    t0 <- (- log(v) / (lambdaT * exp(cbind(X^3,0) %*% beta_T)))^(1 / rho)
    t1 <- (- log(v) / (lambdaT * exp(cbind(X^3,1) %*% beta_T)))^(1 / rho)
    # follow-up times and event indicators

    time <- pmin(event.time, cens.time)
    event <- as.numeric(event.time <= cens.time)
    out <- as.data.frame(cbind(X,RCT,A,time, event.time, cens.time, event, t0, t1))
    colnames(out) <- c(colnames(X),"RCT","A","time", "event.time", "cens.time", "event", "t0", "t1")
    return(out)

}


#' Function to generate simulation results
#'
#' This function was used to generate the simulation results
#' @export
#' @import
#' dplyr
#' survival

simul_methods_causal <- function(N, covs, beta_R, beta_A,
                                 lambdaT, lambdaC, beta_T, beta_C,shapeT, shapeC, rho = 1,
                                 s, interactions = NULL, covs_int = NULL,
                                 prop.RCT.SL.library = c("SL.rpart","SL.glmnet"),
                                 prop.SL.library = c("SL.mean"),
                                 event.SL.library = c("survSL.coxph", "survSL.weibreg", "survSL.rfsrc"),
                                 cens.SL.library = c("survSL.coxph", "survSL.weibreg", "survSL.rfsrc"),
                                 quants = c(0.75,0.5,0.25),
                                 # parameters for interactions #
                                 misspe = "",
                                 method = "dml",
                                 run_type = "correct",
                                 fit.times = seq(0, 10, by=1),
                                 V = 5, linear = T
){

    # library(dplyr)
    # library(survival)
    set.seed(s)


    if(linear){
        dat <- simulWeibDat(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
                            lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
                            shapeT = shapeT, shapeC = shapeC, rho = rho,
                            s = s+1, interactions = NULL, covs_int = NULL)
        # load("~/Causal/CFsurvival/results/counterfactual_data.Rdata")
        load("~/Causal/CFsurvival/results/S0.Rdata")
        load("~/Causal/CFsurvival/results/S1.Rdata")
    }
    else{
        dat <- simulWeibDat_nonlinear(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
                                      lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
                                      shapeT = shapeT, shapeC = shapeC, rho = rho,
                                      s = s+1, interactions = NULL, covs_int = NULL)
        load("~/Causal/CFsurvival/results/counterfactual_data_nonlinear.Rdata")
    }

    # fit.times = seq(0, 10, by=0.1)
    t0 <- dat$t0
    t1 <- dat$t1
    RCT <- dat$RCT
    A <- dat$A
    X <- dat[,grep("X", colnames(dat))]

    # out_S0 <- dat_counterfactual$S_0
    # out_S1 <- dat_counterfactual$S_1
    out_S0 <- out_S0[match(fit.times, names(out_S0))]
    out_S1 <- out_S1[match(fit.times, names(out_S1))]


    X_mod <- cbind(sqrt(abs(X[,1])/2),exp(X[,2] - X[,1]), log(abs(X[,3]+ X[,4] + 1)), sqrt(abs(X[,4]+X[,5])), (X[,5])^2)

    # return(list(dat = dat, out_S0 = out_S0, out_S1 = out_S1) )

    # set up data needed #
    obs.time <- dat$event.time[dat$RCT == 1]
    obs.event <- dat$event[dat$RCT == 1]
    rx <- dat$A[dat$RCT == 1]
    confounders <- dat[dat$RCT == 1,grep("X", colnames(dat))]
    W_c <- dat[dat$RCT == 0,grep("X", colnames(dat))]
    mod_conf <- X_mod[dat$RCT==1,]
    mod_W_c <- X_mod[dat$RCT ==0,]
    # Run methods associated to each of the cases of interest #

    # Ours #
    if(method == "dml"){

        # correct specificiation #
        if(run_type == "correct"){

            fit <- CFsurvival(time = obs.time, event = obs.event, treat = rx,
                              confounders =  confounders, contrasts = NULL,
                              verbose = FALSE, fit.times = fit.times,
                              fit.treat = c(0,1),
                              nuisance.options = list(
                                  prop.RCT.SL.library = prop.RCT.SL.library,
                                  prop.SL.library = prop.SL.library,
                                  event.SL.library = event.SL.library,
                                  cens.SL.library = cens.SL.library, verbose = F, V = V),
                              W_c = W_c
            )
        }

        # misspecified models #

        if(run_type == "incorrect"){
            if(length(misspe) == 1 & "sampling" %in% misspe)
                fit <- CFsurvival(time = obs.time, event = obs.event, treat = rx,
                                  confounders =  confounders, contrasts = NULL,
                                  verbose = FALSE, fit.times = fit.times,
                                  fit.treat = c(0,1),
                                  nuisance.options = list(
                                      prop.RCT.SL.library = prop.RCT.SL.library,
                                      prop.SL.library = prop.SL.library,
                                      event.SL.library = event.SL.library,
                                      cens.SL.library = cens.SL.library, verbose = F, V = V),
                                  W_c = W_c,misspe = misspe, mod_conf = mod_conf,mod_W_c = mod_W_c
                )

            else if (length(misspe) == 1 & "survival" %in% misspe)
                fit <- CFsurvival(time = obs.time, event = obs.event, treat = rx,
                                  confounders =  confounders, contrasts = NULL,
                                  verbose = FALSE, fit.times = fit.times,
                                  fit.treat = c(0,1),
                                  nuisance.options = list(
                                      prop.RCT.SL.library = prop.RCT.SL.library,
                                      prop.SL.library = prop.SL.library,
                                      event.SL.library = event.SL.library,
                                      cens.SL.library = cens.SL.library, verbose = F, V = V),
                                  W_c = W_c,misspe = misspe, mod_conf = mod_conf,mod_W_c = mod_W_c
                )

            else
                fit <- CFsurvival(time = obs.time, event = obs.event, treat = rx,
                                  confounders =  mod_conf, contrasts = NULL,
                                  verbose = FALSE, fit.times = fit.times,
                                  fit.treat = c(0,1),
                                  nuisance.options = list(
                                      prop.RCT.SL.library = prop.RCT.SL.library,
                                      prop.SL.library = prop.SL.library,
                                      event.SL.library = event.SL.library,
                                      cens.SL.library = cens.SL.library, verbose = F, V = V),
                                  W_c = mod_W_c
                )

        }
    }


    #################################################

    ### Naive model with only coxph outcome model ###

    #################################################

    if(method == "naive"){

        if(!is.null(mod_conf))
            mod_conf <- as.data.frame(mod_conf)
        if(!is.null(mod_W_c))
            mod_W_c <- as.data.frame(mod_W_c)

        nuis <- list()
        nuis$eval.times <- sort(unique(c(0,obs.time[obs.time > 0 & obs.time <= max(fit.times)], max(fit.times))))
        n <- length(obs.time)

        S.hats.0 <- lower_0 <- upper_0 <- matrix(0L, nrow = N, ncol = length(fit.times))
        S.hats.1 <- lower_1 <- upper_1 <- matrix(0L, nrow = N, ncol = length(fit.times))

        # lower_0 <- upper_0 <- lower_1 <- upper_1 <- c(0L)

        if(run_type == "correct")
            data_use <- as.data.frame(cbind(obs.time, obs.event, rx, confounders))
        if(run_type == "incorrect")
            data_use <- as.data.frame(cbind(obs.time, obs.event, rx, mod_conf))
        colnames(data_use) <- c("time", "event", "A",paste0("X",1:covs))

        fit <- coxph(Surv(time, event) ~ A + X1 + X2 + X3 + X4 + X5, data = data_use)

        # no trt #
        if(run_type == "correct")
            new_data <- as.data.frame( cbind(A, X))
        if(run_type == "incorrect")
            new_data <- as.data.frame( cbind(A, X_mod))
        colnames(new_data) <- c("A",paste0("X",1:covs))

        predicted_0 <- survfit(fit, newdata = new_data %>% mutate(A = 0),
                               conf.type = "plain")
        S.hats.0 <- S.hats.0 + t(summary(predicted_0, times = fit.times)$surv)
        lower_0 <- lower_0 + t(summary(predicted_0, times = fit.times)$lower)
        upper_0 <- upper_0 + t(summary(predicted_0, times = fit.times)$upper)

        # trt #
        predicted_1 <- survfit(fit, newdata = new_data %>% mutate(A = 1),
                               conf.type = "plain")
        S.hats.1 <- S.hats.1 +  t(summary(predicted_1, times = fit.times)$surv)
        lower_1 <- lower_1 + t(summary(predicted_1, times = fit.times)$lower)
        upper_1 <- upper_1 + t(summary(predicted_1, times = fit.times)$upper)

        ####
        final_pred_0 <- as_tibble(cbind(fit.times, rep(0,length(fit.times)), apply(S.hats.0,2, mean),out_S0,
                                        apply(lower_0,2, mean),apply(upper_0,2, mean)))
        # final_pred_0 <- as_tibble(cbind(fit.times, rep(0,length(fit.times)), apply(S.hats.0,2, mean),out_S0,
        #                                 as.vector(lower_0),as.vector(upper_0)))
        colnames(final_pred_0) <- c('time', 'trt','surv','true.surv','ptwise.lower','ptwise.upper')

        final_pred_1 <- as_tibble(cbind(fit.times, rep(1,length(fit.times)), apply(S.hats.1,2, mean),out_S1,
                                        apply(lower_1,2, mean),apply(upper_1,2, mean)))
        # final_pred_1 <- as_tibble(cbind(fit.times, rep(1,length(fit.times)), apply(S.hats.1,2, mean),out_S1,
        #                                 as.vector(lower_1),as.vector(upper_1)))
        colnames(final_pred_1) <- c('time', 'trt','surv','true.surv','ptwise.lower','ptwise.upper')

        # merge results #
        surv.df <- rbind(final_pred_0, final_pred_1)
        fit <- list(surv.df = surv.df)


    }

    #######################################################

    # SL outcome #

    #######################################################

    if(method == "SL"){

        ######################################################
        ########## MAKE FOLDS SIMILAR TO THE CF ##############
        ######################################################


        if(!is.null(mod_conf))
            mod_conf <- as.data.frame(mod_conf)
        if(!is.null(mod_W_c))
            mod_W_c <- as.data.frame(mod_W_c)

        nuis <- list()
        nuis$eval.times <- sort(unique(c(0,obs.time[obs.time > 0 & obs.time <= max(fit.times)], max(fit.times))))
        n <- length(obs.time)

        nuis <- do.call("CFsurvival.nuisance.options", list(
            prop.RCT.SL.library = prop.RCT.SL.library,
            prop.SL.library = prop.SL.library,
            event.SL.library = event.SL.library,
            cens.SL.library = cens.SL.library, verbose = F))
        fit.treat <- c(0,1)

        S.hats.0 <- matrix(0L, nrow = N, ncol = length(fit.times))
        S.hats.1 <- matrix(0L, nrow = N, ncol = length(fit.times))


        if(run_type == "correct"){
            W <- as.data.frame(confounders)
            newW <- as.data.frame(X)
            colnames(W) <- colnames(newW) <- c( paste0("X", 1:covs))
        }
        if(run_type == "incorrect"){
            W <- as.data.frame(mod_conf)
            newW <- as.data.frame(X_mod)
            colnames(W) <- colnames(newW) <- c( paste0("X", 1:covs))
        }

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
        S.hats.0 <- S.hats.0 + fit$event.pred.0
        # estimates for events under treat = 1 #
        S.hats.1 <- S.hats.1 + fit$event.pred.1

        S.hats.0 <- S.hats.0
        lower0 <- apply(S.hats.0,2, function(x){quantile(x,0.05)})
        upper0 <- apply(S.hats.0,2, function(x){quantile(x,0.95)})
        final_pred_0 <- as_tibble(cbind(fit.times, rep(0,length(fit.times)), apply(S.hats.0,2, mean),out_S0, lower0, upper0))
        colnames(final_pred_0) <- c('time', 'trt','surv','true.surv','ptwise.lower','ptwise.upper')

        S.hats.1 <- S.hats.1
        lower1 <- apply(S.hats.1,2, function(x){quantile(x,0.05)})
        upper1 <- apply(S.hats.1,2, function(x){quantile(x,0.95)})
        final_pred_1 <- as_tibble(cbind(fit.times, rep(1,length(fit.times)), apply(S.hats.1,2, mean),out_S1, lower1, upper1))
        colnames(final_pred_1) <- c('time', 'trt','surv','true.surv','ptwise.lower','ptwise.upper')


        # merge results #
        surv.df <- rbind(final_pred_0, final_pred_1)
        fit <- list(surv.df = surv.df)
    }




    ##########################################################

    # Coxph weigthed #

    ##########################################################

    if(method == "weighted"){

        if(!is.null(mod_conf))
            mod_conf <- as.data.frame(mod_conf)
        if(!is.null(mod_W_c))
            mod_W_c <- as.data.frame(mod_W_c)

        nuis <- list()
        nuis$eval.times <- sort(unique(c(0,obs.time[obs.time > 0 & obs.time <= max(fit.times)], max(fit.times))))
        n <- length(obs.time)

        S.hats.0 <- lower_0 <- upper_0 <- matrix(0L, nrow = N, ncol = length(fit.times))
        S.hats.1 <- lower_1 <- upper_1 <- matrix(0L, nrow = N, ncol = length(fit.times))

        if(run_type == "correct"){
            temp_RCT <- cbind(rep(1, sum(RCT == 1)), confounders)
            temp_obs <- cbind(rep(0, nrow(W_c)), W_c)
            new_data <- as.data.frame( cbind(A, X))
        }
        if(run_type == "incorrect"){

            if(length(misspe) == 1 & "sampling" %in% misspe){
                temp_RCT <- cbind(rep(1, sum(RCT == 1)), mod_conf)
                temp_obs <- cbind(rep(0, nrow(W_c)), mod_W_c)
                new_data <- as.data.frame( cbind(A, X))
            }

            else if (length(misspe) == 1 & "survival" %in% misspe){
                temp_RCT <- cbind(rep(1, sum(RCT == 1)), confounders)
                temp_obs <- cbind(rep(0, nrow(W_c)), W_c)
                new_data <- as.data.frame( cbind(A, X_mod))
            }

            else{
                temp_RCT <- cbind(rep(1, sum(RCT == 1)), mod_conf)
                temp_obs <- cbind(rep(0, nrow(W_c)), mod_W_c)
                new_data <- as.data.frame( cbind(A, X_mod))
            }
        }


        colnames(temp_obs) <- colnames(temp_RCT) <- c("A",paste0("X",1:covs))
        data_use <- as.data.frame(rbind(temp_RCT, temp_obs) )
        # as.data.frame(cbind(RCT, X))
        colnames(data_use) <- c("RCT",paste0("X",1:covs))
        logit_fit <- glm(RCT ~ ., data = data_use)

        weigths_use <- as.numeric(predict(logit_fit))
        weigths_use <- 1/ifelse(weigths_use <= 0, 0.01,weigths_use)[data_use$RCT == 1]

        data_use <- as.data.frame(cbind(obs.time, obs.event, A[RCT == 1],X[RCT == 1,]))
        colnames(data_use) <- c("time", "event", "A",paste0("X",1:covs))
        fit <- coxph(Surv(time, event) ~ A + X1 + X2 + X3 + X4 + X5, data = data_use,weights = weigths_use)

        # no trt #
        colnames(new_data) <- c("A",paste0("X",1:covs))

        predicted_0 <- survfit(fit, newdata = new_data %>% mutate(A = 0),
                               conf.type = "plain")
        S.hats.0 <- S.hats.0 + t(summary(predicted_0, times = fit.times)$surv)
        lower_0 <- lower_0 + t(summary(predicted_0, times = fit.times)$lower)
        upper_0 <- upper_0 + t(summary(predicted_0, times = fit.times)$upper)


        # trt #
        predicted_1 <- survfit(fit, newdata = new_data %>% mutate(A = 1),
                               conf.type = "plain")
        S.hats.1 <- S.hats.1 +  t(summary(predicted_1, times = fit.times)$surv)
        lower_1 <- lower_1 + t(summary(predicted_1, times = fit.times)$lower)
        upper_1 <- upper_1 + t(summary(predicted_1, times = fit.times)$upper)



        final_pred_0 <- as_tibble(cbind(fit.times, rep(0,length(fit.times)), apply(S.hats.0,2, mean),out_S0,
                                        apply(lower_0,2, mean),apply(upper_0,2, mean)))
        colnames(final_pred_0) <- c('time', 'trt','surv','true.surv','ptwise.lower','ptwise.upper')


        final_pred_1 <- as_tibble(cbind(fit.times, rep(1,length(fit.times)), apply(S.hats.1,2, mean),out_S1,
                                        apply(lower_1,2, mean),apply(upper_1,2, mean)))
        colnames(final_pred_1) <- c('time', 'trt','surv','true.surv','ptwise.lower','ptwise.upper')

        # merge results #
        surv.df <- rbind(final_pred_0, final_pred_1)
        fit <- list(surv.df = surv.df)

    }




    ##########################################################

    # SL weigthed Coxph #

    ##########################################################

    if(method == "SLweighted"){


        if(!is.null(mod_conf))
            mod_conf <- as.data.frame(mod_conf)
        if(!is.null(mod_W_c))
            mod_W_c <- as.data.frame(mod_W_c)

        nuis <- list()
        nuis$eval.times <- sort(unique(c(0,obs.time[obs.time > 0 & obs.time <= max(fit.times)], max(fit.times))))
        n <- length(obs.time)


        nuis <- do.call("CFsurvival.nuisance.options", list(
            prop.RCT.SL.library = prop.RCT.SL.library,
            prop.SL.library = prop.SL.library,
            event.SL.library = event.SL.library,
            cens.SL.library = cens.SL.library, verbose = F))

        S.hats.0 <- lower_0 <- upper_0 <- matrix(0L, nrow = N, ncol = length(fit.times))
        S.hats.1 <- lower_1 <- upper_1 <- matrix(0L, nrow = N, ncol = length(fit.times))


        if(run_type == "correct"){
            temp_RCT <- cbind(rep(1, sum(RCT == 1)), confounders)
            temp_obs <- cbind(rep(0, nrow(W_c)), W_c)
            new_data <- as.data.frame( cbind(A, X))
        }
        if(run_type == "incorrect"){

            if(length(misspe) == 1 & "sampling" %in% misspe){
                temp_RCT <- cbind(rep(1, sum(RCT == 1)), mod_conf)
                temp_obs <- cbind(rep(0, nrow(W_c)), mod_W_c)
                new_data <- as.data.frame( cbind(A, X))
            }

            else if (length(misspe) == 1 & "survival" %in% misspe){
                temp_RCT <- cbind(rep(1, sum(RCT == 1)), confounders)
                temp_obs <- cbind(rep(0, nrow(W_c)), W_c)
                new_data <- as.data.frame( cbind(A, X_mod))
            }

            else{
                temp_RCT <- cbind(rep(1, sum(RCT == 1)), mod_conf)
                temp_obs <- cbind(rep(0, nrow(W_c)), mod_W_c)
                new_data <- as.data.frame( cbind(A, X_mod))
            }
        }


        colnames(temp_obs) <- colnames(temp_RCT) <- c("A",paste0("X",1:covs))
        data_use <- as.data.frame(rbind(temp_RCT, temp_obs) )
        # as.data.frame(cbind(RCT, X))
        colnames(data_use) <- c("RCT",paste0("X",1:covs))

        weigths_use <- as.numeric(unlist(.estimate.RCT_cohort.propensity(A=data_use$RCT,
                                                                         W=data_use[,paste0("X",1:covs)],
                                                                         newW=data_use[,paste0("X",1:covs)],
                                                                         SL.library=prop.RCT.SL.library,
                                                                         fit.treat=c(0,1),
                                                                         prop.trunc=nuis$prop.trunc,
                                                                         save.fit = nuis$save.nuis.fits,
                                                                         verbose = nuis$verbose)))

        weigths_use <- 1/ifelse(weigths_use <= 0, 0.01,weigths_use)[data_use$RCT == 1]

        data_use <- as.data.frame(cbind(obs.time, obs.event, A[RCT == 1],X[RCT == 1,]))
        colnames(data_use) <- c("time", "event", "A",paste0("X",1:covs))
        fit <- coxph(Surv(time, event) ~ A + X1 + X2 + X3 + X4 + X5, data = data_use,weights = weigths_use)

        # no trt #
        colnames(new_data) <- c("A",paste0("X",1:covs))

        predicted_0 <- survfit(fit, newdata = new_data %>% mutate(A = 0),
                               conf.type = "plain")
        S.hats.0 <- S.hats.0 + t(summary(predicted_0, times = fit.times)$surv)
        lower_0 <- lower_0 + t(summary(predicted_0, times = fit.times)$lower)
        upper_0 <- upper_0 + t(summary(predicted_0, times = fit.times)$upper)


        # trt #
        predicted_1 <- survfit(fit, newdata = new_data %>% mutate(A = 1),
                               conf.type = "plain")
        S.hats.1 <- S.hats.1 +  t(summary(predicted_1, times = fit.times)$surv)
        lower_1 <- lower_1 + t(summary(predicted_1, times = fit.times)$lower)
        upper_1 <- upper_1 + t(summary(predicted_1, times = fit.times)$upper)



        final_pred_0 <- as_tibble(cbind(fit.times, rep(0,length(fit.times)), apply(S.hats.0,2, mean),out_S0,
                                        apply(lower_0,2, mean),apply(upper_0,2, mean)))
        colnames(final_pred_0) <- c('time', 'trt','surv','true.surv','ptwise.lower','ptwise.upper')


        final_pred_1 <- as_tibble(cbind(fit.times, rep(1,length(fit.times)), apply(S.hats.1,2, mean),out_S1,
                                        apply(lower_1,2, mean),apply(upper_1,2, mean)))
        colnames(final_pred_1) <- c('time', 'trt','surv','true.surv','ptwise.lower','ptwise.upper')

        # merge results #
        surv.df <- rbind(final_pred_0, final_pred_1)
        fit <- list(surv.df = surv.df)


    }

    ####################
    if(method == "dml")
        fit$surv.df$true.surv <- c(out_S0[1:(nrow(fit$surv.df)/2)], out_S1[1:(nrow(fit$surv.df)/2)])
    ###### metrics #####

    ### Overall ###

    ### bias survival ###

    # shift focus to S1 - S0 (difference of the survival curves) #
    Delta_bias <- function(t){
        temp <- fit$surv.df %>% filter(time == t)
        diff_surv_true <- temp$true.surv[temp$trt == 1] - temp$true.surv[temp$trt == 0]
        diff_surv_est <- temp$surv[temp$trt == 1] - temp$surv[temp$trt == 0]
        return(diff_surv_est - diff_surv_true)
    }
    Delta_bias_out <- sapply(X = fit.times, FUN = Delta_bias)


    # Delta_0_bias <- function(t){
    #     temp <- fit$surv.df %>% filter(trt == 0, time == t)
    #     sum(temp$surv - temp$true.surv)
    # }
    # Delta_0_bias_out <- sapply(X = fit.times, FUN = Delta_0_bias)
    # Delta_1_bias <- function(t){
    #     temp <- fit$surv.df %>% filter(trt == 1, time == t)
    #     sum(temp$surv - temp$true.surv)
    # }
    # Delta_1_bias_out <- sapply(X = fit.times, FUN = Delta_1_bias)

    ### MSE survival ###
    # Delta_0_MSE <- function(t){
    #     temp <- fit$surv.df %>% filter(trt == 0, time == t)
    #     sum((temp$surv - temp$true.surv)^2 )
    # }
    # Delta_0_MSE_out <- sapply(X = fit.times, FUN = Delta_0_MSE)
    #
    # Delta_1_MSE <- function(t){
    #     temp <- fit$surv.df %>% filter(trt == 1, time == t)
    #     sum(temp$surv - temp$true.surv)
    # }
    # Delta_1_MSE_out <- sapply(X = fit.times, FUN = Delta_1_MSE)


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


    if(method != "dml")
        return(list(
            # "dat" = dat,
            "fit.surv" = fit$surv.df,
            "bias" = Delta_bias_out,
            # "bias0" = Delta_0_bias_out,
            # "bias1" = Delta_1_bias_out,
            # "RMSE0" = Delta_0_MSE_out,
            # "RMSE1" = Delta_1_MSE_out,
            "quant_bias0" = quant_bias_0,
            "quant_bias1" = quant_bias_1,
            "quant_MSE0" = quant_MSE_0,
            "quant_MSE1" = quant_MSE_1
        ))
    if(method == "dml")
        return(list(
            # "dat" = dat,
            "fit.surv" = fit$surv.df,
            "bias" = Delta_bias_out,
            "EIF0" = fit$IF.vals.0,
            "EIF1" = fit$IF.vals.1,
            # "bias0" = Delta_0_bias_out,
            # "bias1" = Delta_1_bias_out,
            # "RMSE0" = Delta_0_MSE_out,
            # "RMSE1" = Delta_1_MSE_out,
            "quant_bias0" = quant_bias_0,
            "quant_bias1" = quant_bias_1,
            "quant_MSE0" = quant_MSE_0,
            "quant_MSE1" = quant_MSE_1
        ))

}


############################################


# rm(list=ls())
# N = 500
# s = 10
#
# # parameters #
# covs <- 5
# rho <- 1
# rateC = 0.05
# covs <- 5
# V <- 5
#
# beta_R = rep(0.5, covs)
# beta_A = rep(0.5, covs)
# beta_T = c(rep(0.5, covs),-1)
# beta_C = c(rep(-1, covs),-1/5)
#
# shapeC = 0.5
# shapeT = 1
# time_c = NULL
# event_c = NULL
# treat_c = NULL
# conf.level=.95
# conf.band = 0.95
# fit.times = seq(0,10, 1)
# prop.RCT.SL.library = c("SL.glm","SL.rpart","SL.glmnet")
# prop.SL.library = c("SL.mean")
# event.SL.library = c("survSL.coxph", "survSL.weibreg", "survSL.rfsrc")
# cens.SL.library = c("survSL.coxph", "survSL.weibreg", "survSL.rfsrc")
# lambdaC = 0.01
# lambdaT = 0.2
# method = "SLweighted"
# run_type = "correct"
# linear = T


#
# test_correct <- simul_methods_causal(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
#                                      lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
#                                      shapeT = shapeT, shapeC = shapeC,
#                                      s = s, interactions = NULL, covs_int = NULL,V = V,
#                                      prop.RCT.SL.library = prop.RCT.SL.library, prop.SL.library = prop.SL.library,
#                                      event.SL.library = event.SL.library, cens.SL.library = cens.SL.library,
#                                      method = method, run_type = run_type, linear = linear)
#
# test_correct$bias
#
# run_type = "incorrect"
#
#
# test_incorrect <- simul_methods_causal(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
#                                        lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
#                                        shapeT = shapeT, shapeC = shapeC,
#                                        s = s, interactions = NULL, covs_int = NULL,V = V,
#                                        prop.RCT.SL.library = prop.RCT.SL.library, prop.SL.library = prop.SL.library,
#                                        event.SL.library = event.SL.library, cens.SL.library = cens.SL.library,
#                                        method = method, run_type = run_type, linear = linear,misspe = "survival")
#
# test_incorrect$bias
#
#
# test_correct$fit.surv
# test_incorrect$fit.surv
#
# test_correct$fit.surv %>%
#     select(time,trt, surv, true.surv, ptwise.lower, ptwise.upper) %>%
#     # group_by(trt, time) %>%
#     # rowwise() %>%
#     mutate(coverage = ptwise.lower <= true.surv & ptwise.upper >= true.surv)
#
# test_incorrect$fit.surv %>%
#     select(time,trt, surv, true.surv, ptwise.lower, ptwise.upper) %>%
#     # group_by(trt, time) %>%
#     # rowwise() %>%
#     mutate(coverage = ptwise.lower <= true.surv & ptwise.upper >= true.surv)





# # weigthed #
#
# # correct #
# test_correct <- simul_methods_causal(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
#                                      lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
#                                      shapeT = shapeT, shapeC = shapeC,
#                                      s = s, interactions = NULL, covs_int = NULL,V = V,
#                                      prop.RCT.SL.library = prop.RCT.SL.library, prop.SL.library = prop.SL.library,
#                                      event.SL.library = event.SL.library, cens.SL.library = cens.SL.library,
#                                      method = "SLweighted")
#
# # test_correct
#
#
# # sampling #
# test_sampling <- simul_methods_causal(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
#                                       lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
#                                       shapeT = shapeT, shapeC = shapeC,
#                                       s = s, interactions = NULL, covs_int = NULL,
#                                       run_type = "incorrect", misspe = "sampling", V = V,
#                                       prop.RCT.SL.library = prop.RCT.SL.library, prop.SL.library = prop.SL.library,
#                                       event.SL.library = event.SL.library, cens.SL.library = cens.SL.library,
#                                       method = "SLweighted")
#
# # test_sampling
#
#
# # survival #
# test_survival <- simul_methods_causal(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
#                                       lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
#                                       shapeT = shapeT, shapeC = shapeC,
#                                       s = s, interactions = NULL, covs_int = NULL,
#                                       run_type = "incorrect", misspe = "survival", V = V,
#                                       prop.RCT.SL.library = prop.RCT.SL.library, prop.SL.library = prop.SL.library,
#                                       event.SL.library = event.SL.library, cens.SL.library = cens.SL.library,
#                                       method = "SLweighted")
#
# # test_survival
#
#
# # incorrect #
# test_incorrect <- simul_methods_causal(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
#                                        lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
#                                        shapeT = shapeT, shapeC = shapeC,
#                                        s = s, interactions = NULL, covs_int = NULL,
#                                        run_type = "incorrect", misspe = c("sampling","survival"), V = V,
#                                        prop.RCT.SL.library = prop.RCT.SL.library, prop.SL.library = prop.SL.library,
#                                        event.SL.library = event.SL.library, cens.SL.library = cens.SL.library,
#                                        method = "SLweighted")
#
# # test_incorrect
#
# test_correct$bias0
# test_sampling$bias0
# test_survival$bias0
# test_incorrect$bias0
#
#
# test_correct$bias1
# test_sampling$bias1
# test_survival$bias1
# test_incorrect$bias1




#
# # weigthed #
#
# # correct #
# test_correct <- simul_methods_causal(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
#                                      lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
#                                      shapeT = shapeT, shapeC = shapeC,
#                                      s = s, interactions = NULL, covs_int = NULL,V = V,
#                                      prop.RCT.SL.library = prop.RCT.SL.library, prop.SL.library = prop.SL.library,
#                                      event.SL.library = event.SL.library, cens.SL.library = cens.SL.library,
#                                      method = "weighted")
#
# # test_correct
#
#
# # sampling #
# test_sampling <- simul_methods_causal(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
#                                       lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
#                                       shapeT = shapeT, shapeC = shapeC,
#                                       s = s, interactions = NULL, covs_int = NULL,
#                                       run_type = "incorrect", misspe = "sampling", V = V,
#                                       prop.RCT.SL.library = prop.RCT.SL.library, prop.SL.library = prop.SL.library,
#                                       event.SL.library = event.SL.library, cens.SL.library = cens.SL.library,
#                                       method = "weighted")
#
# # test_sampling
#
#
# # survival #
# test_survival <- simul_methods_causal(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
#                                       lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
#                                       shapeT = shapeT, shapeC = shapeC,
#                                       s = s, interactions = NULL, covs_int = NULL,
#                                       run_type = "incorrect", misspe = "survival", V = V,
#                                       prop.RCT.SL.library = prop.RCT.SL.library, prop.SL.library = prop.SL.library,
#                                       event.SL.library = event.SL.library, cens.SL.library = cens.SL.library,
#                                       method = "weighted")
#
# # test_survival
#
#
# # incorrect #
# test_incorrect <- simul_methods_causal(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
#                                        lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
#                                        shapeT = shapeT, shapeC = shapeC,
#                                        s = s, interactions = NULL, covs_int = NULL,
#                                        run_type = "incorrect", misspe = c("sampling","survival"), V = V,
#                                        prop.RCT.SL.library = prop.RCT.SL.library, prop.SL.library = prop.SL.library,
#                                        event.SL.library = event.SL.library, cens.SL.library = cens.SL.library,
#                                        method = "weighted")
#
# # test_incorrect
#
# test_correct$bias0
# test_sampling$bias0
# test_survival$bias0
# test_incorrect$bias0
#
#
# test_correct$bias1
# test_sampling$bias1
# test_survival$bias1
# test_incorrect$bias1

#
# test_correct <- simul_methods_causal(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
#                                      lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
#                                      shapeT = shapeT, shapeC = shapeC,
#                                      s = s, interactions = NULL, covs_int = NULL,V = V,
#                                      prop.RCT.SL.library = prop.RCT.SL.library, prop.SL.library = prop.SL.library,
#                                      event.SL.library = event.SL.library, cens.SL.library = cens.SL.library,
#                                      method = "SL")
#
# test_correct
#
# test_incorrect <- simul_methods_causal(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
#                                      lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
#                                      shapeT = shapeT, shapeC = shapeC,
#                                      s = s, interactions = NULL, covs_int = NULL,V = V,
#                                      prop.RCT.SL.library = prop.RCT.SL.library, prop.SL.library = prop.SL.library,
#                                      event.SL.library = event.SL.library, cens.SL.library = cens.SL.library,
#                                      method = "SL",run_type = "incorrect")
#
# test_incorrect

#
# test_correct <- simul_methods_causal(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
#                                      lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
#                                      shapeT = shapeT, shapeC = shapeC,
#                                      s = s, interactions = NULL, covs_int = NULL,V = V,
#                                      prop.RCT.SL.library = prop.RCT.SL.library, prop.SL.library = prop.SL.library,
#                                      event.SL.library = event.SL.library, cens.SL.library = cens.SL.library,
#                                      method = "naive")
#
# test_correct
#
# test_incorrect <- simul_methods_causal(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
#                                      lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
#                                      shapeT = shapeT, shapeC = shapeC,
#                                      s = s, interactions = NULL, covs_int = NULL,V = V,
#                                      prop.RCT.SL.library = prop.RCT.SL.library, prop.SL.library = prop.SL.library,
#                                      event.SL.library = event.SL.library, cens.SL.library = cens.SL.library,
#                                      method = "naive",run_type = "incorrect")
#
# test_incorrect


# dat <- simulWeibDat(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
#                     lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
#                     shapeT = shapeT, shapeC = shapeC,
#                     s = s+1, interactions = NULL, covs_int = NULL)
# summary(dat$time)
# summary(dat$event)
# # correct #
# test_correct <- simul_methods_causal(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
#                                      lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
#                                      shapeT = shapeT, shapeC = shapeC,
#                                      s = s, interactions = NULL, covs_int = NULL,V = V,
#                                      prop.RCT.SL.library = prop.RCT.SL.library, prop.SL.library = prop.SL.library,
#                                      event.SL.library = event.SL.library, cens.SL.library = cens.SL.library)
#
# # test_correct
#
#
# # sampling #
# test_sampling <- simul_methods_causal(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
#                                       lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
#                                       shapeT = shapeT, shapeC = shapeC,
#                                       s = s, interactions = NULL, covs_int = NULL,
#                                       run_type = "incorrect", misspe = "sampling", V = V,
#                                       prop.RCT.SL.library = prop.RCT.SL.library, prop.SL.library = prop.SL.library,
#                                       event.SL.library = event.SL.library, cens.SL.library = cens.SL.library)
#
# # test_sampling
#
#
# # survival #
# test_survival <- simul_methods_causal(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
#                                       lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
#                                       shapeT = shapeT, shapeC = shapeC,
#                                       s = s, interactions = NULL, covs_int = NULL,
#                                       run_type = "incorrect", misspe = "survival", V = V,
#                                       prop.RCT.SL.library = prop.RCT.SL.library, prop.SL.library = prop.SL.library,
#                                       event.SL.library = event.SL.library, cens.SL.library = cens.SL.library)
#
# # test_survival
#
#
# # incorrect #
# test_incorrect <- simul_methods_causal(N = N, covs = covs, beta_R = beta_R, beta_A = beta_A,
#                                        lambdaT = lambdaT, lambdaC = lambdaC, beta_T = beta_T, beta_C = beta_C,
#                                        shapeT = shapeT, shapeC = shapeC,
#                                        s = s, interactions = NULL, covs_int = NULL,
#                                        run_type = "incorrect", misspe = c("sampling","survival"), V = V,
#                                        prop.RCT.SL.library = prop.RCT.SL.library, prop.SL.library = prop.SL.library,
#                                        event.SL.library = event.SL.library, cens.SL.library = cens.SL.library)
#
# # test_incorrect
#
#
# # results #
# # test_correct$bias0
# # test_correct$bias1
# # test_correct$RMSE0
# # test_correct$RMSE1
# # test_correct$quant_bias0
# # test_correct$quant_bias1
# #
# #
# # test_sampling$bias0
# # test_sampling$bias1
# # test_sampling$RMSE0
# # test_sampling$RMSE1
# # test_sampling$quant_bias0
# # test_sampling$quant_bias1
# #
# #
# # test_survival$bias0
# # test_survival$bias1
# # test_survival$RMSE0
# # test_survival$RMSE1
# # test_survival$quant_bias0
# # test_survival$quant_bias1
# #
# #
# # test_incorrect$bias0
# # test_incorrect$bias1
# # test_incorrect$RMSE0
# # test_incorrect$RMSE1
# # test_incorrect$quant_bias0
# # test_incorrect$quant_bias1
#
#
# test_correct$bias0 - test_sampling$bias0
# test_correct$bias0 - test_survival$bias0
# test_correct$bias0 - test_incorrect$bias0
#
#
# test_correct$bias1 - test_sampling$bias1
# test_correct$bias1 - test_survival$bias1
# test_correct$bias1 - test_incorrect$bias1
#
#
# test_correct$RMSE0 - test_sampling$RMSE0
# test_correct$RMSE0 - test_survival$RMSE0
# test_correct$RMSE0 - test_incorrect$RMSE0
#
#
# test_correct$RMSE1 - test_sampling$RMSE1
# test_correct$RMSE1 - test_survival$RMSE1
# test_correct$RMSE1 - test_incorrect$RMSE1
#
# # Check data stuff ###
# # does the counterfactual values fit expectations #
# # test_correct$dat %>%
# #     filter(event == 1, A == 0) %>%
# #     select(time, t0,t1)
# #
# # test_correct$dat %>%
# #     filter(event == 1, A == 1) %>%
# #     select(time, t0,t1)
# # good #
#
# # do we get the correct-ish betas #
# library(survival)
# dat <- test_correct$dat
#
# # event #
# fit.event <- coxph(Surv(time, event) ~ X1 + X2 + X3 + X4 + X5 + A, data = dat)
# summary(fit.event)
# # yep #
#
# # censoring #
# fit.cens <- coxph(Surv(time, 1 - event) ~ X1 + X2 + X3 + X4 + X5 + A, data = dat)
# summary(fit.cens)
# # Close-sh #
