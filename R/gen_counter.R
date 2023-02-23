#' function to generate counterfactual data
#'
#' Generate counterfactual data for simulations + find the
#' empirical efficiency bound of our method
#' @export

generate_counterfactual <- function(N = 10^6, beta_T, lambdaT, rho,
                                    beta_C, lambdaC,
                                    beta_R, beta_A,
                                    s, covs = 5, fit.times = seq(0,10,1)){

    # generate covariates #
    set.seed(s)
    X <- matrix(rnorm(n = N*covs,0,1), nrow = N, ncol = covs)
    colnames(X) <- paste0("X",1:covs)


    # generate covariates, sampling and treatment #
    expit <- function(x,beta,int = 0) 1/(1 + exp(-x %*% beta)+int)

    gRs <- as.vector(expit(X, beta_R,1))
    RCT <- rbinom(n = N, size = 1, prob = gRs)
    # sum(RCT)/length(RCT)

    g0s <- ifelse(RCT == 1,1/2 ,expit(X, beta_A))
    A <- rbinom(n=N, size=1, prob=g0s)

    # help function #
    H <- function(lambda, t, rho){
        return(lambda * t^rho)
    }

    # make counterfactual event time #
    S_0_temp <- sapply(fit.times, function(t)
        exp(- H(lambdaT, t, rho) * exp(cbind(X,0) %*% beta_T)))
    S_1_temp <- sapply(fit.times, function(t)
        exp(- H(lambdaT, t, rho) * exp(cbind(X,1) %*% beta_T)))
    colnames(S_0_temp) <- fit.times
    colnames(S_1_temp) <- fit.times
    # get mean survival at each time point #
    S_0 <- apply(S_0_temp, 2, mean)
    S_1 <- apply(S_1_temp, 2, mean)
    names(S_0) <- fit.times
    names(S_1) <- fit.times

    # make counterfactual event time #
    G_0_temp <- sapply(fit.times, function(t)
        exp(- H(lambdaC, t, rho) * exp(cbind(X,0) %*% beta_C)))
    G_1_temp <- sapply(fit.times, function(t)
        exp(- H(lambdaC, t, rho) * exp(cbind(X,1) %*% beta_C)))

    # get mean survival at each time point #
    G_0 <- apply(G_0_temp, 2, mean)
    G_1 <- apply(G_1_temp, 2, mean)
    names(G_0) <- fit.times
    names(G_1) <- fit.times

    v <- runif(n=N)
    event.time <- (- log(v) / (lambdaT * exp(cbind(X,A) %*% beta_T)))^(1 / rho)
    k <- runif(n=N)
    cens.time <- (- log(k) / (lambdaC * exp(cbind(X,A) %*% beta_C)))^(1 / rho)

    time <- as.vector(pmin(event.time, cens.time))
    event <- as.numeric(event.time <= cens.time)


    # get the efficiency bound #
    # first compute EIF for each combination of t and X generated #
    EIFs_0 <- sapply(fit.times, function(t){

        # print(t)
        # get approx value of S(y|x,a) and G(y|x,a) for each observation #
        S_times <- S_0[sapply(time, function(x){which.min(abs(x - fit.times))})]
        G_times <- G_0[sapply(time, function(x){which.min(abs(x - fit.times))})]

        # get approx for integral #
        S_int <- sapply(sapply(time, function(x){
            S_0[1:which.min(abs(x - t))]
        }), sum)
        G_int <- sapply(sapply(time, function(x){
            G_0[1:which.min(abs(x - t))]
        }), sum)


        # NEED TO FIGURE OUT THIS INTEGRAL !!!! #
        int.vals <- 1 - 1/(S_int * G_int)# diff(1/(S_int))*1/G_int

        EIF <- S_0_temp[,as.character(t)]*(1 - ((A == 1)*(RCT == 1))/(gRs * g0s) *
                                               (((time <= t)*(event == 1))/( S_times * G_times ) -
                                                    int.vals
                                               )
        )

        return(EIF)
    })


    EIFs_1 <- sapply(fit.times, function(t){

        # get approx value of S(y|x,a) and G(y|x,a) for each observation #
        S_times <- S_1[sapply(time, function(x){which.min(abs(x - fit.times))})]
        G_times <- G_1[sapply(time, function(x){which.min(abs(x - fit.times))})]

        # get approx for integral #
        S_int <- sapply(sapply(time, function(x){
            S_1[1:which.min(abs(x - t))]
        }), sum)
        G_int <- sapply(sapply(time, function(x){
            G_1[1:which.min(abs(x - t))]
        }), sum)

        int.vals <- 1 - 1/(S_int * G_int)# diff(1/(S_int))*1/G_int

        EIF <- S_1_temp[,as.character(t)]*(1 - ((A == 1)*(RCT == 1))/(gRs * g0s) *
                                               (((time <= t)*(event == 1))/( S_times * G_times ) -
                                                    int.vals
                                               )
        )
        return(EIF)
    })

    ## get the efficiency bound by taking the variance at each time point ##
    EB_0 <- apply(EIFs_0, 2, var)
    EB_1 <- apply(EIFs_1, 2, var)

    return(list(S_0 = S_0, S_1 = S_1,
                G_0 = G_0, G_1 = G_1,
                X = X, event.time = event.time, cens.time = cens.time,
                EIFs_0 = EIFs_0, EIFs_1 = EIFs_1, EB_0 = EB_0, EB_1 = EB_1))
}

# N = 500
# covs = 5
# lambdaT = 0.2
# beta_T = c(rep(0.5, covs),-1)
# s = 1
# rho = 1
# fit.times = seq(0,10,0.1)
# beta_C = c(rep(-1, covs),-1/5)
# lambdaC = 0.01
# beta_R = rep(0.5, covs)
# beta_A = rep(0.5, covs)
#
# dat_counterfactual <- generate_counterfactual(N = N, covs = covs, lambdaT = lambdaT, beta_T = beta_T,
#                                 lambdaC = lambdaC, beta_C = beta_C, beta_R = beta_R, beta_A = beta_A,
#                                 s = s, rho = rho, fit.times = fit.times)
#
# # save(dat_counterfactual, file = "./results/counterfactual_data.Rdata")
# library(ggplot2)
# library(reshape2)
# library(dplyr)
# toy_ex <- generate_counterfactual(N = N, covs = covs, lambdaT = lambdaT, beta_T = beta_T,
#                                               lambdaC = lambdaC, beta_C = beta_C, beta_R = beta_R, beta_A = beta_A,
#                                               s = s, rho = rho, fit.times = fit.times)
# toy_ex <- as.data.frame(cbind(c(test$S_0, test$S_1), rep(fit.times, 2),
#       c(rep("S_0", length(fit.times)),
#         rep("S_1", length(fit.times)))))
# colnames(toy_ex) <- c("Surv", "t", "type")
# toy_ex %>%
#     mutate(Surv = as.numeric(as.character(Surv)),
#            t = as.numeric(as.character(t))) %>%
#     ggplot(aes(x = t, y = Surv, group = type, color = type)) + geom_line()




#' function to generate counterfactual data
#'
#' Generate counterfactual data for simulations + find the
#' empirical efficiency bound of our method
#' @export

generate_counterfactual_nonlinear <- function(N = 10^6, beta_T, lambdaT, rho,
                                              beta_C, lambdaC,
                                              beta_R, beta_A,
                                              s, covs = 5, fit.times = seq(0,10,1)){

    # generate covariates #
    set.seed(s)
    X <- matrix(rnorm(n = N*covs,0,1), nrow = N, ncol = covs)
    colnames(X) <- paste0("X",1:covs)


    # generate covariates, sampling and treatment #
    expit <- function(x,beta,int = 0) 1/(1 + exp(-x %*% beta)+int)

    gRs <- as.vector(expit(X, beta_R,1))
    RCT <- rbinom(n = N, size = 1, prob = gRs)
    # sum(RCT)/length(RCT)

    g0s <- ifelse(RCT == 1,1/2 ,expit(X, beta_A))
    A <- rbinom(n=N, size=1, prob=g0s)

    # help function #
    H <- function(lambda, t, rho){
        return(lambda * t^rho)
    }

    # make counterfactual event time #
    # S_0_temp <- sapply(fit.times, function(t)
    #     exp(- H(lambdaT, t, rho) * exp(cbind(exp(X),0) %*% beta_T)))
    # S_1_temp <- sapply(fit.times, function(t)
    #     exp(- H(lambdaT, t, rho) * exp(cbind(exp(X),1) %*% beta_T)))
    S_0_temp <- sapply(fit.times, function(t)
        exp(- H(lambdaT, t, rho) * exp(cbind(X^3,0) %*% beta_T)))
    S_1_temp <- sapply(fit.times, function(t)
        exp(- H(lambdaT, t, rho) * exp(cbind(X^3,1) %*% beta_T)))
    colnames(S_0_temp) <- fit.times
    colnames(S_1_temp) <- fit.times
    # get mean survival at each time point #
    S_0 <- apply(S_0_temp, 2, mean)
    S_1 <- apply(S_1_temp, 2, mean)
    names(S_0) <- fit.times
    names(S_1) <- fit.times

    # make counterfactual event time #
    G_0_temp <- sapply(fit.times, function(t)
        exp(- H(lambdaC, t, rho) * exp(cbind(X^3,0) %*% beta_C)))
    G_1_temp <- sapply(fit.times, function(t)
        exp(- H(lambdaC, t, rho) * exp(cbind(X^3,1) %*% beta_C)))

    # get mean survival at each time point #
    G_0 <- apply(G_0_temp, 2, mean)
    G_1 <- apply(G_1_temp, 2, mean)
    names(G_0) <- fit.times
    names(G_1) <- fit.times

    v <- runif(n=N)
    event.time <- (- log(v) / (lambdaT * exp(cbind(X^3,A) %*% beta_T)))^(1 / rho)
    k <- runif(n=N)
    cens.time <- (- log(k) / (lambdaC * exp(cbind(X^3,A) %*% beta_C)))^(1 / rho)

    time <- as.vector(pmin(event.time, cens.time))
    event <- as.numeric(event.time <= cens.time)


    # get the efficiency bound #
    # first compute EIF for each combination of t and X generated #
    EIFs_0 <- sapply(fit.times, function(t){

        # print(t)
        # get approx value of S(y|x,a) and G(y|x,a) for each observation #
        S_times <- S_0[sapply(time, function(x){which.min(abs(x - fit.times))})]
        G_times <- G_0[sapply(time, function(x){which.min(abs(x - fit.times))})]

        # get approx for integral #
        S_int <- sapply(sapply(time, function(x){
            S_0[1:which.min(abs(x - t))]
        }), sum)
        G_int <- sapply(sapply(time, function(x){
            G_0[1:which.min(abs(x - t))]
        }), sum)


        # NEED TO FIGURE OUT THIS INTEGRAL !!!! #
        int.vals <- 1 - 1/(S_int * G_int)# diff(1/(S_int))*1/G_int

        EIF <- S_0_temp[,as.character(t)]*(1 - ((A == 1)*(RCT == 1))/(gRs * g0s) *
                                               (((time <= t)*(event == 1))/( S_times * G_times ) -
                                                    int.vals
                                               )
        )

        return(EIF)
    })


    EIFs_1 <- sapply(fit.times, function(t){

        # get approx value of S(y|x,a) and G(y|x,a) for each observation #
        S_times <- S_1[sapply(time, function(x){which.min(abs(x - fit.times))})]
        G_times <- G_1[sapply(time, function(x){which.min(abs(x - fit.times))})]

        # get approx for integral #
        S_int <- sapply(sapply(time, function(x){
            S_1[1:which.min(abs(x - t))]
        }), sum)
        G_int <- sapply(sapply(time, function(x){
            G_1[1:which.min(abs(x - t))]
        }), sum)

        int.vals <- 1 - 1/(S_int * G_int)# diff(1/(S_int))*1/G_int

        EIF <- S_1_temp[,as.character(t)]*(1 - ((A == 1)*(RCT == 1))/(gRs * g0s) *
                                               (((time <= t)*(event == 1))/( S_times * G_times ) -
                                                    int.vals
                                               )
        )
        return(EIF)
    })

    ## get the efficiency bound by taking the variance at each time point ##
    EB_0 <- apply(EIFs_0, 2, var)
    EB_1 <- apply(EIFs_1, 2, var)

    return(list(S_0 = S_0, S_1 = S_1,
                G_0 = G_0, G_1 = G_1,
                X = X, event.time = event.time, cens.time = cens.time,
                EIFs_0 = EIFs_0, EIFs_1 = EIFs_1, EB_0 = EB_0, EB_1 = EB_1))
}



# N = 10^3
# covs = 5
# lambdaT = 0.2
# beta_T = c(rep(0.5, covs),-1)
# s = 1
# rho = 1
# fit.times = seq(0,10,0.1)
# beta_C = c(rep(-1, covs),-1/5)
# lambdaC = 0.01
# beta_R = rep(0.5, covs)
# beta_A = rep(0.5, covs)
#
# # dat_counterfactual <- generate_counterfactual_nonlinear(N = N, covs = covs, lambdaT = lambdaT, beta_T = beta_T,
# #                                 lambdaC = lambdaC, beta_C = beta_C, beta_R = beta_R, beta_A = beta_A,
# #                                 s = s, rho = rho, fit.times = fit.times)
# #
# # save(dat_counterfactual, file = "./results/counterfactual_data_nonlinear.Rdata")
# library(ggplot2)
# library(reshape2)
# library(dplyr)
# test <- generate_counterfactual_nonlinear(N = N, covs = covs, lambdaT = lambdaT, beta_T = beta_T,
#                                           lambdaC = lambdaC, beta_C = beta_C, beta_R = beta_R, beta_A = beta_A,
#                                           s = s, rho = rho, fit.times = fit.times)
# toy_ex <- as.data.frame(cbind(c(test$S_0, test$S_1), rep(fit.times, 2),
#                               c(rep("S_0", length(fit.times)),
#                                 rep("S_1", length(fit.times)))))
# colnames(toy_ex) <- c("Surv", "t", "type")
# toy_ex %>%
#     mutate(Surv = as.numeric(as.character(Surv)),
#            t = as.numeric(as.character(t))) %>%
#     ggplot(aes(x = t, y = Surv, group = type, color = type)) + geom_line()
