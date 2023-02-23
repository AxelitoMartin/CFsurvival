#' Estimate counterfactual survival functions
#'
#' This function estimates counterfactual survival functions and contrasts from right-censored data subject to potential confounding.
#'
#' @param data \code{n x 1} numeric vector of observed right-censored follow-up times; i.e. the minimum of the event and censoring times.
#' @export
#'

summarise_sim <- function(data, get_cov = F, get_EB = F, N){


    library(dplyr)
    # get the different scenarios of interest #
    scenarios <- names(data[[1]])
    iter <- length(data)
    # N <- nrow(data[[1]][["correct"]][["dat"]])
    time_points <- length(data[[1]][["correct"]][["bias"]])
    bias <- RMSE <-
        as.data.frame(matrix(0L, nrow = length(scenarios), ncol = time_points))
    quants_0 <- quants_1 <- as.data.frame(matrix(0L, nrow = length(scenarios), ncol = length(data[[1]][["correct"]][["quant_MSE0"]])))
    rownames(bias) <- rownames(RMSE) <- rownames(quants_0) <- rownames(quants_1) <-  scenarios
    colnames(bias) <- colnames(RMSE) <- paste0("t",1:time_points)
    colnames(quants_0) <- colnames(quants_1) <- names(data[[1]][["correct"]][["quant_MSE0"]])
    coverage <- list()
    EB <- list()

    for(i in scenarios){

        # bias #
        bias_mean <- cumsum(abs(apply(sapply(data, function(x){
            return(x[[i]][["bias"]])
        }), 1, function(y){mean(y)})))
        bias[i,] <- bias_mean
        # bias0 <- cumsum(abs(apply(sapply(data, function(x){
        #     return(x[[i]][["bias0"]])
        # }), 1, function(y){mean(y)})))
        # bias1 <- cumsum(abs(apply(sapply(data, function(x){
        #     return(x[[i]][["bias1"]])
        # }), 1, function(y){mean(y)})))
        # # mean_bias <- apply(rbind(bias0,bias1), 2, mean)
        # bias_0[i,] <- bias0
        # bias_1[i,] <- bias1

        # RMSE #
        RMSE0 <- cumsum(sqrt(apply(sapply(data, function(x){
            return(x[[i]][["bias"]]^2)
        }),1,mean)))
        RMSE[i,] <- RMSE0
        # RMSE0 <- cumsum(sqrt(apply(sapply(data, function(x){
        #     return(x[[i]][["bias0"]])
        # })^2,1,mean)))
        # RMSE1 <- cumsum(sqrt(apply(sapply(data, function(x){
        #     return(x[[i]][["bias1"]])
        # })^2,1,mean)))
        # # mean_RMSE <- apply(rbind(RMSE0, RMSE1), 2, mean)
        # RMSE_0[i,] <- RMSE0
        # RMSE_1[i,] <- RMSE1

        # quantiles #
        quant0 <- apply(sapply(data, function(x){
            return(abs(x[[i]][["quant_bias0"]]))
        }), 1, mean)
        quant1 <- apply(sapply(data, function(x){
            return(abs(x[[i]][["quant_bias1"]]))
        }), 1, mean)
        # mean_quant <- apply(rbind(quant0, quant1), 2, mean)
        # quants[i,] <- mean_quant
        quants_0[i,] <- quant0
        quants_1[i,] <- quant1

        # Add coverage here #

        # Coverage counts how many times the true counterfactual value
        # falls within the bounds of our 95% constructed CI

        if(get_cov)
            coverage[[i]] <- as_tibble(do.call('rbind',lapply(data, function(x){

                temp <- x[[i]]$fit.surv

                # temp <- data[[1]][["correct"]]$fit.surv

                temp_bool <- temp %>%
                    select(time,trt, true.surv, ptwise.lower, ptwise.upper) %>%
                    # group_by(trt, time) %>%
                    # rowwise() %>%
                    mutate(coverage = ptwise.lower <= true.surv & ptwise.upper >= true.surv) %>%
                    select(time,trt, coverage)

                return(temp_bool)

            }))) %>%
            group_by(trt, time) %>%
            summarise(Cov = sum(coverage)/iter) %>%
            ungroup() %>%
            mutate(scenario = i)

        else
            coverage <- NULL


        #######################
        if(get_EB){

            # EB_0 <- apply(do.call('rbind',lapply(data, function(x){
            #
            #     temp_0 <- x[[i]]$EIF0
            #     return(apply(temp_0,2,var))
            # })), 2, mean)



           EB_0 <- apply(do.call('rbind',lapply(data, function(x){

               temp_0 <- x[[i]]$EIF0
           })), 2, var)

           EB_1 <- apply(do.call('rbind',lapply(data, function(x){

               temp_0 <- x[[i]]$EIF1
           })), 2, var)

           EB[[i]] <- list(EB_0, EB_1)
        }
    }



    ####################

    return(list(
        "bias" = bias,
        "RMSE" = RMSE,
        # "bias0" = bias_0,
        # "bias1" = bias_1,
        # "RMSE0" = RMSE_0,
        # "RMSE1" = RMSE_1,
        "quants0" = quants_0,
        "quants1" = quants_1,
        "coverage" = coverage,
        "EB" = EB
    ))
}


# load("./results/DML_1000_results_test.Rdata")
# DML_500 <- DML_results
# test_DML <- summarise_sim(data = DML_500, get_cov = T, get_EB = T)
# test_DML$bias
# test_DML$RMSE
# test_DML$coverage[["correct"]]
# test_DML$coverage[["sampling"]]
# test_DML$coverage[["survival"]]
# test_DML$coverage[["incorrect"]]


# #
# load("./results/Naive_500_results.Rdata")
# Naive_500 <- Naive_results
# test_Naive <- summarise_sim(data = Naive_500)
# test_Naive$bias0
# test_Naive$bias1
# test_Naive$RMSE0
# test_Naive$RMSE1

#
#
# load("./results/DML_2000_results.Rdata")
# DML_2000 <- DML_results
# test_DML <- summarise_sim(data = DML_2000)
# test_DML$bias
#
# load("./results/Naive_2000_results.Rdata")
# Naive_2000 <- Naive_results
# test_Naive <- summarise_sim(data = Naive_2000)
# test_Naive$bias

########################################

#
# summarise_sim <- function(CF_sim_obs, t = NULL){
#
#     # remove failed #
#     rm <- c()
#     for(i in 1:length(CF_sim_obs)){
#         if(!is.null(CF_sim_obs[[i]]$fit$nuisance$prop.pred.RCT))
#             if(min(CF_sim_obs[[i]]$fit$nuisance$prop.pred.RCT) == 0)
#                 rm <- c(rm, i)
#     }
#     if(length(rm) > 0)
#         CF_sim_obs <- CF_sim_obs[-rm]
#
#     out <- lapply(t, function(temp){
#         # delta #
#         ## Bias ##
#         # print(temp)
#         bias0 <- mean(sapply(CF_sim_obs, "[[", "Delta_0_bias")[temp,], na.rm = T)
#         bias1 <- mean(sapply(CF_sim_obs, "[[", "Delta_1_bias")[temp,], na.rm = T)
#         ## MSE ##
#         MSE0 <- mean(sapply(CF_sim_obs, "[[", "Delta_0_MSE")[temp,], na.rm = T)
#         MSE1 <- mean(sapply(CF_sim_obs, "[[", "Delta_1_MSE")[temp,], na.rm = T)
#
#         # Quantiles #
#         ## Bias ##
#         quant_bias0 <- apply((sapply(CF_sim_obs, "[[","quant_bias_0")),1, function(x){mean(abs(x), rm.na = T)})
#         quant_bias1 <- apply((sapply(CF_sim_obs, "[[","quant_bias_1")),1, function(x){mean(abs(x), rm.na = T)})
#         ## MSE ##
#         quant_MSE0 <- apply((sapply(CF_sim_obs, "[[","quant_MSE_0")),1, function(x){mean(abs(x), rm.na = T)})
#         quant_MSE1 <- apply((sapply(CF_sim_obs, "[[","quant_MSE_1")),1, function(x){mean(abs(x), rm.na = T)})
#
#         return(list(
#             "bias0" = bias0,
#             "bias1" = bias1,
#             "MSE0" = MSE0,
#             "MSE1" = MSE1,
#             "quant_bias0" = quant_bias0,
#             "quant_bias1" = quant_bias1,
#             "quant_MSE0" = quant_MSE0,
#             "quant_MSE1" = quant_MSE1
#         ))
#     })
#     return(out)
# }



