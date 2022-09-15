#' Estimate counterfactual survival functions
#'
#' This function estimates counterfactual survival functions and contrasts from right-censored data subject to potential confounding.
#'
#' @param time \code{n x 1} numeric vector of observed right-censored follow-up times; i.e. the minimum of the event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of whether an event was observed.
#'
#'

summarise_sim <- function(CF_sim_obs, t = NULL){

    # remove failed #
    rm <- c()
    for(i in 1:length(CF_sim_obs)){
        if(!is.null(CF_sim_obs[[i]]$fit$nuisance$prop.pred.RCT))
            if(min(CF_sim_obs[[i]]$fit$nuisance$prop.pred.RCT) == 0)
                rm <- c(rm, i)
    }
    if(length(rm) > 0)
        CF_sim_obs <- CF_sim_obs[-rm]

    out <- lapply(t, function(temp){
        # delta #
        ## Bias ##
        # print(temp)
        bias0 <- mean(sapply(CF_sim_obs, "[[", "Delta_0_bias")[temp,], na.rm = T)
        bias1 <- mean(sapply(CF_sim_obs, "[[", "Delta_1_bias")[temp,], na.rm = T)
        ## MSE ##
        MSE0 <- mean(sapply(CF_sim_obs, "[[", "Delta_0_MSE")[temp,], na.rm = T)
        MSE1 <- mean(sapply(CF_sim_obs, "[[", "Delta_1_MSE")[temp,], na.rm = T)

        # Quantiles #
        ## Bias ##
        quant_bias0 <- apply((sapply(CF_sim_obs, "[[","quant_bias_0")),1, function(x){mean(abs(x), rm.na = T)})
        quant_bias1 <- apply((sapply(CF_sim_obs, "[[","quant_bias_1")),1, function(x){mean(abs(x), rm.na = T)})
        ## MSE ##
        quant_MSE0 <- apply((sapply(CF_sim_obs, "[[","quant_MSE_0")),1, function(x){mean(abs(x), rm.na = T)})
        quant_MSE1 <- apply((sapply(CF_sim_obs, "[[","quant_MSE_1")),1, function(x){mean(abs(x), rm.na = T)})

        return(list(
            "bias0" = bias0,
            "bias1" = bias1,
            "MSE0" = MSE0,
            "MSE1" = MSE1,
            "quant_bias0" = quant_bias0,
            "quant_bias1" = quant_bias1,
            "quant_MSE0" = quant_MSE0,
            "quant_MSE1" = quant_MSE1
        ))
    })
    return(out)
}
