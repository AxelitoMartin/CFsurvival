# # remotes::install_github("AxelitoMartin/CFsurvival")
# rm(list=ls())
# library(CFsurvival)
# library(dplyr)
# library(haven)
# library(tibble)
# library(gtsummary)
# library(fastDummies)
#
# # load WHI data #
# whi_ct <- as_tibble(read_dta("~/Causal/WHI_CT.dta")) %>% mutate(RCT = 1)
# whi_os <- as_tibble(read_dta("~/Causal/WHI_OS.dta")) %>% mutate(RCT = 0)
# keep <- intersect(colnames(whi_ct),colnames(whi_os))
# whi_dat <- rbind(whi_ct %>%  select(keep),
#                  whi_os %>%  select(keep)) %>%
#     column_to_rownames("id")
#
# sum(whi_ct$chd_failure)
# dim(whi_dat)
# head(whi_dat)
# colnames(whi_dat)
#
# apply(whi_dat, 2, function(x) sum(is.na(x)))
# whi_dat <- whi_dat %>%
#     select(-one_of("syst","hyptpill","cardrest","hip55","agefbir"))
# whi_dat_comp <- whi_dat[complete.cases(whi_dat),]
# dim(whi_dat_comp)
#
# # whi_dat %>%
# #   tbl_summary()
#
# # change categorical variables written as numeric factors to dummy variables #
# dums <- c("ethnic")
# tmp <- whi_dat_comp %>%
#     select(dums) %>%
#     mutate_all(as.character) %>%
#     fastDummies::dummy_cols(remove_first_dummy = T) %>%
#     select(-one_of(dums))
# whi_dat_comp <- as.data.frame(cbind(
#     whi_dat_comp %>% select(-one_of(dums)),
#     tmp
# ) %>% mutate_all(as.character) %>%
#     mutate_all(as.numeric)
# )
#
#
# whi_dat_comp %>%
#     tbl_summary(by = RCT) %>%
#     add_p() %>%
#     bold_labels()
#
# # make data for the CF survival run #
#
# # data #
# obs.time <- whi_dat_comp$time_to_chd[whi_dat_comp$RCT == 1]
# obs.event <- whi_dat_comp$chd_failure[whi_dat_comp$RCT == 1]
# rx <- whi_dat_comp$HT_intervention[whi_dat_comp$RCT == 1]
# confounders <- as.data.frame(
#     whi_dat_comp %>%
#         filter(RCT == 1) %>%
#         select(-one_of(c("time_to_chd","chd_failure","HT_intervention","RCT","cvd")))
# )
#
# W_c <- as.data.frame(
#     whi_dat_comp %>%
#         filter(RCT == 0) %>%
#         select(-one_of(c("time_to_chd","chd_failure","HT_intervention","RCT","cvd")))
# )
#
# # data summary #
# summary(obs.time) # lots of observations at the boundary
# summary(obs.event) # low amount of events
# summary(rx) # is random good
# summary(obs.time[rx == 1]) # lots of observations at the boundary
# summary(obs.event[rx == 1]) # low amount of events
# summary(rx) # is random good
#
# # parameters #
# prop.RCT.SL.library = c("SL.glm","SL.rpart","SL.glmnet")
# prop.SL.library = c("SL.mean")
# event.SL.library = c( "survSL.coxph", "survSL.weibreg","survSL.rfsrc")
# cens.SL.library = c("survSL.coxph", "survSL.weibreg","survSL.rfsrc")
# # summary(obs.time)
# fit.times <- seq(0,2150,150)
#
#
# fit <- CFsurvival(time = obs.time, event = obs.event, treat = rx,
#                   confounders =  confounders, contrasts = NULL,
#                   verbose = FALSE, fit.times = fit.times,
#                   fit.treat = c(0,1),
#                   nuisance.options = list(
#                       prop.RCT.SL.library = prop.RCT.SL.library,
#                       prop.SL.library = prop.SL.library,
#                       event.SL.library = event.SL.library,
#                       cens.SL.library = cens.SL.library, verbose = T),
#                   W_c = W_c
# )
#
# library(ggplot2)
# fit$surv.df
# ggplot(fit$surv.df) +
#     # geom_line(aes(time, true.surv, group=trt), color='black') +
#     geom_step(aes(time, surv, color=as.factor(trt), group=trt)) +
#     # geom_step(aes(time, ptwise.lower, color=as.factor(trt), group=trt), linetype=2) +
#     # geom_step(aes(time, ptwise.upper, color=as.factor(trt), group=trt), linetype=2) +
#     geom_step(aes(time, unif.ew.lower, color=as.factor(trt), group=trt), linetype=3) +
#     geom_step(aes(time, unif.ew.upper, color=as.factor(trt), group=trt), linetype=3) +
#     scale_color_discrete("Treatment") +
#     xlab("Time") +
#     ylab("Survival") +
#     coord_cartesian(xlim=c(0,2150), ylim=c(0.98,1))
# save(fit, file = "~/Causal/WHI_fit.Rdata")
