rm(list=ls())

if(!require("pacman")) install.packages("pacman")
pacman::p_load(MASS, dplyr, ggplot2, truncnorm)

CreateContConfounder <- function(d0, d1, var) {
    sapply(var, function(T) {
      return(ifelse(T == 0, rnorm(n = 1, mean = d0, sd = 1), rnorm(n=1, mean=d1, sd=1)))
    })
}

createRandomModel <- function(
  N = 100,
  total_confounders = 10,
  probT = 0.5,
  probKnown = 0.5,
  true_effect_size = 0,
  error_mean = 0,
  error_sd = 5,
  lo_mean=0.3,
  hi_mean=0.7,
  lo_es=0.5,
  hi_es=5,
  prob_pos_confounder = 0.6
){

  set.seed(NULL)
  dat <- data.frame(T = rbinom(n = N, size = 1, prob = probT))

  for(i in 1:(total_confounders / 2)){

    k <- ifelse(runif(1,0,1) < probKnown, "K", "U")

    p1 <- rtruncnorm(n=1, a=0, b=1, mean=lo_mean, sd=1)
    p2 <- rtruncnorm(n=1, a=0, b=1, mean=hi_mean, sd=1)

    p3 <- rnorm(n=1, mean=lo_mean, sd=1)
    p4 <- rnorm(n=1, mean=hi_mean, sd=1)

    dat[, paste0(k, "B", i)] <- ifelse(
      dat$T == 0,
      rbinom(n = nrow(dat), size = 1, prob = p1),
      rbinom(n = nrow(dat), size = 1, prob = p2)
    )
    dat[, paste0(k, "C", i)] <- CreateContConfounder(p3, p4, dat$T)
  }


  params <- c(true_effect_size, (runif(ncol(dat) - 1, lo_es / 10, hi_es / 10) * 10) * ifelse(runif(1,0,1)<prob_pos_confounder, 1, -1))
  dat$Y <- (dat %>% as.matrix()) %*% (params %>% as.matrix()) + rnorm(nrow(dat), mean=error_mean, sd=error_sd)
  dat$Y <- dat$Y / sd(dat$Y)

  total_model <- lm(Y ~ T + ., data = dat)
  known_model <- lm(Y ~ T + ., data = dat %>% dplyr::select(Y, T, starts_with("K")))

  sum_total_model <- summary(total_model)
  sum_known_model <- summary(known_model)
  ci_total_model <- confint(total_model)
  ci_known_model <- confint(known_model)

  total_b <- sum_total_model$coefficients[2,1]
  total_lo95 <- ci_total_model[2, 2]
  total_hi95 <- ci_total_model[2, 1]
  total_p <- sum_total_model$coefficients[2,4]
  known_b <- sum_known_model$coefficients[2,1]
  known_lo95 <- ci_known_model[2, 2]
  known_hi95 <- ci_known_model[2, 1]
  known_p <- sum_known_model$coefficients[2,4]

  return(rbind(c(
    "total_b" = total_b,
    "total_lo95" = total_lo95,
    "total_hi95" = total_hi95,
    "total_p" = total_p,
    "known_b" = known_b,
    "known_lo95" = known_lo95,
    "known_hi95" = known_hi95,
    "known_p" = known_p,
    "known_confounders" = dat %>% dplyr::select(starts_with("K")) %>% ncol()
  )))

}

es <- 0
total_confounders <- 20
probKnown <- 0.3
lo_mean = 0.4
hi_mean = 0.6
lo_es <- 0.5
hi_es <- 5
N <- 1000
prob_pos_confounder <- 0.5
probT <- 0.5

sims <- replicate(
  n = 100,
  createRandomModel(
    N=N,
    total_confounders = total_confounders,
    probKnown = probKnown,
    probT = probT,
    true_effect_size = es,
    error_sd = 10,
    lo_mean=lo_mean,
    hi_mean=hi_mean,
    lo_es=lo_es,
    hi_es=hi_es,
    prob_pos_confounder = prob_pos_confounder

  ),
  simplify = F
)

sims_df <- do.call(rbind.data.frame, sims)
sims_df$run <- 1:nrow(sims_df)
sims_df$rank <- rank(sims_df$known_b)

sims_df_long <- reshape(sims_df %>% dplyr::select(-known_confounders), direction='long',
        varying=sims_df %>% dplyr::select(-run, -rank, -known_confounders) %>% names(),
        timevar='var',
        times=c('total', 'known'),
        v.names=c('b', 'lo95', 'hi95', 'p'),
        idvar=c('run', "rank")) %>% tibble()

sims_df_long$group <- as.vector(as.character(sims_df_long$var))
sims_df_long <- sims_df_long %>% mutate(
  group = case_when(
    group == "known" ~ "Only Known Confounders Included\n(Observational Methods)",
    group == "total" ~ "All Confounders Included\n(Randomized Methods)"
  ),
  group = factor(group,
    levels=c("Only Known Confounders Included\n(Observational Methods)", "All Confounders Included\n(Randomized Methods)")
  ),
  col = case_when(
    es >= lo95 & es <= hi95 ~ "Accurate",
    T ~ "Inaccurate"
  )
)

cols <- c("darkgreen", "darkred")
p <- ggplot(data = sims_df_long, aes(x=b, y=rank, xmin=lo95, xmax=hi95)) +
        geom_pointrange(aes(colour = col)) +
        geom_vline(xintercept=es) +
        facet_wrap(~group) +
        labs(
          x = "Estimated Effect Size (95%CI)",
          y = "Simulation",
          title = sprintf("Simulated Studies (N=%s) with Known and Unknown Confounders", N),
          subtitle=sprintf("%s Confounders (%s%% Known) with Group Means of %s and %s and ES ~+/-Unif(%s, %s)", total_confounders, probKnown*100, lo_mean, hi_mean, lo_es, hi_es),
          col = "Study Accuracy"
        ) +
        theme_bw() +
        scale_color_manual(values=cols, labels=c("Accurate", "Inaccurate"))

ggsave(sprintf("confounder-simulation-%sconfounders-%sknown.png", total_confounders, probKnown*100), p, width=8, height=8)
