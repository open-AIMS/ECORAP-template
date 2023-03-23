## ---- loadFunctions
source('../R/functions.R')
## ----end

## Coral Growth *****************************************************************

## ---- CoralGrowth readData
data <- read_csv(paste0(DATA_PATH,
                        "primary/All_size_2021-2022_ML.csv"),
                 trim_ws=TRUE)
## ----end

## ---- CoralGrowth glimpse
data %>% glimpse() 
## ----end

## ---- CoralGrowth report table
data %>% report::report_table() 
## ----end

## Process data =================================================================
## ---- CoralGrowth process
data <- 
    data %>% 
    ## Replace Acyt and Ahya with Atab in all Class* fields
    mutate(across(matches('Class_.*'),
                  ~ forcats::fct_recode(.x,
                                        'Atab' = 'Acyt',
                                        'Atab' = 'Ahya'))) %>%
    ## Exclude colonies that were dead or missing in 2022
    filter(!Note_2022 %in% c(",,,d", ",,,m")) %>%
    ## Only keep colonies that have good images in 2021 and 2022
    filter(to_use_for_growth == "yes") %>%
    ## Drop the extraneous columns
    dplyr::select(-AreaChange, -to_use_for_growth, -to_use_2021, -Survival,
                  -matches("Note_.*")) %>%
    ## Add a colony ID
    mutate(Colony_ID = factor(1:n())) %>%
    ## Pivote longer
    pivot_longer(cols = matches(".*_[0-9]{4}"),#c(Class_2021, Class_2022, ID_2021, ID_2022),
    ## pivot_longer(cols = c(Class_2021, Class_2022),
                 names_to = c(".value","Year"),
                 names_pattern = "(.+)_(.+)",
                 ## names_sep = "_",
                 values_to = 'Class') %>%
    ## Express area in cm2
    ## mutate(Area = Area * 10000) %>%
    ## Define unique hierarchical levels
    mutate(
      Reef_Site = factor(paste(Reef, Site, sep = '_')),
      Reef_Site_Zone = factor(paste(Reef_Site, Zone, sep = '_')),
      Plot = str_replace_all(Plot_ID, "[^_]*_", ""),
      Reef_Site_Zone_Plot = factor(paste(Reef_Site_Zone, Plot, sep = '_')),
      Colony = factor(paste(Plot, Colony_ID, sep = "_"))
      ) %>%
    ## Define Site_code as each unique Reef/Site/Zone combination
    mutate(Site_code = factor(paste(Reef, Site, Zone, sep = '_'))) %>%
    mutate(Site_type = factor(str_sub(Site, 1, 2)))

data %>% glimpse()
saveRDS(data, file = paste0(DATA_PATH, "processed/data1.RData"))
## ----end

## ---- CoralGrowth process tests
data <- readRDS(file = paste0(DATA_PATH, "processed/data1.RData"))
##1. There should be exactly two records per colony,
##   highlight those for which this is not the case
data %>%
    group_by(Plot_ID, Colony_ID, Class) %>%
    count() %>%
    filter(n != 2) -> tests
tests

## Address this
data <- data %>%
    filter(!Colony_ID %in% unique(tests$Colony_ID))
saveRDS(data, file = paste0(DATA_PATH, "processed/data.RData"))
## ----end

## EDA ==========================================================================
## ---- CoralGrowth EDA
g <- data %>%
    ggplot(aes(y = Area, x = Site_Type, color = Year)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.02, dodge.width = 0.9),
               alpha = 0.1) +
    geom_violin(fill = NA) +
    facet_wrap(~Class, scales='free_y') +
    theme_bw()

ggsave(filename = paste0(FIGS_PATH, "/EDA1.png"),
       g,
       width = 10,
       height = 5,
       dpi = 300)
## ----end

## ---- CoralGrowth EDA Plots
g <- data %>%
  filter(Class == "Atab") %>%
  ggplot(aes(y = Area, x = Year, group = Colony, colour = Zone)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  facet_grid(Zone ~ Site_Type, scales='free') +
  theme_bw()
ggsave(filename = paste0(FIGS_PATH, "/EDA3.png"),
       g,
       width = 10,
       height = 6,
       dpi = 300)
## ----end

## ---- CoralGrowth processing further
data.sub <- data %>% 
    filter(Class == 'Atab') %>%
    droplevels() %>% 
    mutate(Treatment = paste(Site_Type, Zone, sep = '_'))


## Note that the first level has only a single observation - not a
## good choice as the intercept
data.sub %>% pull(Treatment) %>% levels()
data.sub %>% pull(Year) %>% unique()

data.sub <- data.sub %>%
  mutate(Treatment = forcats::fct_relevel(Treatment, 'BA_S')) 
saveRDS(data.sub, file = paste0(DATA_PATH, "processed/data.sub.RData"))
## ----end


## Modelling ====================================================================

## ---- CoralGrowth brms Model 1.1 priors
data.sub <- readRDS(file = paste0(DATA_PATH, "processed/data.sub.RData"))
model <- 'mod1.1'
data.sub %>%
    group_by(Reef) %>%
    summarise(across(Area,  list(Mean = ~ mean(log(.x)),
                                 Median = ~ median(log(.x)),
                                 SD = ~ sd(log(.x)),
                                 MAD = ~ mad(log(.x)),
                                 N = ~ length(.x)))) 
data.sub %>%
    group_by(Reef, Site, Plot, Colony) %>%
    summarise(across(Area,  list(Mean = ~ mean(log(.x)),
                                 Median = ~ median(log(.x)),
                                 SD = ~ sd(log(.x)),
                                 MAD = ~ mad(log(.x)),
                                 N = ~ length(.x)))) %>%
    ungroup() %>% 
    summarise(SD = sd(Area_SD, na.rm = TRUE))
    
data.sub %>%
    summarise(across(Area,  list(Mean = ~ mean(log(.x)),
                                 Median = ~ median(log(.x)),
                                 SD = ~ sd(log(.x)),
                                 MAD = ~ mad(log(.x)),
                                 N = ~ length(.x))))
form <- bf(log(Area) ~ 1 + 
               (1 | Reef/Site/Zone/Plot/Colony),
           family = gaussian())
get_prior(form, data = data.sub)
## ----end
## ---- CoralGrowth brms Model 1.1 fit
priors <- set_prior("student_t(3, -2, 1)", class = "Intercept") +
    set_prior("student_t(3, 0, 1)", class = "sd") +
    set_prior("student_t(3, 0, 1)", class = "sigma")

mod <- brm(form,
           data = data.sub,
           prior = priors,
           sample_prior = "yes",
           cores = 3,
           chains = 3,
           warmup = 1000,
           thin = 5,
           iter = 3000,
           control = list(adapt_delta = 0.99, max_treedepth = 17),
           backend = 'cmdstanr'
           )
save(priors, mod, file = paste0(DATA_PATH, "modelled/",model,".RData"))
## ----end
## ---- CoralGrowth brms Model 1.1 posterior_prior
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
g <- mod %>% SUYR_prior_and_posterior() +
    theme_classic() +
    theme(legend.position = c(1,0), legend.justification = c(1,0))
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Posterior_prior_',model,'.png'),
       g,
       width = 10,
       height = 8,
       dpi = 100
       )
## ----end
## ---- CoralGrowth brms Model 1.1 MCMC diagnostics
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Trace_',model,'.png'),
       rstan::stan_trace(mod$fit),
       width = 12,
       height = 8,
       dpi = 300
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ac_',model,'.png'),
       rstan::stan_ac(mod$fit),
       width = 12,
       height = 8,
       dpi = 300
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Rhat_',model,'.png'),
       rstan::stan_rhat(mod$fit),
       width = 12,
       height = 8,
       dpi = 300
       )
## ----end
## ---- CoralGrowth brms Model 1.1 ppc
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc1_',model,'.png'),
       mod %>% pp_check(type = "dens_overlay", ndraws = 100),
       width = 5,
       height = 4,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc3_',model,'.png'),
       mod %>% pp_check(type = 'loo_pit_qq'),
       width = 5,
       height = 4,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc4_',model,'.png'),
       mod %>% pp_check(type = 'loo_pit_overlay'),
       width = 5,
       height = 4,
       dpi = 100
       )
## ----end
## ---- CoralGrowth brms Model 1.1 DHARMa
brms_dharma_res <- make_brms_dharma_res(
    mod, integerResponse = FALSE
) 
ggsave(filename = paste0(OUTPUT_PATH, 'figures/DHARMa1_',model,'.png'),
       wrap_elements(~testUniformity(brms_dharma_res)) +
       wrap_elements(~plotResiduals(brms_dharma_res, form = factor(rep(1, nrow(data.sub))))) +
       wrap_elements(~testDispersion(brms_dharma_res)),
       width = 12,
       height = 4,
       dpi = 100)
## ----end

## ---- CoralGrowth brms Model 3.1 priors
data.sub <- readRDS(file = paste0(DATA_PATH, "processed/data.sub.RData"))
model <- 'mod3.1'
data.sub %>%
    group_by(Site_Type, Zone) %>%
    summarise(across(Area,  list(Mean = ~ mean(log(.x)),
                                 Median = ~ median(log(.x)),
                                 SD = ~ sd(log(.x)),
                                 MAD = ~ mad(log(.x)),
                                 N = ~ length(.x)))) 

data.sub <- readRDS(file = paste0(DATA_PATH, "processed/data.sub.RData"))

data.sub1 <- data.sub %>%
    mutate(Zone = factor(Zone, levels = c('S', 'D')))
form <- bf(log(Area) ~ Site_Type*Zone + 
               (1 | Reef/Site/Zone/Plot/Colony),
           family = gaussian())
get_prior(form, data = data.sub1)
## ----end
## ---- CoralGrowth brms Model 3.1 fit
priors <- set_prior("student_t(3, -2, 1)", class = "Intercept") +
    set_prior("student_t(3, 0, 1)", class = "b") +
    set_prior("student_t(3, 0, 1)", class = "sd") +
    set_prior("student_t(3, 0, 1)", class = "sigma")

mod <- brm(form,
           data = data.sub1,
           prior = priors,
           sample_prior = "yes",
           cores = 3,
           chains = 3,
           warmup = 1000,
           thin = 5,
           iter = 3000,
           control = list(adapt_delta = 0.99, max_treedepth = 17),
           backend = 'cmdstanr'
           )
save(priors, mod, file = paste0(DATA_PATH, "modelled/",model,".RData"))
## ----end
## ---- CoralGrowth brms Model 3.1 posterior_prior
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
g <- mod %>% SUYR_prior_and_posterior() +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal")
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Posterior_prior_',model,'.png'),
       g,
       width = 10,
       height = 8,
       dpi = 100
       )
## ----end
## ---- CoralGrowth brms Model 3.1 MCMC diagnostics
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Trace_',model,'.png'),
       rstan::stan_trace(mod$fit),
       width = 12,
       height = 8,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ac_',model,'.png'),
       rstan::stan_ac(mod$fit),
       width = 12,
       height = 8,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Rhat_',model,'.png'),
       rstan::stan_rhat(mod$fit),
       width = 12,
       height = 8,
       dpi = 100
       )
## ----end
## ---- CoralGrowth brms Model 3.1 ppc
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc1_',model,'.png'),
       mod %>% pp_check(type = "dens_overlay", ndraws = 100),
       width = 5,
       height = 4,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc2_',model,'.png'),
       mod %>% pp_check(group = 'Site_Type', type = 'intervals_grouped'),
       width = 5,
       height = 4,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc3_',model,'.png'),
       mod %>% pp_check(type = 'loo_pit_qq'),
       width = 5,
       height = 4,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc4_',model,'.png'),
       mod %>% pp_check(type = 'loo_pit_overlay'),
       width = 5,
       height = 4,
       dpi = 100
       )
## ----end
## ---- CoralGrowth brms Model 3.1 DHARMa
brms_dharma_res <- make_brms_dharma_res(
    mod, integerResponse = FALSE
) 
ggsave(filename = paste0(OUTPUT_PATH, 'figures/DHARMa1_',model,'.png'),
       wrap_elements(~testUniformity(brms_dharma_res)) +
       wrap_elements(~plotResiduals(brms_dharma_res, form = factor(rep(1, nrow(data.sub))))) +
       wrap_elements(~testDispersion(brms_dharma_res)),
       width = 12,
       height = 4,
       dpi = 100)
ggsave(filename = paste0(OUTPUT_PATH, 'figures/DHARMa2_',model,'.png'),
       wrap_elements(~testUniformity(brms_dharma_res)) +
       wrap_elements(~plotResiduals(brms_dharma_res)) +
       wrap_elements(~testDispersion(brms_dharma_res)),
       width = 12,
       height = 4,
       dpi = 100)
## ----end
## ---- CoralGrowth brms Model 3.1 Partial plots derived means
em.sum <- mod %>% emmeans(~Site_Type*Zone) %>%
    tidybayes::tidy_draws() %>%
    exp() %>%
    posterior::summarise_draws(median,
                               HDInterval::hdi,
                               rhat,
                               ess_bulk,
                               ess_tail,
                               length) %>%
    mutate(variable = str_replace_all(variable, "Site_Type |Zone ", "")) %>%
    separate(variable, into = c('Site_Type', 'Zone'), sep = ', ')
## ----end
## ---- CoralGrowth brms Model 3.1 Partial plots derived means figure
g <- 
    em.sum %>%  
    ggplot(aes(y = median, x = Site_Type, colour = Zone)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.2)) +
    scale_y_log10() +
    theme_classic()
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Partial.type1_',model,'.png'),
       g, 
       width = 6, height = 6,
       dpi = 100
       )
## OR
g <- mod %>% emmeans(~Site_Type*Zone) %>%
    tidybayes::gather_emmeans_draws() %>%
    mutate(.value = exp(.value)) %>%
    ggplot(aes(x = .value)) +
    ggridges::geom_density_ridges_gradient(aes(y = interaction(Site_Type,Zone),
                                               fill = stat(quantile)),
                                           ## colour = "orange",
                                           quantile_lines = TRUE,
                                           quantiles = c(0.025, 0.05, 0.95, 0.975),
                                           show.legend = FALSE) +
    ## scale_fill_gradient(low = "white", high = "orange") +
    scale_fill_manual(values = c('#ffa600', '#7a5195', '#003f5c', '#7a5195', '#ffa600')) +
    scale_y_discrete('Intercept') + 
    scale_x_log10('Area (cm²)') +
    theme_classic()
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Partial.type_',model,'.png'),
       g, 
       width = 6, height = 6,
       dpi = 100
       )
## ----end

## ---- CoralGrowth brms Model 5.1 priors
data.sub <- readRDS(file = paste0(DATA_PATH, "processed/data.sub.RData"))
model <- 'mod5.1'
data.sub %>%
    group_by(Treatment, Year) %>%
    summarise(across(Area,  list(Mean = ~ mean(log(.x)),
                                 Median = ~ median(log(.x)),
                                 SD = ~ sd(log(.x)),
                                 MAD = ~ mad(log(.x)),
                                 N = ~ length(.x)))) 
form <- bf(log(Area) ~ Treatment*Year + 
               (1 | Reef/Site/Zone/Plot/Colony),
           family = gaussian())
get_prior(form, data = data.sub)
## ----end
## ---- CoralGrowth brms Model 5.1 fit
priors <- set_prior("student_t(3, -2, 1)", class = "Intercept") +
    set_prior("student_t(3, 0, 1)", class = "b") +
    set_prior("student_t(3, 0, 1)", class = "sd") +
    set_prior("student_t(3, 0, 1)", class = "sigma")

mod <- brm(form,
           data = data.sub,
           prior = priors,
           sample_prior = "yes",
           cores = 3,
           chains = 3,
           warmup = 1000,
           thin = 5,
           iter = 3000,
           control = list(adapt_delta = 0.99, max_treedepth = 17),
           backend = 'cmdstanr'
           )
save(priors, mod, file = paste0(DATA_PATH, "modelled/",model,".RData"))
## ----end
## ---- brms Model 5.1 posterior_prior
g <- mod %>% SUYR_prior_and_posterior() +
    theme_classic() +
    theme(legend.position = "bottom")
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Posterior_prior_',model,'.png'),
       g,
       width = 10,
       height = 10,
       dpi = 100
       )
## ----end
## ---- brms Model 5.1 MCMC diagnostics
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Trace_',model,'.png'),
       rstan::stan_trace(mod$fit),
       width = 12,
       height = 8,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ac_',model,'.png'),
       rstan::stan_ac(mod$fit),
       width = 12,
       height = 8,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Rhat_',model,'.png'),
       rstan::stan_rhat(mod$fit),
       width = 12,
       height = 8,
       dpi = 100
       )
## ----end
## ---- brms Model 5.1 ppc
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc1_',model,'.png'),
       mod %>% pp_check(type = "dens_overlay", ndraws = 100),
       width = 5,
       height = 4,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc2_',model,'.png'),
       mod %>% pp_check(group = 'Treatment', type = 'intervals_grouped'),
       width = 5,
       height = 4,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc3_',model,'.png'),
       mod %>% pp_check(type = 'loo_pit_qq'),
       width = 5,
       height = 4,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc4_',model,'.png'),
       mod %>% pp_check(type = 'loo_pit_overlay'),
       width = 5,
       height = 4,
       dpi = 100
       )
## ----end
## ---- brms Model 5.1 DHARMa
brms_dharma_res <- make_brms_dharma_res(
    mod, integerResponse = FALSE
) 
ggsave(filename = paste0(OUTPUT_PATH, 'figures/DHARMa1_',model,'.png'),
       wrap_elements(~testUniformity(brms_dharma_res)) +
       wrap_elements(~plotResiduals(brms_dharma_res)) +
       wrap_elements(~testDispersion(brms_dharma_res)),
       width = 12,
       height = 4,
       dpi = 100)
## ----end

## ---- CoralGrowth brms Model 6.1 priors
data.sub <- readRDS(file = paste0(DATA_PATH, "processed/data.sub.RData"))
model <- 'mod6.1'
data.sub1 <- data.sub %>%
    filter(Reef != 'PAOR') %>% droplevels()

data.sub1 %>%
    group_by(Reef) %>%
    summarise(across(Area,  list(Mean = ~ mean(log(.x)),
                                 Median = ~ median(log(.x)),
                                 SD = ~ sd(log(.x)),
                                 MAD = ~ mad(log(.x)),
                                 N = ~ length(.x)))) 
data.sub1 %>%
    group_by(Reef, Site, Plot, Colony) %>%
    summarise(across(Area,  list(Mean = ~ mean(log(.x)),
                                 Median = ~ median(log(.x)),
                                 SD = ~ sd(log(.x)),
                                 MAD = ~ mad(log(.x)),
                                 N = ~ length(.x)))) %>%
    ungroup() %>% 
    summarise(SD = sd(Area_SD, na.rm = TRUE))
    
data.sub1 %>%
    summarise(across(Area,  list(Mean = ~ mean(log(.x)),
                                 Median = ~ median(log(.x)),
                                 SD = ~ sd(log(.x)),
                                 MAD = ~ mad(log(.x)),
                                 N = ~ length(.x))))
data.sub1 %>%
    filter(Reef != 'PAOR') %>% droplevels() %>% 
    group_by(Treatment, Year) %>%
    summarise(across(Area,  list(Mean = ~ mean(log(.x)),
                                 Median = ~ median(log(.x)),
                                 SD = ~ sd(log(.x)),
                                 MAD = ~ mad(log(.x)),
                                 N = ~ length(.x)))) 
form <- bf(log(Area) ~ Treatment*Year + 
               (1 | Reef/Site/Zone/Plot/Colony),
           family = gaussian())
get_prior(form, data = data.sub1)
## ----end
## ---- CoralGrowth brms Model 6.1 fit
priors <- set_prior("student_t(3, -2, 1)", class = "Intercept") +
    set_prior("student_t(3, 0, 1)", class = "b") +
    set_prior("student_t(3, 0, 1)", class = "sd") +
    set_prior("student_t(3, 0, 0.5)", class = "sigma")

mod <- brm(form,
           data = data.sub1,
           prior = priors,
           sample_prior = "yes",
           cores = 3,
           chains = 3,
           warmup = 1000,
           thin = 5,
           iter = 3000,
           control = list(adapt_delta = 0.99, max_treedepth = 17),
           backend = 'cmdstanr'
           )
save(priors, mod, file = paste0(DATA_PATH, "modelled/",model,".RData"))
## ----end
## ---- CoralGrowth brms Model 6.1 posterior_prior
g <- mod %>% SUYR_prior_and_posterior() +
    theme_classic() +
    theme(legend.position = "bottom")
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Posterior_prior_',model,'.png'),
       g,
       width = 10,
       height = 10,
       dpi = 100
       )
## ----end
## ---- CoralGrowth brms Model 6.1 MCMC diagnostics
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Trace_',model,'.png'),
       rstan::stan_trace(mod$fit),
       width = 12,
       height = 8,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ac_',model,'.png'),
       rstan::stan_ac(mod$fit),
       width = 12,
       height = 8,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Rhat_',model,'.png'),
       rstan::stan_rhat(mod$fit),
       width = 12,
       height = 8,
       dpi = 100
       )
## ----end
## ---- CoralGrowth brms Model 6.1 ppc
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc1_',model,'.png'),
       mod %>% pp_check(type = "dens_overlay", ndraws = 100),
       width = 5,
       height = 4,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc2_',model,'.png'),
       mod %>% pp_check(group = 'Treatment', type = 'intervals_grouped'),
       width = 5,
       height = 4,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc3_',model,'.png'),
       mod %>% pp_check(type = 'loo_pit_qq'),
       width = 5,
       height = 4,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc4_',model,'.png'),
       mod %>% pp_check(type = 'loo_pit_overlay'),
       width = 5,
       height = 4,
       dpi = 100
       )
## ----end
## ---- CoralGrowth brms Model 6.1 DHARMa
brms_dharma_res <- make_brms_dharma_res(
    mod, integerResponse = FALSE
) 
ggsave(filename = paste0(OUTPUT_PATH, 'figures/DHARMa1_',model,'.png'),
       wrap_elements(~testUniformity(brms_dharma_res)) +
       wrap_elements(~plotResiduals(brms_dharma_res)) +
       wrap_elements(~testDispersion(brms_dharma_res)),
       width = 12,
       height = 4,
       dpi = 100)
## ----end
## ---- CoralGrowth brms Model 6.1 Partial plots derived means
em.sum <- mod %>% emmeans(~Treatment|Year) %>%
    tidybayes::tidy_draws() %>%
    exp() %>%
    posterior::summarise_draws(median,
                               HDInterval::hdi,
                               rhat,
                               ess_bulk,
                               ess_tail,
                               length) %>%
    mutate(variable = str_replace_all(variable, "Treatment ", "")) %>%
    mutate(variable = str_replace_all(variable, ", Year ", "_")) %>%
    separate(variable, into = c('Site_Type', 'Zone', 'Year'), sep = '_')
## ----end
## ---- CoralGrowth brms Model 6.1 Partial plots derived means figure
g <- 
    em.sum %>%  
    ggplot(aes(y = median, x = Site_Type, colour = Year)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.2)) +
    scale_y_log10() +
    facet_grid(~Zone) + 
    theme_classic()
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Partial.type1_',model,'.png'),
       g, 
       width = 6, height = 6,
       dpi = 100
       )
## OR
g <- mod %>% emmeans(~Treatment|Year) %>%
    tidybayes::gather_emmeans_draws() %>%
    mutate(.value = exp(.value)) %>%
    separate(Treatment, into = c('Site_Type', 'Zone'), sep = '_') %>% 
    ggplot(aes(x = .value)) +
    ggridges::geom_density_ridges_gradient(aes(y = Site_Type,
                                               fill = stat(quantile)),
                                           ## colour = "orange",
                                           quantile_lines = TRUE,
                                           quantiles = c(0.025, 0.05, 0.95, 0.975),
                                           show.legend = FALSE) +
    ## scale_fill_gradient(low = "white", high = "orange") +
    scale_fill_manual(values = c('#ffa600', '#7a5195', '#003f5c', '#7a5195', '#ffa600')) +
    scale_y_discrete('Intercept') + 
    scale_x_log10('Area (cm²)') +
    facet_grid(Zone~Year) + 
    theme_classic()
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Partial.type_',model,'.png'),
       g, 
       width = 6, height = 6,
       dpi = 100
       )
## ----end
## ---- CoralGrowth brms Model 6.1 Partial plots derived contrasts
## Margninal means of Site Type and Zone
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
data.sub1 %>%
    pull(Treatment) %>%
    levels()
## - BA mean
## - FR mean
## - FL mean
## - S mean
## - D mean
cmat <- list('Site Type' = cbind('BA' = c(0.5, 0.5, 0, 0, 0),
                                 'FR' = c(0, 0, 0.5, 0, 0.5),
                                 'FL' = c(0,0,0,1,0)),
             'Zone' = cbind('S' = c(1/3, 0, 0, 1/3, 1/3),
                            'D' = c(0, 1/2, 1/2, 0, 0))
             )
mod %>% emmeans(~Treatment) %>%
    contrast(method = list(Treatment = cmat)) %>%
    knitr::kable()
## ----end
## ---- CoralGrowth brms Model 6.1 Partial plots derived contrasts 1
mod.sum <-
    mod %>% emmeans(~Treatment) %>%
    contrast(method = list(Treatment = cmat)) %>%
    tidybayes::tidy_draws() %>%
    exp() %>%
    `*`(10000) %>%
    posterior::summarise_draws(median,
                               HDInterval::hdi
                               ) %>%
    mutate(variable = str_replace_all(variable, "contrast Treatment.(Site Type|Zone).", ""))
mod.sum
## ----end
## ---- CoralGrowth brms Model 6.1 Partial plots derived contrasts figure
g <- mod.sum %>%
      ggplot(aes(y = median, x = variable)) +
      geom_pointrange(aes(ymin = lower, ymax = upper)) +
      scale_y_log10("Area") +
      theme_classic() 
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Partial_contrast.Site_type1_',model,'.png'),
       g, 
       width = 6, height = 4,
       dpi = 100
       )

## ----end
## ---- CoralGrowth brms Model 6.1 Partial plots derived contrasts effects
## Margninal effects of Site Type and Zone
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
data.sub1 %>%
    pull(Treatment) %>%
    levels()
cmat <- list('BA vs FR (S)' = c(1, 0, 0, 0, -1),
             'BA vs FR (D)' = c(0, 1, -1, 0, 0),
             'BA vs FR' = c(0.5, 0.5, -0.5, 0, -0.5),
             'BA vs FL (S)' = c(1, 0, 0, -1, 0),
             'FR vs FL (S)' = c(0, 0, 0, 1, -1),
             'FR vs FL' = c(0, 0, 0.5, -1, 0.5),
             'BA vs FL' = c(0.5, 0.5, 0, -1, 0))
mod %>% emmeans(~Treatment) %>%
    contrast(method = list(Treatment = cmat)) %>%
    knitr::kable()
## ----end
## ---- CoralGrowth brms Model 6.1 Partial plots derived contrasts effects1
mod.sum <-
    mod %>% emmeans(~Treatment) %>%
    contrast(method = list(Treatment = cmat)) %>%
    tidybayes::tidy_draws() %>%
    exp() %>%
    posterior::summarise_draws(median,
                               HDInterval::hdi,
                               `P(>1)` = ~ sum(.x>1)/length(.x),
                               `P(>10%)` = ~ sum(.x>1.1)/length(.x),
                               `P(<1)` = ~ sum(.x<1)/length(.x),
                               `P(<10%)` = ~ sum(.x<0.909)/length(.x)
                               ) %>%
    mutate(variable = str_replace_all(variable, "contrast Treatment.", "")) 
mod.sum %>% knitr::kable()
## ----end
## ---- CoralGrowth brms Model 6.1 Partial plots derived contrasts effects1 figure
g <- mod.sum %>%
      ggplot(aes(x = median, y = variable)) +
      geom_vline(xintercept = 1, linetype = 'dashed') +
      geom_pointrange(aes(xmin = lower, xmax = upper)) +
      scale_x_continuous("Effect (ratio)", trans = scales::log2_trans()) +
      theme_classic() 
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Partial_contrast.Site_type2_',model,'.png'),
       g, 
       width = 6, height = 4,
       dpi = 100
       )

## ----end

