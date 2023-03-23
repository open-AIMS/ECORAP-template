## ---- loadFunctions
source('../R/functions.R')
## ----end

## Fish Abundance ***************************************************************

## ---- FishAbundance readData
data <- read_csv(paste0(DATA_PATH,
                        "primary/small_fish_complexity_extract.csv"),
                 trim_ws=TRUE)
## ----end

## ---- FishAbundance glimpse
data %>% glimpse() 
## ----end

## ---- FishAbundance process
## The number of plot/transects
data %>% dplyr::select(plot_id, transect) %>% distinct() %>% dim()

## The expected number of plot/transects
data.full <- data %>%
    tidyr::expand(plot_id, transect)
data.full %>% dim()

data <- 
    data %>% 
    mutate(Pomacentrus = str_detect(genus_species, "Pomacentrus")) %>%
    ## aggregate to plot/transect level
    group_by(plot_id, transect) %>%
    summarise(Abundance = sum(abund, na.rm = TRUE),
              Abundance.spp1 = if_else(first(Pomacentrus),
                                       sum(abund),
                                       0)) %>%
    ungroup() %>%
    ## ensure that missing values are included
    ## full_join(data.full) %>%
    mutate(Abundance = ifelse(is.na(Abundance), 0, Abundance),
           Abundance.spp1 = ifelse(is.na(Abundance.spp1), 0, Abundance.spp1)) %>%
    ## extract Reef, Site and Zone info from the Site code
    mutate(
        Region = factor(str_replace(plot_id, "(..).*", "\\1")),
        Reef = factor(str_replace(plot_id, "([^_]*)_.*", "\\1")),
        Site = factor(str_replace(plot_id, "[^_]*_(...).*", "\\1")),
        Site_Type = factor(str_replace(Site, "(.*)[0-9]$", "\\1")),
        Zone = factor(str_replace(plot_id, "[^_]*_...(.).*", "\\1")),
        Plot = factor(str_replace(plot_id, "[^_]*_.*_(.*)", "\\1")),
        ) %>%
    mutate(Site = factor(paste(Reef, Site, sep = "_")),
           Site_Zone = factor(paste(Site, Zone, sep = "_")),
           Plot = factor(paste(Site, Zone, Plot, sep = "_")),
           Transect = factor(paste(Plot, transect, sep = "_"))
           ) %>%
    dplyr::select(-plot_id, -transect) 

data %>% glimpse()
saveRDS(data, file = paste0(DATA_PATH, "processed/data_fish.RData"))
## ----end

## EDA ==========================================================================
## ---- FishAbundance EDA
g <- data %>%
    ggplot(aes(y = Abundance, x = Site_Type, color = Zone)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.02, dodge.width = 0.9),
               alpha = 0.1) +
    geom_violin(fill = NA) +
    facet_wrap(~Region, scales='free_y') +
    theme_bw() +
    scale_y_continuous(trans = scales::pseudo_log_trans())
ggsave(filename = paste0(FIGS_PATH, "/EDA1_fish.png"),
       g,
       width = 10,
       height = 5,
       dpi = 100)

g <- data %>%
    ggplot(aes(y = Abundance.spp1, x = Site_Type, color = Zone)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.02, dodge.width = 0.9),
               alpha = 0.1) +
    geom_violin(fill = NA) +
    facet_wrap(~Region, scales='free_y') +
    theme_bw() +
    scale_y_continuous(trans = scales::pseudo_log_trans())
ggsave(filename = paste0(FIGS_PATH, "/EDA2_fish.png"),
       g,
       width = 10,
       height = 5,
       dpi = 100)
## ----end

## ---- FishAbundance processing further
data.sub <- data %>% 
    droplevels() %>% 
    mutate(
        Treatment = interaction(Site_Type, Zone, sep = '_')
    ) %>% 
    droplevels() 

## Note that the first level has only a single observation - not a
## good choice as the intercept
data.sub %>% pull(Treatment) %>% levels()

data.sub <- data.sub %>%
  mutate(Treatment = forcats::fct_relevel(Treatment, 'BA_S')) 
saveRDS(data.sub, file = paste0(DATA_PATH, "processed/data.sub_fish.RData"))
## ----end

## Modelling ====================================================================

## ---- FishAbundance brms Model 1.1 priors
data.sub <- readRDS(file = paste0(DATA_PATH, "processed/data.sub_fish.RData"))
model <- 'mod1.1_fish'
# SD between regions
data.sub %>%
    group_by(Region) %>%
    summarise(across(Abundance,  list(Mean = ~ log(mean(.x)),
                                         Median = ~ log(median(.x)),
                                         SD = ~ log(sd(.x)),
                                         MAD = ~ log(mad(.x)),
                                         N = ~ length(.x)))) %>%
    ungroup() %>% 
    summarise(SD = sd(Abundance_Median, na.rm = TRUE))

# SD between reefs
data.sub %>%
    group_by(Region, Reef) %>%
    summarise(across(Abundance,  list(Mean = ~ log(mean(.x)),
                                         Median = ~ log(median(.x)),
                                         SD = ~ log(sd(.x)),
                                         MAD = ~ log(mad(.x)),
                                         N = ~ length(.x)))) %>%
    ungroup() %>% 
    summarise(SD = sd(Abundance_Median, na.rm = TRUE))

# SD between sites
data.sub %>%
    group_by(Region, Reef, Site) %>%
    summarise(across(Abundance,  list(Mean = ~ log(mean(.x)),
                                         Median = ~ log(median(.x)),
                                         SD = ~ log(sd(.x)),
                                         MAD = ~ log(mad(.x)),
                                         N = ~ length(.x)))) %>%
    ungroup() %>% 
    summarise(SD = sd(Abundance_Median, na.rm = TRUE))

# SD between plots
data.sub %>%
    group_by(Region, Reef, Site, Plot) %>%
    summarise(across(Abundance,  list(Mean = ~ log(mean(.x)),
                                         Median = ~ log(median(.x)),
                                         SD = ~ log(sd(.x)),
                                         MAD = ~ log(mad(.x)),
                                         N = ~ length(.x)))) %>%
    ungroup() %>% 
    summarise(SD = sd(Abundance_Median, na.rm = TRUE))

data.sub %>%
    group_by(Region, Reef, Site, Plot) %>%
    summarise(across(Abundance,  list(Mean = ~ log(mean(.x)),
                                         Median = ~ log(median(.x)),
                                         SD = ~ log(sd(.x)),
                                         MAD = ~ log(mad(.x)),
                                         N = ~ length(.x)))) %>%
    ungroup() %>% 
    summarise(Median = median(Abundance_Median, na.rm = TRUE),
              SD = sd(Abundance_Median, na.rm = TRUE))

form <- bf(Abundance ~ 1 + 
               (1 | Reef/Site/Zone/Plot),
           family = poisson())
get_prior(form, data = data.sub)
## ----end
## ---- FishAbundance brms Model 1.1 fit
priors <- set_prior("student_t(3, 4.3, 1)", class = "Intercept") +
    set_prior("student_t(3, 0, 1)", class = "sd") 

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
## ---- FishAbundance brms Model 1.1 posterior_prior
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
g <- mod %>% SUYR_prior_and_posterior() +
    theme_classic() +
    theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Posterior_prior_',model,'.png'),
       g,
       width = 10,
       height = 8,
       dpi = 100
       )
## ----end
## ---- FishAbundance brms Model 1.1 MCMC diagnostics
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
## ---- FishAbundance brms Model 1.1 ppc
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc1_',model,'.png'),
       mod %>% pp_check(type = "dens_overlay", ndraws = 100) + scale_x_log10(),
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
## ---- FishAbundance brms Model 1.1 DHARMa
brms_dharma_res <- make_brms_dharma_res(
    mod, integerResponse = TRUE
) 
ggsave(filename = paste0(OUTPUT_PATH, 'figures/DHARMa1_',model,'.png'),
       wrap_elements(~testUniformity(brms_dharma_res)) +
       wrap_elements(~plotResiduals(brms_dharma_res, form = factor(rep(1, nrow(data.sub))))) +
       wrap_elements(~testDispersion(brms_dharma_res)),
       width = 12,
       height = 4,
       dpi = 100)
## ----end

## ---- FishAbundance brms Model 2.1 priors
data.sub <- readRDS(file = paste0(DATA_PATH, "processed/data.sub_fish.RData"))
model <- 'mod2.1_fish'

data.sub %>%
    group_by(Site_Type, Zone) %>%
    summarise(across(Abundance,  list(Mean = ~ log(mean(.x)),
                                         Median = ~ log(median(.x)),
                                         SD = ~ log(sd(.x)),
                                         MAD = ~ log(mad(.x)),
                                         N = ~ length(.x)))) 

form <- bf(Abundance ~ Site_Type*Zone + 
               (1 | Reef/Site/Plot),
           family = poisson())
get_prior(form, data = data.sub)
## ----end
## ---- FishAbundance brms Model 2.1 fit
priors <- set_prior("student_t(3, 4.3, 1)", class = "Intercept") +
    set_prior("student_t(3, 0, 1)", class = "b") +
    set_prior("student_t(3, 0, 1)", class = "sd") 

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
## ---- FishAbundance brms Model 2.1 posterior_prior
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
g <- mod %>% SUYR_prior_and_posterior() +
    theme_classic() +
    theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Posterior_prior_',model,'.png'),
       g,
       width = 10,
       height = 8,
       dpi = 100
       )
## ----end
## ---- FishAbundance brms Model 2.1 MCMC diagnostics
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
## ---- FishAbundance brms Model 2.1 ppc
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc1_',model,'.png'),
       mod %>% pp_check(type = "dens_overlay", ndraws = 100) + scale_x_log10(),
       width = 5,
       height = 4,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc2_',model,'.png'),
       mod %>% pp_check(group = "Site_Type", type = "intervals_grouped"),
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
## ---- FishAbundance brms Model 2.1 DHARMa
brms_dharma_res <- make_brms_dharma_res(
    mod, integerResponse = TRUE
) 
ggsave(filename = paste0(OUTPUT_PATH, 'figures/DHARMa1_',model,'.png'),
       wrap_elements(~testUniformity(brms_dharma_res)) +
       wrap_elements(~plotResiduals(brms_dharma_res, form = factor(rep(1, nrow(data.sub))))) +
       wrap_elements(~plotResiduals(brms_dharma_res)) +
       wrap_elements(~testDispersion(brms_dharma_res)),
       width = 12,
       height = 8,
       dpi = 100)
## ----end
## ---- FishAbundance brms Model 2.1 Partial plots derived means
em.sum <- mod %>% emmeans(~Site_Type|Zone) %>%
    tidybayes::tidy_draws() %>%
    exp() %>%
    posterior::summarise_draws(median,
                               HDInterval::hdi,
                               rhat,
                               ess_bulk,
                               ess_tail,
                               length) %>%
    mutate(variable = str_replace_all(variable, "Site_Type ", "")) %>%
    mutate(variable = str_replace_all(variable, ", Zone ", "_")) %>%
    separate(variable, into = c('Site_Type', 'Zone'), sep = '_')
## ----end
## ---- FishAbundance brms Model 2.1 Partial plots derived means figure
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
g <- mod %>% emmeans(~Site_Type | Zone) %>%
    tidybayes::gather_emmeans_draws() %>%
    mutate(.value = exp(.value)) %>%
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
    facet_grid(~Zone) + 
    theme_classic()
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Partial.type_',model,'.png'),
       g, 
       width = 6, height = 6,
       dpi = 100
       )
## ----end
## ---- FishAbundance brms Model 2.1 Partial plots derived contrasts
## Margninal means of Site Type and Zone
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
mod %>% emmeans(~Site_Type | Zone, type = "response") %>%
    pairs() %>%
    knitr::kable()

mod.sum <- 
    mod %>% emmeans(~Site_Type | Zone) %>%
    pairs() %>%
    tidybayes::tidy_draws() %>%
    exp() %>%
    posterior::summarise_draws(median,
                               HDInterval::hdi,
                               `P(>1)` = ~ mean(.x>1),
                               `P(>10%)` = ~ mean(.x>1.1),
                               `P(<1)` = ~ mean(.x<1),
                               `P(<10%)` = ~ mean(.x<0.909),
                               length) %>%
    mutate(variable = str_replace_all(variable, "contrast ", "")) %>% 
    mutate(variable = str_replace_all(variable, " Zone ", "")) %>% 
    separate(variable, into = c('variable', 'Zone'), sep = ",") 
mod.sum %>% knitr::kable()
## ----end
## ---- FishAbundance brms Model 2.1 Partial plots derived contrasts figure
g <- mod.sum %>%
    ggplot(aes(x = median, y = variable)) +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    geom_pointrange(aes(xmin = lower, xmax = upper)) +
    scale_x_continuous("Effect (ratio)", trans = scales::log2_trans()) +
    facet_grid(~Zone) + 
    theme_classic() 
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Partial_contrast.Site_type1_',model,'.png'),
       g, 
       width = 6, height = 4,
       dpi = 100
       )

## ----end

## ---- FishAbundance brms Model 3.1 priors
data.sub <- readRDS(file = paste0(DATA_PATH, "processed/data.sub_fish.RData"))
model <- 'mod3.1_fish'

data.sub %>%
    group_by(Site_Type, Zone) %>%
    summarise(across(Abundance.spp1,  list(Mean = ~ log(mean(.x)),
                                         Median = ~ log(median(.x)),
                                         SD = ~ log(sd(.x)),
                                         MAD = ~ log(mad(.x)),
                                         N = ~ length(.x)))) 

form <- bf(Abundance.spp1 ~ Site_Type*Zone + 
               (1 | Reef/Site/Plot),
           family = poisson())
get_prior(form, data = data.sub)
## ----end
## ---- FishAbundance brms Model 3.1 fit
priors <- set_prior("student_t(3, 4, 1)", class = "Intercept") +
    set_prior("student_t(3, 0, 1)", class = "b") +
    set_prior("student_t(3, 0, 1)", class = "sd") 

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
## ---- FishAbundance brms Model 3.1 posterior_prior
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
g <- mod %>% SUYR_prior_and_posterior() +
    theme_classic() +
    theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Posterior_prior_',model,'.png'),
       g,
       width = 10,
       height = 8,
       dpi = 100
       )
## ----end
## ---- FishAbundance brms Model 3.1 MCMC diagnostics
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
## ---- FishAbundance brms Model 3.1 ppc
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc1_',model,'.png'),
       mod %>% pp_check(type = "dens_overlay", ndraws = 100) + scale_x_log10(),
       width = 5,
       height = 4,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc2_',model,'.png'),
       mod %>% pp_check(group = "Site_Type", type = "intervals_grouped"),
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
## ---- FishAbundance brms Model 3.1 DHARMa
brms_dharma_res <- make_brms_dharma_res(
    mod, integerResponse = TRUE
) 
ggsave(filename = paste0(OUTPUT_PATH, 'figures/DHARMa1_',model,'.png'),
       wrap_elements(~testUniformity(brms_dharma_res)) +
       wrap_elements(~plotResiduals(brms_dharma_res, form = factor(rep(1, nrow(data.sub))))) +
       wrap_elements(~plotResiduals(brms_dharma_res)) +
       wrap_elements(~testDispersion(brms_dharma_res)) +
       wrap_elements(~testZeroInflation(brms_dharma_res)),
       width = 12,
       height = 8,
       dpi = 100)
## ----end

## ---- FishAbundance brms Model 4.1 priors
data.sub <- readRDS(file = paste0(DATA_PATH, "processed/data.sub_fish.RData"))
model <- 'mod4.1_fish'

data.sub %>%
    group_by(Site_Type, Zone) %>%
    summarise(across(Abundance.spp1,  list(Mean = ~ log(mean(.x)),
                                         Median = ~ log(median(.x)),
                                         SD = ~ log(sd(.x)),
                                         MAD = ~ log(mad(.x)),
                                         N = ~ length(.x)))) 

form <- bf(Abundance.spp1 ~ Site_Type*Zone + 
               (1 | Reef/Site/Plot),
           zi = ~ 1,
           family = zero_inflated_poisson())
get_prior(form, data = data.sub)
## ----end
## ---- FishAbundance brms Model 4.1 fit
priors <- set_prior("student_t(3, 4, 1)", class = "Intercept") +
    set_prior("student_t(3, 0, 1)", class = "b") +
    set_prior("student_t(3, 0, 1)", class = "sd") +
    set_prior("logistic(0, 1)", class = "Intercept", dpar = "zi") 
    
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
## ---- FishAbundance brms Model 4.1 posterior_prior
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
g <- mod %>% SUYR_prior_and_posterior() +
    theme_classic() +
    theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Posterior_prior_',model,'.png'),
       g,
       width = 10,
       height = 8,
       dpi = 100
       )
## ----end
## ---- FishAbundance brms Model 4.1 MCMC diagnostics
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
## ---- FishAbundance brms Model 4.1 ppc
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc1_',model,'.png'),
       mod %>% pp_check(type = "dens_overlay", ndraws = 100) + scale_x_log10(),
       width = 5,
       height = 4,
       dpi = 100
       )
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Ppc2_',model,'.png'),
       mod %>% pp_check(group = "Site_Type", type = "intervals_grouped"),
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
## ---- FishAbundance brms Model 4.1 DHARMa
brms_dharma_res <- make_brms_dharma_res(
    mod, integerResponse = TRUE
) 
ggsave(filename = paste0(OUTPUT_PATH, 'figures/DHARMa1_',model,'.png'),
       wrap_elements(~testUniformity(brms_dharma_res)) +
       wrap_elements(~plotResiduals(brms_dharma_res, form = factor(rep(1, nrow(data.sub))))) +
       wrap_elements(~plotResiduals(brms_dharma_res)) +
       wrap_elements(~testDispersion(brms_dharma_res)) +
       wrap_elements(~testZeroInflation(brms_dharma_res)),
       width = 12,
       height = 8,
       dpi = 100)
## ----end

## ---- FishAbundance brms Model 4.1 Partial plots derived means
em.sum <- mod %>% emmeans(~Site_Type|Zone) %>%
    tidybayes::tidy_draws() %>%
    exp() %>%
    posterior::summarise_draws(median,
                               HDInterval::hdi,
                               rhat,
                               ess_bulk,
                               ess_tail,
                               length) %>%
    mutate(variable = str_replace_all(variable, "Site_Type ", "")) %>%
    mutate(variable = str_replace_all(variable, ", Zone ", "_")) %>%
    separate(variable, into = c('Site_Type', 'Zone'), sep = '_')
## ----end
## ---- FishAbundance brms Model 4.1 Partial plots derived means figure
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
g <- mod %>% emmeans(~Site_Type | Zone) %>%
    tidybayes::gather_emmeans_draws() %>%
    mutate(.value = exp(.value)) %>%
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
    facet_grid(~Zone) + 
    theme_classic()
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Partial.type_',model,'.png'),
       g, 
       width = 6, height = 6,
       dpi = 100
       )
## ----end
## ---- FishAbundance brms Model 4.1 Partial plots derived contrasts
## Margninal means of Site Type and Zone
load(file = paste0(DATA_PATH, "modelled/",model,".RData"))
mod %>% emmeans(~Site_Type | Zone, type = "response") %>%
    pairs() %>%
    knitr::kable()

mod.sum <- 
    mod %>% emmeans(~Site_Type | Zone) %>%
    pairs() %>%
    tidybayes::tidy_draws() %>%
    exp() %>%
    posterior::summarise_draws(median,
                               HDInterval::hdi,
                               `P(>1)` = ~ mean(.x>1),
                               `P(>10%)` = ~ mean(.x>1.1),
                               `P(<1)` = ~ mean(.x<1),
                               `P(<10%)` = ~ mean(.x<0.909),
                               length) %>%
    mutate(variable = str_replace_all(variable, "contrast ", "")) %>% 
    mutate(variable = str_replace_all(variable, " Zone ", "")) %>% 
    separate(variable, into = c('variable', 'Zone'), sep = ",") 
mod.sum %>% knitr::kable()
## ----end
## ---- FishAbundance brms Model 4.1 Partial plots derived contrasts figure
g <- mod.sum %>%
    ggplot(aes(x = median, y = variable)) +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    geom_pointrange(aes(xmin = lower, xmax = upper)) +
    scale_x_continuous("Effect (ratio)", trans = scales::log2_trans()) +
    facet_grid(~Zone) + 
    theme_classic() 
ggsave(filename = paste0(OUTPUT_PATH, 'figures/Partial_contrast.Site_type1_',model,'.png'),
       g, 
       width = 6, height = 4,
       dpi = 100
       )

## ----end
