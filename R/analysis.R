## ---- loadFunctions
source('../R/functions.R')
## ----end

## ---- CoralGrowth readData
data <- read_csv(paste0(DATA_PATH,
                        "primary/All_size_2021-2022_ML.csv"),
                 trim_ws=TRUE)
## ----end

## ---- CoralGrowth glimpse
data %>% glimpse() 
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
        Reef_Site_Zone_Plot = factor(paste(Reef_Site_Zone, Plot, sep = '_'))) %>%
    ## Define Site_code as each unique Reef/Site/Zone combination
    mutate(Site_code = factor(paste(Reef, Site, Zone, sep = '_'))) %>%
    mutate(Site_type = factor(str_sub(Site, 1, 2)))

data %>% glimpse()
saveRDS(data, file = paste0(DATA_PATH, "processed/data.RData"))
## ----end
