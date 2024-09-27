# Juvenile mangrove growth is restricted by belowground competition with mature tree roots
# Jack W Hill, Vicki Bennion, Emer T Cunningham, John M Dwyer, Catherine E Lovelock

##### 0. set up workspace #####

library(betareg)
library(broom)
library(DHARMa)
library(emmeans)
library(glmmTMB)
library(patchwork)
library(ragg)
library(readxl)
library(survival)
library(survminer)

# last to avoid conflicts
library(brms) 
library(tidyverse)

set.seed(070922) # for position_jitter()

bag_wgt_coles_snack <- 1.745 # sample bag weight

# calculate basal area from raw circumference data
calc_ba <- function (dat) {
  max_circs <- max(str_count(dat$circ, ","), 
                   na.rm = TRUE) + 1
  
  circ_cols <- paste0("circ_", seq(1:max_circs))
  
  ba_dat <- dat %>% 
    separate(col = circ,
             into = circ_cols,
             sep = ", ",
             fill = "right") %>% 
    pivot_longer(cols = starts_with("circ")) %>% 
    group_by(tree_id) %>% 
    mutate(value = as.numeric(value),
           value = pi * ((value / pi / 2)^2)) %>% 
    summarise(ba = sum(value, na.rm = TRUE)) %>% 
    ungroup()
  
  return(ba_dat)
}

# custom ggplot theme
theme_root_cut <- function() {
  theme_classic(base_size = 12, base_family = "sans") %+replace%
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(size = rel(0.75)),
      axis.text.x = element_text(margin = margin(2, 0, 3, 0)),
      axis.text.y = element_text(margin = margin(0, 2, 0, 1)),
      axis.ticks = element_line(colour = "black"),
      axis.line = element_line(colour = "black"),
      legend.position = "none",
      plot.tag = element_text(face = "bold")
    )
}

# feature colours for plots
triL_col <- "#7D31ED" # purple
triM_col <- "#ED7D31" # orange
triR_col <- "#13d863" # green

#
##### 1. import and tidy data #####

raw_height <- read_excel("data/raw/root-cutting-tinchi.xlsx",
                         sheet = "heights")

raw_juves <- read_excel("data/raw/root-cutting-tinchi.xlsx",
                        sheet = "juveniles") %>% 
  mutate(treatment = case_when(is.na(treatment) == TRUE ~ "control",
                               .default = treatment))

juves <- raw_juves %>% 
  select(juve_id, treatment)

height <- raw_height %>% 
  mutate(height = case_when(!is.na(ruler_sink) ~ 
                              height - ruler_sink,
                            .default = height)) %>% 
  left_join(., juves, 
            by = "juve_id")

dat <- height %>% 
  group_by(juve_id) %>% 
  arrange(week,
          .by_group = TRUE) %>% 
  mutate(pw_change = case_when(week == 0 ~ NA_real_,
                               week > 0 ~ (height - lag(height, 1)) / 
                                 (week - lag(week, 1)),
                               .default = 10000), # error value
         pc_pw_change = pw_change/first(height) * 100,
         dead = case_when(is.na(height) ~ TRUE,
                          juve_id != 40 & str_detect(notes, "dead") ~ TRUE,
                          .default = FALSE)) %>%
  ungroup()

# remove any measurements following the first dead week
dead_weeks <- dat %>%
  group_by(juve_id) %>% 
  filter(dead == TRUE, 
         .preserve = TRUE) %>% # juves that survived the whole period get `Inf`
  summarise(first_dead_week = min(week))

reduc_dat <- dat %>% 
  left_join(., dead_weeks,
            by = "juve_id") %>% 
  filter(week <= first_dead_week)

# pcq = point-centred quarter
pcq_dat_raw <- read_excel("data/raw/root-cutting-tinchi.xlsx",
                          sheet = "pcq")

ba_dat <- pcq_dat_raw %>% 
  mutate(tree_id = paste0(juve_id, "-q", quarter)) %>% 
  calc_ba() %>% 
  filter(ba > 0) %>% 
  mutate(juve_id = as.numeric(str_split_i(tree_id, "-q", 1)),
         ba = ba) %>% # cm2
  group_by(juve_id) %>% 
  summarise(ba = mean(ba)) %>% 
  ungroup()

pcq_only_dat <- pcq_dat_raw %>% 
  filter(!is.na(dist_from)) %>% # because some quarters have missingness
  left_join(., ba_dat,
            by = "juve_id") %>% 
  group_by(juve_id) %>% 
  summarise(density = mean(dist_from),
            ba = mean(ba),
            tree_height = mean(height))
# not enough points to do traditional PCQ analysis (i.e. scaling up to 
# per hectare values)

pcq_dat <- reduc_dat %>% 
  left_join(., pcq_only_dat,
            by = "juve_id")

#
##### 2. analysis -- juvenile growth #####

mod_juve_hgt <- glmmTMB(height ~ week + treatment + week:treatment + (week|juve_id), 
                        data = reduc_dat %>% filter(!is.na(height)))
summary(mod_juve_hgt)
# result: week and treatment interaction, and main effects

emmeans(mod_juve_hgt, pairwise ~ treatment|week, 
        at = list("week" = 52)) 
# result: heights at week 52

emtrends(mod_juve_hgt, pairwise ~ treatment, "week")
# result: weekly growth rates

emmeans_output_full_year <- emmeans(mod_juve_hgt, pairwise ~ treatment|week, 
                                    at = list("week" = -1:53)) 
emmeans_full_year <- tidy(emmeans_output_full_year$emmeans)

#
##### 3. plotting -- juvenile growth #####

plot_reduc_dat <- reduc_dat %>% 
  group_by(week, treatment) %>% 
  filter(!is.na(height)) %>% 
  summarise(mean_hgt = mean(height),
            se_hgt = sd(height)/sqrt(n())) %>% 
  ungroup() %>% 
  # offset cut points slightly on x-axis for aesthetics
  mutate(week = case_when(treatment == "cut" ~ week + 0.7,
                          TRUE ~ week))

plot_growth <- ggplot(data = emmeans_full_year,
                       aes(x = week)) +
  geom_ribbon(aes(ymin = estimate - std.error,
                  ymax = estimate + std.error,
                  group = treatment),
              fill = "grey90",
              colour = "white") +
  geom_errorbar(data = plot_reduc_dat,
                aes(ymin = mean_hgt - se_hgt,
                    ymax = mean_hgt + se_hgt,
                    colour = treatment),
                width = 0) +
  geom_point(data = plot_reduc_dat,
             aes(y = mean_hgt,
                 colour = treatment)) +
  geom_line(aes(y = estimate,
                colour = treatment),
            size = 0.75) +
  scale_colour_manual(values = c(triR_col,
                                 triM_col)) +
  scale_y_continuous(name = "Juvenile height (cm)",
                     expand = expansion()) +
  scale_x_continuous(name = "Sampling week",
                     expand = expansion())+
  coord_cartesian(xlim = c(-1, 52.5)) +
  theme_root_cut()
plot_growth

agg_png("fig-output/figure-1.png",
        width = 8, height = 6, units = "cm",
        res = 1080, scaling = 0.7)
plot_growth
dev.off()

#
##### 4. analysis -- juvenile survival #####

surv_dat <- reduc_dat %>% 
  group_by(juve_id) %>% 
  distinct(first_dead_week, .keep_all = TRUE) %>% 
  ungroup() %>% 
  mutate(status = case_when(is.infinite(first_dead_week) ~ 0,
                            # 0 indicates right-censoring (juvenile did not die)
                            # 1 indicates that juvenile died at some point
                            .default = 1), 
         first_dead_week = case_when(is.infinite(first_dead_week) ~ 51,
                                     # week 51 indicates the last 'follow-up visit' 
                                     # for juves which did not die
                                     .default = first_dead_week)) %>% 
  select(juve_id, treatment, first_dead_week, status)

model_cox <- coxph(Surv(time = first_dead_week, event = status) ~ treatment,
                   data = surv_dat)
summary(model_cox)
# result: survival rates

#
##### 5. plotting -- juvenile survival #####

longevity_dat <- dat %>% 
  left_join(., dead_weeks,
            by = "juve_id") %>% 
  group_by(juve_id) %>% 
  mutate(dead = case_when(week < first_dead_week ~ FALSE,
                          .default = TRUE)) %>% 
  group_by(week, treatment) %>% # replaces previous grouping
  count(dead) %>% 
  ungroup() %>% 
  filter(dead == FALSE) %>% 
  mutate(n = n / 20) # sample size per treatment

plot_longevity <- ggplot(data = longevity_dat,
                         aes(x = week,
                             y = n,
                             colour = treatment)) +
  geom_line(size = 0.75) +
  geom_point(size = 2) +
  annotate("text",
           label = "control juveniles",
           x = 18, y = 0.88,
           hjust = 0, vjust = 0,
           colour = triR_col,
           size = 3,
           fontface = "bold") +
  annotate("text",
           label = "cut juveniles",
           x = 30, y = 0.54,
           hjust = 0, vjust = 0,
           colour = triM_col,
           size = 3,
           fontface = "bold") +
  scale_colour_manual(values = c(triR_col,
                                 triM_col)) +
  scale_y_continuous(name = "Surviving juveniles (proportion)",
                     limits = c(0, 1),
                     expand = expansion(add = c(0, 0.05))) +
  scale_x_continuous(name = "Sampling week",
                     expand = expansion(add = c(0, 0.53))) +
  theme_root_cut()
plot_longevity

agg_png("fig-output/figure-2.png",
        width = 8, height = 6, units = "cm",
        res = 1080, scaling = 0.7)
plot_longevity
dev.off()

#
##### 6. analysis -- pneumatophore quadrats #####

pneumat_dat <- raw_juves %>% 
  mutate(pneumats_total = pneumats_dead + pneumats_alive,
         pneumats_dead_prop = pneumats_dead / pneumats_total,
         pneumats_alive_prop = 1 - pneumats_dead_prop)

# pneumatophore density
pneumat_dens_control <- pneumat_dat %>% 
  filter(treatment == "control") %>% 
  pull(pneumats_total)

pneumat_dens_cut <- pneumat_dat %>% 
  filter(treatment == "cut") %>% 
  pull(pneumats_total)

t.test(pneumat_dens_control, pneumat_dens_cut)

# zero-one-inflated beta regression required for pneumatophore mortality data
# because there are a number of 1s and 0s
zoib_pneumats <- brm(
  bf(pneumats_alive_prop ~ treatment, 
     phi ~ treatment, 
     zoi ~ treatment,  
     coi ~ treatment),
  data = pneumat_dat %>% 
    filter(!is.na(pneumats_alive_prop)),
  family = zero_one_inflated_beta(),
  file = "data/output/zoib_pneumat_model"
)

pneumat_emms <- emmeans(zoib_pneumats, ~ treatment, 
                        type = "response", epred = TRUE) 
pneumat_emms_summary_dat <- as_tibble(pneumat_emms)

contrast(emmeans(zoib_pneumats, ~ treatment, 
                 type = "response", epred = TRUE),
         method = "pairwise", point.est = mean)
# result: difference in proportion of live pneumatophores
        
#
##### 7. plotting -- pneumatophore quadrats #####

pneumat_summary_dat <- pneumat_dat %>% 
  filter(!is.na(pneumats_total)) %>% 
  group_by(treatment) %>% 
  summarise(mean_total = mean(pneumats_total),
            se_total = sd(pneumats_total) / sqrt(n()),
            mean_alive = mean(pneumats_alive_prop),
            se_alive = sd(pneumats_alive_prop) / sqrt(n()))

pneumat_density_xlab <- expression("Pneumatophore density (# 90 cm"^-2*")")

plot_pneumat_density <- ggplot(mapping = aes(x = treatment)) +
  geom_point(data = pneumat_dat,
             aes(y = pneumats_total),
             colour = "grey80",
             position = position_jitter(width = 0.05)) +
  geom_pointrange(data = pneumat_summary_dat,
                  aes(y = mean_total,
                      ymin = mean_total - se_total,
                      ymax = mean_total + se_total),
                  size = 0.3) +
  scale_y_continuous(name = pneumat_density_xlab,
                     limits = c(0, NA),
                     expand = expansion(add = c(0, 3))) +
  scale_x_discrete(name = "Root treatment",
                   labels = c("Control", "Cut")) +
  theme_root_cut()
plot_pneumat_density

plot_pneumat_mortality <- ggplot(mapping = aes(x = treatment)) +
  geom_point(data = pneumat_dat,
             aes(y = pneumats_alive_prop),
             colour = "grey80",
             position = position_jitter(width = 0.05)) +
  geom_pointrange(data = pneumat_emms_summary_dat,
                  aes(y = emmean,
                      ymin = lower.HPD,
                      ymax = upper.HPD),
                  size = 0.3) +
  scale_y_continuous(name = "Live pneumatophores (proportion)",
                     expand = expansion(),
                     limits = c(0, 1.02)) +
  scale_x_discrete(name = "Root treatment",
                   labels = c("Control", "Cut")) +
  theme_root_cut()
plot_pneumat_mortality

agg_png("fig-output/figure-3.png",
        width = 15, height = 8, units = "cm",
        res = 1080, scaling = 0.9)
plot_pneumat_density + plot_pneumat_mortality
dev.off()

#
##### 8. analysis -- covariates #####

# initial juvenile height
start_hgt_control <- dat %>% 
  filter(week == 0, treatment == "control") %>% 
  pull(height)

start_hgt_cut <- dat %>% 
  filter(week == 0,
         treatment == "cut") %>% 
  pull(height)

t.test(start_hgt_control, start_hgt_cut)

# canopy cover (densiometer dots)
densi_dat <- raw_juves %>% 
  mutate(across(starts_with("densiometer"),
                ~ 100 - . * 1.04), # gives canopy cover (%)
         # average of three sampled cardinal directions, 
         # and convert % -> proportion
         canopy_cover = (densiometer_e + densiometer_n + densiometer_w) / 3 / 100) %>% 
  filter(!is.na(canopy_cover))

betamod_canopy <- betareg(canopy_cover ~ treatment, 
                          data = densi_dat)
summary(betamod_canopy)

# PCQ variables
pcq_tests_dat <- juves %>% 
  left_join(., pcq_only_dat,
            by = "juve_id") %>% 
  filter(!is.na(density))

# neighbourhood adult tree stem density
dens_control <- pcq_tests_dat %>% 
  filter(treatment == "control") %>% 
  pull(density)

dens_cut <- pcq_tests_dat %>% 
  filter(treatment == "cut") %>% 
  pull(density)

t.test(dens_control, dens_cut)

# neighbourhood adult tree height
hgt_control <- pcq_tests_dat %>% 
  filter(treatment == "control") %>% 
  pull(tree_height)

hgt_cut <- pcq_tests_dat %>% 
  filter(treatment == "cut") %>% 
  pull(tree_height)

t.test(hgt_control, hgt_cut)

# neighbourhood adult tree basal area
ba_control <- pcq_tests_dat %>% 
  filter(treatment == "control") %>% 
  pull(ba)

ba_cut <- pcq_tests_dat %>% 
  filter(treatment == "cut") %>% 
  pull(ba)

t.test(ba_control, ba_cut)

# soil plug core variables
plug_dat <- raw_juves %>% 
  filter(!is.na(plug_core_length)) %>% 
  mutate(across(c(plug_core_dry_wgt, plug_core_wet_wgt),
                ~ . -bag_wgt_coles_snack),
         moisture_content = (plug_core_wet_wgt - plug_core_dry_wgt) / 
           plug_core_wet_wgt,
         soil_volume = (pi * (15^2) * plug_core_length) / 1000,
         bulk_dens = plug_core_dry_wgt / soil_volume) 

# soil bulk density
bd_control <- plug_dat %>% 
  filter(treatment == "control") %>% 
  pull(bulk_dens)

bd_cut <- plug_dat %>% 
  filter(treatment == "cut") %>% 
  pull(bulk_dens)

t.test(bd_control, bd_cut)

# soil moisture content
betamod_soil_moist <- betareg(moisture_content ~ treatment, 
                              data = plug_dat)
summary(betamod_soil_moist)

# porewater salinity
sal_control <- raw_juves %>% 
  filter(treatment == "control") %>% 
  pull(salinity)

sal_cut <- raw_juves %>% 
  filter(treatment == "cut") %>% 
  pull(salinity)

t.test(sal_control, sal_cut)

#