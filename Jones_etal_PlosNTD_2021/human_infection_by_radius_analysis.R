library(here)
library(tidyverse)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(bbmle)
library(visreg)

# data has been summarized to protect patient privacy
# to discuss the possibility of obtaining de-identified data, please inquire with the senior author
# Giulio De Leo, deleo@stanford.edu

# use a loop to run model on all radii of interest
bands <- c(1,5,45,60,75,90,105,120)
names <- c("one", "five", "forty-five", "sixty", "seventy-five", "ninety", "one0five", "onetwenty")

# Schistosoma haematobium models
# first for infection presence
mods <- vector(mode="list", length=length(bands))
for (i in 1:length(bands)) {
  fit <- glmer(Sh ~ Class + sex + scale(Pop) + factor(year)
               + scale(FloatingVeg)
               + scale(length)
               + scale(width_shore) + scale(width_open)
               + circ_score*LakeYN
               + (1|VillageCode/ID),
               family = "binomial",
               control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
               data = data[data$distance==bands[[i]],]) # data_subset
  mods[[i]] <- fit
  names(mods) <- names
}

# aic comparison
AICtab(mods)

# now for egg burden
mods.egg <- vector(mode="list", length=length(bands))
for (i in 1:length(bands)) {
  fit <- glmmTMB(ShW ~ Class + sex + scale(Pop) + factor(year)
                 + scale(FloatingVeg)
                 + scale(length) 
                 + scale(width_shore) + scale(width_open)
                 + circ_score*LakeYN 
                 + (1|VillageCode/ID),
                 family = "nbinom2",
                 data = data[data_frame()$distance==bands[[i]],])
  mods.egg[[i]] <- fit
  names(mods.egg) <- names
}

AICtab(mods.egg)

# Schistosoma mansoni models
# first for infection presence
mods.SM <- vector(mode="list", length=length(bands))
for (i in 1:length(bands)) {
  fit <- glmer(Sm ~ Class + sex + scale(Pop) + factor(year)
               + scale(FloatingVeg)
               + scale(length)
               + scale(width_shore) + scale(width_open)
               + circ_score*LakeYN
               + (1|VillageCode/ID),
               family = "binomial",
               control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
               data = data[data_frame_()$distance==bands[[i]],]) # data_subset
  mods.SM[[i]] <- fit
  names(mods.SM) <- names
}

# aic comparison
AICtab(mods.SM)

# now for egg burden
mods.egg.SM <- vector(mode="list", length=length(bands))
for (i in 1:length(bands)) {
  fit <- glmmTMB(SmW ~ Class + sex + scale(Pop) + factor(year)
                 + scale(FloatingVeg)
                 + scale(length) 
                 + scale(width_shore) + scale(width_open)
                 + circ_score*LakeYN 
                 + (1|VillageCode/ID),
                 family = "nbinom2",
                 data = data_noflow_peak_150m[data_noflow_peak_150m$distance==bands[[i]],])
  mods.egg.SM[[i]] <- fit
  names(mods.egg.SM) <- names
}

AICtab(mods.egg.SM)

# AIC plot
############ plot AIC or log likelihood ############ 
distance <- data.frame(bands, names)
colnames(distance) <- c("distance", "Modnames")

sh_aic <- as.data.frame(aictab(mods)) %>% 
  mutate(species = "haematobium",
         outcome = "infection prevalence")
sh_aic <- left_join(sh_aic, distance, by="Modnames")

shW_aic <- as.data.frame(aictab(mods.egg)) %>% 
  mutate(species = "haematobium",
         outcome = "egg burden")
shW_aic <- left_join(shW_aic, distance, by="Modnames")

sm_aic <- as.data.frame(aictab(mods.SM)) %>% 
  mutate(species = "mansoni",
         outcome = "infection prevalence")
sm_aic <- left_join(sm_aic, distance, by="Modnames")

smW_aic <- as.data.frame(aictab(mods.egg.SM)) %>% 
  mutate(species = "mansoni",
         outcome = "egg burden")
smW_aic <- left_join(smW_aic, distance, by="Modnames")

plot_prev_aic <- rbind(sh_aic, sm_aic)
plot_egg_aic <- rbind(shW_aic, smW_aic)

plot_prev_aic$distance[plot_prev_aic$distance==5] <- 15
# add 5 to all distances for Sm so that overlap is more interpretable
plot_prev_aic$distance[plot_prev_aic$species=="mansoni"] <- plot_prev_aic$distance[plot_prev_aic$species=="mansoni"] + 5
# now for IRR
plot_egg_aic$distance[plot_egg_aic$distance==5] <- 15
plot_egg_aic$distance[plot_egg_aic$species=="mansoni"] <- plot_egg_aic$distance[plot_egg_aic$species=="mansoni"] + 5

head(plot_prev_aic)
dodger = position_dodge(width = 7)
p_aic <- ggplot(plot_prev_aic, aes(distance, as.numeric(Delta_AICc), color=species, shape=species)) +  
  geom_hline(yintercept = 0, linetype = "dotted", size = 2) +
  geom_point(position = dodger, size = 6) +
  scale_color_grey(start=0.2, end=0.5) +
  #ylim(c(-1,15)) +
  scale_x_continuous(breaks = c(1, 5, 45, 60, 75, 90, 105, 120), minor_breaks = NULL) +
  facet_grid(species~.) +
  theme_cowplot(); p_aic

# now for egg burden
p_egg_aic <- ggplot(plot_egg_aic, aes(distance, as.numeric(Delta_AICc), color=species, shape=species)) +  
  geom_hline(yintercept = 0, linetype = "dotted", size = 2) +
  geom_point(position = dodger, size = 6) +
  scale_color_grey(start=0.2, end=0.5) +
  #ylim(c(-5,45)) +
  scale_x_continuous(breaks = c(1, 5, 45, 60, 75, 90, 105, 120), minor_breaks = NULL) +
  facet_grid(species~.) +
  theme_cowplot(); p_egg_aic

# Plot S. haematobium and S. mansoni egg burden vs. area of non-emergent vegetation measured within a 90 m sampling radius
# raw data for plots
pdat <- data[data$distance==90 & !is.na(data$sex),]

# S. haematobium
# get prediction line
float.90.egg.f <- ggpredict(mods.egg[[6]], terms = c("FloatingVeg", "sex [F]", "LakeYN [0]", "year [2017]"), interval="confidence")
# plot with raw egg burden data
ggplot(float.90.egg.f, aes(x, predicted)) + 
  geom_jitter(data = pdat, aes(FloatingVeg, ShW+1, group=FloatingVeg, color=ifelse(sex=="F", "darkgrey", "lightgrey")), width=0.03, height=0.1, size=5) +
  scale_color_identity() +
  geom_line(size=5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, group=group), alpha = .1) + 
  scale_y_log10() +
  scale_x_log10() +
  labs(title="90 m snail habitat sampling radius") +
  theme_classic()

# S. mansoni
# plot raw egg burden data
# no prediction line because no statistically significant association between egg burden and non-emergent vegetation at 90m m
ggplot(data = pdat, aes(FloatingVeg, SmW+1, group=FloatingVeg, color=ifelse(sex=="F", "darkgrey", "lightgrey"))) +
  geom_jitter(width=0.03, height=0.1, size = 5) +
  scale_color_identity() +
  scale_x_log10() +
  scale_y_log10() +
  labs(title="90 m snail habitat sampling radius") +
  theme_classic()

# plot prediction lines for infection presence versus water access site shoreline width, and circumscription on lake versus river settings
# water access site shoreline width
# S. haematobium
visreg(mods[[6]], "width_shore", line=list(col="black"), points=list(size=2), gg=TRUE)
# S. mansoni 
visreg(mods.SM[[6]], "width_shore", line=list(col="black"), points=list(size=2), gg=TRUE)

# circumscription 
# S. haematobium
visreg(mods[[6]], "circ_score", by = "LakeYN", line=list(col="black"), points=list(size=2), gg=TRUE)
# S. mansoni 
visreg(mods.SM[[6]], "circ_score", by = "LakeYN", line=list(col="black"), points=list(size=2), gg=TRUE)

#################### #################### #################### #################### 
# Interaction model, to test for differences between S. haematobium and S. mansoni 
# get data in long format
data.sh <- data %>% 
  select(-Sm, -SmW) %>% 
  rename(schisto.prev = Sh,
         schisto.burd = ShW) %>% 
  mutate(spp = "haematobium")

data.sm <- data %>% 
  select(-Sh, -ShW) %>% 
  rename(schisto.prev = Sm,
         schisto.burd = SmW) %>% 
  mutate(spp = "mansoni")

data.long <- rbind(data.sh, data.sm)
data.long$year <- as.factor(data.long$year)

bands <- c(1,5,45,60,75,90,105,120)
names <- c("one", "five", "forty-five", "sixty", "seventy-five", "ninety", "one0five", "onetwenty")

# first for infection presence
mods <- vector(mode="list", length=length(bands))
for (i in 1:length(bands)) {
  fit <- glmer(schisto.prev ~ Class + sex + scale(Pop) + factor(year) 
               + scale(FloatingVeg)*spp #*flow_section
               + scale(length) 
               + scale(width_shore)*spp + scale(width_open)
               + circ_score*LakeYN*spp 
               + (1|VillageCode/ID),
               family = "binomial",
               control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
               data = data.long[data.long$distance==bands[[i]],])
  mods[[i]] <- fit
  names(mods) <- names
}

# aic comparison
AICtab(mods)

# now for egg burden
mods.egg <- vector(mode="list", length=length(bands))
for (i in 1:length(bands)) {
  fit <- glmmTMB(schisto.burd ~ Class + sex + scale(Pop) + factor(year) 
                 + scale(FloatingVeg)*spp #*flow_section
                 + scale(length) 
                 + scale(width_shore)*spp + scale(width_open)
                 + circ_score*LakeYN*spp 
                 + (1|VillageCode/ID),
                 family = "nbinom2",
                 data = data.long[data.long$distance==bands[[i]],])
  mods.egg[[i]] <- fit
  names(mods.egg) <- names
}
AICtab(mods.egg)




