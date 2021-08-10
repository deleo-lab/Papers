library(here)
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(cowplot)

# download data
data <- readxl::read_xlsx(here("deep_water_snail_data.xlsx"))
head(data)

# run model for Bulinus spp. snails 
bul.test <- glmmTMB(COUNT_bulinus ~ distance + LakeYN + (1|village/site_code) + offset(log(quad_area)), family = "nbinom2", data = data)
summary(bul.test) 
# model diagnostics
res = simulateResiduals(bul.test)
plot(res, rank = T)

# run model for Biomphalaria pfeifferi snails
biomph.test <- glmmTMB(COUNT_biomph ~ distance + LakeYN + (1|village/site_code) + offset(log(quad_area)), family = "nbinom2", data = data)
summary(biomph.test)  # slightly lower AIC without village as a random effect in the level above site_code
# model diagnostics
res2 = simulateResiduals(biomph.test)
plot(res2, rank = T) 

# comparisons and plots
# estimate marginal means for snail counts in each distance category, 
# then perform pairwise comparisons (contrasts)

# Bulinus spp.
bulin <- emmeans(bul.test, "distance", type = "response") 
bulin.plot <- as.data.frame(bulin)
bulin.plot$Species <- "Bulinus spp."
contrast(bulin,  adjust = "bonferroni", method = "pairwise", type = "response")

# plot
ggplot(bulin.plot, aes(x=distance, y=response)) +
  geom_bar(stat = "identity") +
  geom_errorbar(width = 0.5, aes(ymin=lower.CL, ymax=upper.CL)) +
  #scale_x_discrete(labels = c("...")) +
  ylab(expression(italic("Bulinus")~'spp. density (square m)')) +
  xlab("") +
  theme_cowplot(font_size = 16)

# Biomphalaria pfeifferi
biomph <- emmeans(biomph.test, "distance", type = "response")
biomph.plot <- as.data.frame(biomph)
biomph.plot$Species <- "Biomphalaria pfeifferi"
contrast(biomph, adjust = "bonferroni", method = "pairwise", type = "response")

# plot
ggplot(biomph.plot, aes(x=distance, y=response)) +
  geom_bar(stat = "identity") +
  geom_errorbar(width = 0.5, aes(ymin=lower.CL, ymax=upper.CL)) +
  #scale_x_discrete(labels = c("...")) +
  ylab(expression(italic("Biomphalaria")~"density (square m)")) +
  xlab("") +
  theme_cowplot(font_size = 16)





