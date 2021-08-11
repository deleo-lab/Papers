library(here)
library(tidyverse)
library(glmmTMB)
library(emmeans)
library(cowplot)

# download data
data <- readxl::read_xlsx(here("snail_habitat_data.xlsx"))
head(data)
data$Depth <- as.numeric(data$Depth)

# run model for Bulinus spp. snails 
hist(data, breaks=100)
# remove quadrats with snail counts > 150
bulin.nb <- glmmTMB(Bulinus ~ microhabitat*factor(LakeYN) + scale(Depth) + (1|Village/Site) + (1|FM),
                    ziformula = ~1,
                    data = data[data$Bulinus<150,],
                    family = "nbinom2")
summary(bulin.nb)

# run model for Biomphalaria pfeifferi snails
biomph.nb <- glmmTMB(Biomph ~ microhabitat*factor(LakeYN) + scale(Depth) + (1|Village/Site) + (1|FM),
                     ziformula = ~1,
                     data = data,
                     family = "nbinom2")
summary(biomph.nb)

# comparisons and plots
# estimate marginal means for snail counts in each microhabitat category, 
# then perform pairwise comparisons (contrasts)

# Bulinus spp.
bulin <- emmeans(bulin.nb, "microhabitat", by = "LakeYN", type = "response")
bulin.plot <- as.data.frame(bulin)
bulin.plot$Species <- "Bulinus spp."
contrast(bulin, method = "pairwise", type = "response", interaction = TRUE)

# Biomphalaria pfeifferi
biomph <- emmeans(biomph.nb, "microhabitat", by = "LakeYN", type = "response")
biomph.plot <- as.data.frame(biomph)
biomph.plot$Species <- "Biomphalaria pfeifferi"
contrast(biomph, method = "pairwise", type = "response", interaction = TRUE)

# plot
plot.data <- rbind(bulin.plot, biomph.plot)
plot.data <- plot.data %>% 
  mutate(ymin = response-SE,
         ymax = response+SE)
plot.data$ymin[plot.data$ymin<0] <- 0

ggplot(plot.data[plot.data$Species=="Bulinus spp.",], aes(x=reorder(microhabitat, -response), y=response/0.3677)) +
  geom_bar(stat = "identity") +
  geom_errorbar(width = 0.5, aes(ymin=ymin/0.3677, ymax=ymax/0.3677)) +
  scale_x_discrete(labels = c("Floating/\nSubmerged", "Emergent", "Water/\nMud")) +
  ylab(expression(italic("Bulinus")~'spp. density (square m)')) +
  xlab("") +
  facet_grid(LakeYN~., scales = "free", labeller = labeller(LakeYN = c("0" = "River", "1" = "Lake"))) +
  theme_cowplot(font_size = 18)

ggplot(plot.data[plot.data$Species=="Biomphalaria pfeifferi",], aes(x=reorder(microhabitat, -response), y=response/0.3677)) +
  geom_bar(stat = "identity") +
  geom_errorbar(width = 0.5, aes(ymin=ymin/0.3677, ymax=ymax/0.3677)) +
  scale_x_discrete(labels = c("Floating/\nSubmerged", "Emergent", "Water/\nMud")) +
  ylab(expression(italic("Biomphalaria pfeifferi")~'spp. density (square m)')) +
  xlab("") +
  facet_grid(LakeYN~., scales = "free", labeller = labeller(LakeYN = c("0" = "River", "1" = "Lake"))) +
  theme_cowplot(font_size = 18)



