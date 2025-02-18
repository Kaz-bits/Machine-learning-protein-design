# Packages
library(ggplot2)

# Load data
df_cahs <- read.csv(file = "D:/MASTER_FILES/DATA/CAHS/all_biosensors_gala.csv", 
                 header = TRUE)

# Load data
df_pcahs <- read.csv(file = "D:/MASTER_FILES/DATA/CAHS/IDRBS_CAHS_variants_sparrow_pSVM.csv", 
                     header = TRUE)

# Load data
df_sparrow <- read.csv(file = "D:/MASTER_FILES/DATA/IDRBS_library_sparrow_188.csv", 
                       header = TRUE)

# List of IDRs to remove
temp <- c("D7", "D8", "D9", "S1")

# Order data from df_cahs
df_cahs <- df_cahs[order(df_cahs$construct), ]

# Order data from df_pcahs
df_pcahs <- df_pcahs[order(df_pcahs$Construct), ]

# Remove IDRs from df_cahs
df_cahs <- df_cahs[!df_cahs$construct %in% temp, ]

# Add prediction from df_pcahs to df_cahs
df_cahs$Prediction <- df_pcahs$Prediction

# Obtain the quantiles from the 188-IDR based library
i <- quantile(df_sparrow$delta)[2]
h <- quantile(df_sparrow$delta)[4]

# Compute the median
j <- median(df_sparrow$delta)

# Plot
ggplot() +
  geom_vline(xintercept = i, lty = 2) +
  geom_vline(xintercept = j, lty = 2, color = "blue") +
  geom_vline(xintercept = h, lty = 2) +
  geom_point(data = df_cahs, aes(x = mean_delta, y = as.factor(construct), 
                                 colour = Prediction), 
             size = 4, shape = 21, 
             fill = "white", stroke = 1.5) +
  geom_errorbar(data = df_cahs, aes(x = mean_delta, 
                                       y = as.factor(construct), 
                                       xmin = mean_delta - sd_delta, 
                                       xmax = mean_delta + sd_delta, 
                                       color = Prediction), width = 0.2) +
  labs(x = expression(Delta * "FRET"), 
       y = "IDR") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank())  + 
  coord_cartesian(xlim = c(-0.1, 2.1)) +
  scale_color_manual(values = c("#FFA500", "#00D0D6"))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/delta_cahs_variants.pdf", device = "pdf", 
       width = 4.5, height = 3, units = "in", dpi = 400)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/delta_cahs_variants.png", device = "png", 
       width = 4.5, height = 3, units = "in", dpi = 400)

