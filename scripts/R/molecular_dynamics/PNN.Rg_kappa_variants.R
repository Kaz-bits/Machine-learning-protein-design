# Packages
library(ggplot2)
library(dplyr)

# load data
df_var_kappa <- read.csv(file = "D:/MASTER_FILES/DATA/pDN_var_kappa_sparrow_pSVM.csv", 
                         header = TRUE)

# De novo IDRs variants ----
# Generate a vector with the number of sequence variants used
vec_names <- c("var1", "var2", "var3", "var4", "var5", "var6", "var7", "var8")

# Generate a vector with the NaCl concentrations used
vec_conc <- c("150", "200", "400", "600", "800", "1000", "1500")

# Create an empty data frame
temp_all <- data.frame(matrix(nrow = 0, ncol = 6))
temp_all_1 <- data.frame(matrix(nrow = 0, ncol = 6))

# Change column names
names(temp_all) <- c("frame", "rg", "re", "Delta", "S", "ID")
names(temp_all_1) <- c("frame", "rg", "re", "Delta", "S", "ID")

## FCR = 0.4 ----
for (j in vec_conc) {
  for (i in vec_names) {
    # Read file one by one
    temp <- read.csv(file = paste0("D:/MASTER_FILES/DATA/calvados_100ns/", "variants/Rep1/", i, "_FCR_0.4_", j, "_IDRLab/", i, "_FCR_0.4_", j, "_IDRLab/", 
                                   "time_series_Rg_Ree_Delta_S.csv"), header = TRUE)
    
    # Add an ID for the concentration
    temp$treatment <- j
    
    # Add an ID for the sequence
    temp$ID <- i
    
    # Change names of some columns
    names(temp)[1] <- "frame" 
    names(temp)[2] <- "rg" 
    names(temp)[3] <- "re"
    
    # Appen the files
    temp_all <- rbind(temp_all, temp)
    
  }
}

# Compute the mean Rg by each sequence variant
df_mean_rg <- temp_all %>% group_by(ID, treatment) %>% summarise(Mean_Rg = mean(rg), Mean_Re = mean(re))

# Create vector with kappa values obtained by GOOSE
vec_kappa <- c(0.0291, 0.1315, 0.2690, 0.4053, 0.5054, 0.6513, 0.7847, 0.9915)

# Add values of kappa to each sequence
df_mean_rg$kappa <- rep(vec_kappa, each = 7)

# Order data 
df_mean_rg$treatment <- factor(df_mean_rg$treatment, levels = c("150", "200", "400", "600", "800", "1000", "1500"))

# Plot
ggplot(data = df_mean_rg, aes(x = kappa, y = Mean_Rg, group = treatment, 
                              color = treatment)) +
  geom_point(size = 1.8) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  labs(x = "kappa", y = "Mean Rg (nm)", title = "FCR = 0.4") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank(),
        legend.key.size = unit(0.48, "cm")) +
  coord_cartesian(ylim = c(1.5, 3.5)) + 
  scale_y_continuous(breaks = seq(1.5, 3.5, 0.5)) +
  scale_color_manual(name = "[NaCl] (mM)", 
                     values = c("#42b540", "#5daa35", "#789e2b", 
                                "#939320", "#ad8815", "#c87c0b", 
                                "#e37100"))


# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/kappa_variants_calvados_rg_FCR_0.4.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/kappa_variants_calvados_rg_FCR_0.4.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)



# Plot
ggplot(data = df_mean_rg, aes(x = kappa, y = Mean_Re, group = treatment, 
                              color = treatment)) +
  geom_point(size = 1.8) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  labs(x = "kappa", y = "Mean Re (nm)", title = "FCR = 0.4") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank(),
        legend.key.size = unit(0.48, "cm")) +
  coord_cartesian(ylim = c(1.5, 9.5)) + 
  scale_y_continuous(breaks = seq(1.5, 9.5, 2)) +
  scale_color_manual(name = "[NaCl] (mM)", 
                     values = c("#42b540", "#5daa35", "#789e2b", 
                                "#939320", "#ad8815", "#c87c0b", 
                                "#e37100"))


# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/kappa_variants_calvados_re_FCR_0.4.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/kappa_variants_calvados_re_FCR_0.4.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)







## FCR = 0.8 ----
for (j in vec_conc) {
  for (i in vec_names) {
    # Read file one by one
    temp <- read.csv(file = paste0("D:/MASTER_FILES/DATA/calvados_100ns/", "variants/Rep1/", i, "_FCR_0.8_", j, "_IDRLab/", i, "_FCR_0.8_", j, "_IDRLab/", 
                                   "time_series_Rg_Ree_Delta_S.csv"), header = TRUE)
    
    # Add an ID for the concentration
    temp$treatment <- j
    
    # Add an ID for the sequence
    temp$ID <- i
    
    # Change names of some columns
    names(temp)[1] <- "frame" 
    names(temp)[2] <- "rg" 
    names(temp)[3] <- "re"
    
    # Appen the files
    temp_all_1 <- rbind(temp_all_1, temp)
    
  }
}

# Compute the mean Rg by each sequence variant
df_mean_rg1 <- temp_all_1 %>% group_by(ID, treatment) %>% summarise(Mean_Rg = mean(rg), Mean_Re = mean(re))

# Create vector with kappa values obtained by GOOSE
vec_kappa <- c(0.0186, 0.1413, 0.2770, 0.3820, 0.4921, 0.6635, 0.7782, 0.9930)

# Add values of kappa to each sequence
df_mean_rg1$kappa <- rep(vec_kappa, each = 7)

# Order data 
df_mean_rg1$treatment <- factor(df_mean_rg$treatment, levels = c("150", "200", "400", "600", "800", "1000", "1500"))

# Plot
ggplot(data = df_mean_rg1, aes(x = kappa, y = Mean_Rg, group = treatment, 
                              color = treatment)) +
  geom_point(size = 1.8) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  labs(x = "kappa", y = "Mean Rg (nm)", title = "FCR = 0.8") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank(),
        legend.key.size = unit(0.48, "cm")) +
  coord_cartesian(ylim = c(1.5, 3.5)) + 
  scale_y_continuous(breaks = seq(1.5, 3.5, 0.5)) +
  scale_color_manual(name = "[NaCl] (mM)", 
                     values = c("#42b540", "#5daa35", "#789e2b", 
                                "#939320", "#ad8815", "#c87c0b", 
                                "#e37100"))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/kappa_variants_calvados_rg_FCR_0.8.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/kappa_variants_calvados_rg_FCR_0.8.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)



# Plot
ggplot(data = df_mean_rg1, aes(x = kappa, y = Mean_Re, group = treatment, 
                              color = treatment)) +
  geom_point(size = 1.8) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  labs(x = "kappa", y = "Mean Re (nm)", title = "FCR = 0.8") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank(),
        legend.key.size = unit(0.48, "cm")) +
  coord_cartesian(ylim = c(1.5, 9.5)) + 
  scale_y_continuous(breaks = seq(1.5, 9.5, 2)) +
  scale_color_manual(name = "[NaCl] (mM)", 
                     values = c("#42b540", "#5daa35", "#789e2b", 
                                "#939320", "#ad8815", "#c87c0b", 
                                "#e37100"))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/kappa_variants_calvados_re_FCR_0.8.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/kappa_variants_calvados_re_FCR_0.8.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)





# SVM prediction ----
## FCR = 0.4 ----
# Add the prediction feature to the df_mean_rg data set for FCR = 0.4
df_mean_rg$Prediction <- rep(df_var_kappa[df_var_kappa$FCR == 0.4, ]$Prediction, each = 7)

# Plot
ggplot(data = df_mean_rg, aes(x = kappa, y = Mean_Rg, group = treatment, 
                              color = Prediction)) +
  geom_point(size = 1.8) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  labs(x = "kappa", y = "Mean Rg (nm)", title = "FCR = 0.4") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank(),
        legend.key.size = unit(0.48, "cm")) +
  coord_cartesian(ylim = c(1.5, 3.5)) + 
  scale_y_continuous(breaks = seq(1.5, 3.5, 0.5)) +
  scale_color_manual(name = "Prediction", 
                     values = c("#FFA500", "#00D0D6"))


# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/kappa_variants_calvados_rg_FCR_0.4_pSVM.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/kappa_variants_calvados_rg_FCR_0.4_pSVM.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)



## FCR = 0.8 ----
# Add the prediction feature to the df_mean_rg data set for FCR = 0.8
df_mean_rg1$Prediction <- rep(df_var_kappa[df_var_kappa$FCR == 0.8, ]$Prediction, each = 7)

# Plot
ggplot(data = df_mean_rg1, aes(x = kappa, y = Mean_Rg, group = treatment, 
                               color = Prediction)) +
  geom_point(size = 1.8) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  labs(x = "kappa", y = "Mean Rg (nm)", title = "FCR = 0.8") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank(),
        legend.key.size = unit(0.48, "cm")) +
  coord_cartesian(ylim = c(1.5, 3.5)) + 
  scale_y_continuous(breaks = seq(1.5, 3.5, 0.5)) +
  scale_color_manual(name = "Prediction", 
                     values = c("#FFA500", "#00D0D6"))


# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/kappa_variants_calvados_rg_FCR_0.8_pSVM.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/kappa_variants_calvados_rg_FCR_0.8_pSVM.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)






# Natural IDRs variants ----
# Load all library data set
df_library <- read.csv(file = "D:/Documentos/MASTER_FILES/DATA/all_fret_data.csv", header = TRUE)

# Order data set by kappa values
df_library <- df_library[order(df_library$kappa), ]

# Extract some IDRs with different values of kappa (0 to 1)
temp_kappa <- df_library[c(1, 30, 105, 160, 180, 184, 185, 186), ]

# Extract FRET data by each concentration
vec_fret <- c(temp_kappa$NaCl_0_P1, temp_kappa$NaCl_200, temp_kappa$NaCl_400, temp_kappa$NaCl_600, 
              temp_kappa$NaCl_800, temp_kappa$NaCl_1000, temp_kappa$NaCl_1500)

# Create vector with NaCl concentratios used in the experiment
vec_conc_fret <- as.factor(rep(c(0, 200, 400, 600, 800, 1000, 1500), each = nrow(temp_kappa)))

# Add NaCl columns in a single column
df_mean_fret <- data.frame("Construct" = temp_kappa$Entry, "FRET_re" = vec_fret, 
                           "kappa" = temp_kappa$kappa, "treatment" = vec_conc_fret)


# Extract data for 0 NaCl mM and 1500 NaCl mM
temp_mean_fret <- df_mean_fret[df_mean_fret$treatment == 0 | df_mean_fret$treatment == 1500, ]

# Plot
ggplot(data = temp_mean_fret, aes(x = kappa, y = FRET_re, group = treatment, color = treatment)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(x = "kappa", y = "DxAm/DxDm") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(0.5, 2), xlim = c(0.01, 0.8)) + 
  scale_y_continuous(breaks = seq(0.5, 2, 0.5)) +
  scale_color_manual(name = "[NaCl] (mM)", 
                     values = c( "#42B540", "#e37100"))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/kappa_fret_experimental.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/kappa_fret_experimental.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)



# Plot
ggplot(data = df_mean_fret, aes(x = kappa, y = FRET_re, group = treatment, color = treatment)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(x = "kappa", y = "DxAm/DxDm") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank(), 
        legend.key.size = unit(0.48, "cm")) +
  coord_cartesian(ylim = c(0.3, 2), xlim = c(0.01, 0.8)) + 
  scale_y_continuous(breaks = seq(0.5, 2, 0.5)) +
  scale_color_manual(name = "[NaCl] (mM)", 
                     values = c("#42b540", "#5daa35", "#789e2b", 
                                "#939320", "#ad8815", "#c87c0b", 
                                "#e37100"))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/kappa_fret_experimental_all_conc.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/kappa_fret_experimental_all_conc.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)


# Round kappa values
df_mean_fret$kappa <- round(df_mean_fret$kappa, 2)

# Create a factor with the kappa values
df_mean_fret$kappa <- as.factor(df_mean_fret$kappa)

# Plot
ggplot(data = df_mean_fret, aes(x = kappa, y = FRET_re, fill = kappa)) +
       geom_boxplot(width = 0.5, show.legend = FALSE) +
       geom_point(alpha = 0.4, show.legend = FALSE) +
       theme_bw() +
       labs(x = "kappa", y = "DxAm/DxDm") +
       theme(axis.title = element_text(size = 14), 
             axis.text = element_text(size = 12),
             panel.grid = element_blank()) +
       coord_cartesian(ylim = c(0.5, 2)) + 
       scale_y_continuous(breaks = seq(0.5, 2, 0.5)) +
       scale_fill_manual(values = c("#42b540", "#59ab37", "#70a22e", "#879825", 
                                    "#9e8e1b", "#b58412", "#cc7b09", "#e37100"))

