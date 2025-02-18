# Packages
library(ggplot2)

# Load data set
df_delta <- read.csv(file = "D:/MASTER_FILES/DATA/all_fret_data_no_normalized.csv", header = TRUE)

# Load SPARROW data set
df_sparrow <- read.csv(file = "D:/MASTER_FILES/DATA/IDRBS_library_sparrow_188.csv", header = TRUE)

# Filter data sets to get the sequences to be analyzed
df_delta <- df_delta[df_delta$Construct %in% df_sparrow$Construct, ]

# Arrange data sets by construct name
df_delta <- df_delta[order(df_delta$Construct), ]
df_sparrow <- df_sparrow[order(df_sparrow$Construct), ]

# Remove column fo construct from df_delta
df_delta <- df_delta[, -1]

# Merge data sets 
df_all <- cbind(df_sparrow, df_delta)

# Plot
ggplot(data = df_all, aes(x = NaCl_1500, y = NaCl_0, color = Response)) +
  geom_point() +
  theme_bw() +
  labs(x = "Final FRET", y = "Basal FRET") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        panel.grid = element_blank())


# Compute pearson correlation
a <- cor.test(df_all$NaCl_1500, df_all$NaCl_0, method = "pearson")[[4]]

# Compute correlation t value
b <- cor.test(df_all$NaCl_1500, df_all$NaCl_0, method = "pearson")[[3]]

# Plot
ggplot(data = df_all, aes(x = NaCl_1500, y = NaCl_0)) +
  geom_point() +
  geom_smooth(lty = 2, se = FALSE, linewidth = 0.8) +
  theme_bw() +
  labs(x = bquote("FRET"["NaCl 1.5M"]), y = bquote("FRET"["NaCl 0M"]), 
       title = paste("r = ", round(a, 3), ";", "p = ", signif(b, 3))) +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(0.5, 2)) +
  scale_y_continuous(breaks = seq(0.5, 2, 0.5))

# Save plots
ggsave(filename = "D:/MASTER_FILES/PLOTS/final_vs_basal_fret.pdf", 
       device = "pdf", width = 3.5, height = 3.5, units = "in", dpi = 450)

ggsave(filename = "D:/MASTER_FILES/PLOTS/final_vs_basal_fret.png", 
       device = "png", width = 3.5, height = 3.5, units = "in", dpi = 450)



# Extract data from NaCl of each condition
temp_nacl <- c(df_all[, "NaCl_0"], 
               df_all[, "NaCl_200"], 
               df_all[, "NaCl_400"], 
               df_all[, "NaCl_600"], 
               df_all[, "NaCl_800"], 
               df_all[, "NaCl_1000"], 
               df_all[, "NaCl_1500"]) 

# Create a vector with the NaCl concentrations used in the experiment
temp_cond <- rep(c("0", "200", "400", "600", "800", "1000", "1500"), each = nrow(df_all))

# Create a new data frame to transform the data to a large format
df_nacl <- data.frame("Construct" = as.factor(df_all$Construct), 
                      "DxDm.DxAm" = temp_nacl, 
                      "Condition" = as.factor(temp_cond))

# Order levels in the construct column
df_nacl$Condition <- factor(df_nacl$Condition, levels = c("0", "200", "400", "600", "800", "1000", "1500"))

# Plot
ggplot(data = df_nacl, aes(x = Condition, y = DxDm.DxAm, group = Construct)) +
  geom_point(alpha = 0.3) +
  geom_line(alpha = 0.3) +
  theme_bw() +
  labs(x = "[NaCl] (mM)", y = "DxDm/DxAm") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        panel.grid = element_blank())


# Plot
ggplot(data = df_nacl[df_nacl$Condition == 0 | df_nacl$Condition == 1500, ], 
       aes(x = Condition, y = DxDm.DxAm, group = Construct)) +
  geom_point(alpha = 0.3) +
  geom_line(alpha = 0.3) +
  theme_bw() +
  labs(x = "[NaCl] (mM)", y = "DxDm/DxAm") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        panel.grid = element_blank())




# Load data set from AGADIR of the library
df_agadir <- read.delim(file = "D:/MASTER_FILES/DATA/IDRBS_Library_200_agadir.txt", 
                        header = FALSE)

# Order df_all based on the sequences
df_all <- df_all[order(df_all$sequence, decreasing = FALSE), ]

# Order df_agadir based on the sequences
df_agadir <- df_agadir[order(df_agadir$V2, decreasing = FALSE), ]

# Extract data from df_agadir
df_agadir <- df_agadir[df_agadir$V2 %in% df_all$sequence, ]

# Change name of the first column from df_agadir
names(df_agadir)[3] <- "Agadir"

# Merge data sets
df_all <- cbind(df_all, df_agadir[3])

# Save data set
write.csv(x = df_all, file = "all_delta_agadir_library_188.csv")


# Plot
ggplot(data = df_all, aes(x = Agadir, y = delta)) +
  geom_point()

