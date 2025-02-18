# Libraries
library(ggplot2)
library(ggridges)
library(dplyr)

# Functions
bins <- function(x) {
  if ((2 * IQR(x, na.rm = TRUE) / length(x)^(1/3)) > 0.05) {
    round(2 * IQR(x, na.rm = TRUE) / length(x)^(1/3), digits = 2)
  } else {
    round(2 * IQR(x, na.rm = TRUE) / length(x)^(1/3), digits = 2)
  }
}


# Name of the sequence folder
seq_name <- "lea45"

# List of concentration 
c1 <- "150"
c2 <- "200"
c3 <- "400"
c4 <- "600"
c5 <- "800"
c6 <- "1000"
c7 <- "1500"

# Combine vector of concentrations
vec_conc <- c(c1, c3, c7)

# Create an empty data frame
temp_all <- data.frame(matrix(nrow = 0, ncol = 6))

# Change column names
names(temp_all) <- c("frame", "rg", "re", "Delta", "S", "ID")

# Read multiple files
for (i in vec_conc) {
  # Read file one by one
  temp <- read.csv(file = paste0("D:/MASTER_FILES/DATA/calvados_1200ns/", seq_name, "/Rep1/", seq_name, "_NaCl_", i, "_IDRLab/", seq_name, "_NaCl_", i, "_IDRLab/", 
                                 "time_series_Rg_Ree_Delta_S.csv"), header = TRUE)
  
  
  # Add an ID for the concentration
  temp$treatment <- i
  
  # Add an ID for the sequence
  temp$ID <- seq_name
  
  
  # Change names of some columns
  names(temp)[1] <- "frame" 
  names(temp)[2] <- "rg" 
  names(temp)[3] <- "re"

  # Appen the files
  temp_all <- rbind(temp_all, temp)
  
}

# Data for 150 mM ----
# Extract data for 150 mM of NaCl
temp_150 <- temp_all[temp_all$treatment == "150", ]

# Compute the mean of Rg
temp_mean <- mean(temp_150$rg)

# Number of bins for histogram
temp_bins <- bins(temp_150$rg)

# Plot of distribution of Rg
ggplot(data = temp_150, aes(x = rg, y = ..ncount..)) +
  geom_histogram(fill = "#e37100", alpha = 0.8, color = "#000000", 
                 linewidth = 0.3, binwidth = temp_bins) +
  geom_vline(xintercept = temp_mean, lty = 2, linewidth = 0.7) +
  theme_bw() +
  labs(x = "Rg (nm)", y = "Frequency", title = "NaCl 0 mM") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(xlim = c(0, 8), ylim = c(0, 1.1)) +
  scale_y_continuous(expand = c(0, 0))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/lea45_150_rg_distribution.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/lea45_150_rg_distribution.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)


# Remove some frames
temp_150 <- temp_150[seq(1, nrow(temp_150), 75), ]

# Plot Rg series for each condition
ggplot(data = temp_150, aes(x = frame/12.5, y = rg)) +
  geom_line(color = "#0082cd", alpha = 0.8) +
  geom_hline(yintercept = temp_mean, lty = 2) +
  theme_bw() +
  labs(x = "Time (ns)", y = "Rg (nm)") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(xlim = c(0, 1200), ylim = c(0, 8)) +
  scale_x_continuous(expand = c(0, 0))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/lea45_150_rg_time_serie.pdf", device = "pdf", 
       width = 4, height = 1.5, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/lea45_150_rg_time_serie.png", device = "png", 
       width = 4, height = 1.5, units = "in", dpi = 450)


# Data for 400 mM ----
# Extract data for 400 mM of NaCl
temp_400 <- temp_all[temp_all$treatment == "400", ]

# Compute the mean of Rg
temp_mean <- mean(temp_400$rg)

# Number of bins for histogram
temp_bins <- bins(temp_400$rg)

# Plot of distribution of Rg
ggplot(data = temp_400, aes(x = rg, y = ..ncount..)) +
  geom_histogram(fill = "#e37100", alpha = 0.8, color = "#000000", 
                 linewidth = 0.3, binwidth = temp_bins) +
  geom_vline(xintercept = temp_mean, lty = 2, linewidth = 0.7) +
  theme_bw() +
  labs(x = "Rg (nm)", y = "Frequency", title = "NaCl 400 mM") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(xlim = c(0, 8), ylim = c(0, 1.1)) +
  scale_y_continuous(expand = c(0, 0))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/lea45_400_rg_distribution.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/lea45_400_rg_distribution.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)



# Remove some frames
temp_400 <- temp_400[seq(1, nrow(temp_400), 75), ]

# Plot Rg series for each condition
ggplot(data = temp_400, aes(x = frame/12.5, y = rg)) +
  geom_line(color = "#0082cd", alpha = 0.8) +
  geom_hline(yintercept = temp_mean, lty = 2) +
  theme_bw() +
  labs(x = "Time (ns)", y = "Rg (nm)") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(xlim = c(0, 1200), ylim = c(0, 8)) +
  scale_x_continuous(expand = c(0, 0))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/lea45_400_rg_time_serie.pdf", device = "pdf", 
       width = 4, height = 1.5, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/lea45_400_rg_time_serie.png", device = "png", 
       width = 4, height = 1.5, units = "in", dpi = 450)






# Data for 1500 mM ----
# Extract data for 1500 mM of NaCl
temp_1500 <- temp_all[temp_all$treatment == "1500", ]

# Compute the mean of Rg
temp_mean <- mean(temp_1500$rg)

# Number of bins for histogram
temp_bins <- bins(temp_1500$rg)

# Plot of distribution of Rg
ggplot(data = temp_1500, aes(x = rg, y = ..ncount..)) +
  geom_histogram(fill = "#e37100", alpha = 0.8, color = "#000000", 
                 linewidth = 0.3, binwidth = temp_bins) +
  geom_vline(xintercept = temp_mean, lty = 2, linewidth = 0.7) +
  theme_bw() +
  labs(x = "Rg (nm)", y = "Frequency", title = "NaCl 1500 mM") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(xlim = c(0, 8), ylim = c(0, 1.1)) +
  scale_y_continuous(expand = c(0, 0))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/lea45_1500_rg_distribution.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/lea45_1500_rg_distribution.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)



# Remove some frames
temp_1500 <- temp_1500[seq(1, nrow(temp_1500), 75), ]

# Plot Rg series for each condition
ggplot(data = temp_1500, aes(x = frame/12.5, y = rg)) +
  geom_line(color = "#0082cd", alpha = 0.8) +
  geom_hline(yintercept = temp_mean, lty = 2) +
  theme_bw() +
  labs(x = "Time (ns)", y = "Rg (nm)") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(xlim = c(0, 1200), ylim = c(0, 8)) +
  scale_x_continuous(expand = c(0, 0))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/lea45_1500_rg_time_serie.pdf", device = "pdf", 
       width = 4, height = 1.5, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/lea45_1500_rg_time_serie.png", device = "png", 
       width = 4, height = 1.5, units = "in", dpi = 450)




# Create scatterplot
ggplot(data = temp_md1, aes(x = Rg..nm., y = Ree..nm.)) + 
  geom_point(color  = "#8B00DF", alpha = 0.3) +
  theme_bw() + 
  labs(x = "Rg (nm)", y = "Ree (nm)") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        panel.grid = element_blank())




# Load data of CALVADOS at 200 mM
temp_md2 <- read.csv(file = dir400, header = TRUE)


# Add a new column with the NaCl concentration used
temp_md2$NaCl <- "400"

# Add temp_md with temp_all
temp_all <- rbind(temp_all, temp_md2)

# Create scatterplot
ggplot(data = temp_md2, aes(x = Rg..nm., y = Ree..nm.)) + 
  geom_point(color  = "#8B00DF", alpha = 0.8) +
  theme_bw() + 
  labs(x = "Rg (nm)", y = "Ree (nm)") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        panel.grid = element_blank())





# Load data of CALVADOS at 200 mM
temp_md6 <- read.csv(file = dir1500, header = TRUE)


# Add a new column with the NaCl concentration used
temp_md6$NaCl <- "1500"

# Add temp_md with temp_all
temp_all <- rbind(temp_all, temp_md6)

# Create scatterplot
ggplot(data = temp_md6, aes(x = Rg..nm., y = Ree..nm.)) + 
  geom_point(color  = "#8B00DF", alpha = 0.8) +
  theme_bw() + 
  labs(x = "Rg (nm)", y = "Ree (nm)") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        panel.grid = element_blank())





# Change name of columns of temp_all
names(temp_all) <- c("Frame", "rg", "re", "delta", "s", "Treatment")

# Obtain the mean of tge Rg from eahc condition
temp_lea45 <- temp_all %>% group_by(Treatment) %>% summarise(Mean_Rg = mean(rg), Mean_Ree = mean(re))

# Add an ID for the sequence
temp_lea45$seq <- "LEA4-5"

# Change order of columns of temp_sed1
temp_lea45 <- temp_lea45[, c(4, 1, 2, 3)]

# Change name of the columns
names(temp_lea45) <- c("ID", "Treatment", "rg", "re")

# Plot 
ggplot(data = temp_lea45, aes(x = Treatment, y = re, color = ID, group = ID)) +
  geom_point(size = 1.8) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  labs(x = "NaCl (mM)", y = "Rg (nm)") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank())



# Extract data for specific frames
temp_all$Frame[seq(0, nrow(temp_all), 150)]


# Plot Rg series for each condition
ggplot(data = temp_all[temp_all$Treatment == 150, ], aes(x = Frame, y = rg)) +
  geom_line(color = "#0082cd", alpha = 0.8) +
  theme_bw() +
  labs(x = "Frames", y = "Rg (nm)") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(2, 8))

# Plot Rg series for each condition
ggplot(data = temp_all[temp_all$Treatment == 1500, ], aes(x = Frame, y = rg)) +
  geom_line(color = "#0082cd", alpha = 0.8) +
  theme_bw() +
  labs(x = "Frames", y = "Rg (nm)") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(2, 8))



# Distribution of rg
ggplot(data = temp_all[temp_all$Treatment == 150, ], aes(x = rg)) +
  geom_density(color = "#0011df", linewidth = 0.6) +
  theme_bw() +
  labs(x = "Rg (nm)", y = "") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank())






# LEA4-5 full tail (80.8 ns) ---- 
# Name of the sequence folder
seq_name <- "lea45"

# Generate a vector with the NaCl concentrations used
vec_conc <- c("150", "200", "400", "600", "800", "1000", "1500") 

# Create an empty data frame
temp_all <- data.frame(matrix(nrow = 0, ncol = 6))

# Change column names
names(temp_all) <- c("frame", "rg", "re", "Delta", "S", "ID")

# Read multiple files
for (i in vec_conc) {
  # Read file one by one
  temp <- read.csv(file = paste0("D:/MASTER_FILES/DATA/calvados_autons/", seq_name, "/Rep1/", seq_name, "_NaCl_", i, 
                                 "_IDRLab/", seq_name, "_NaCl_", i, "_IDRLab/", 
                                 "time_series_Rg_Ree_Delta_S.csv"), header = TRUE)
  
  # Add an ID for the concentration
  temp$treatment <- i
  
  # Add an ID for the sequence
  temp$ID <- seq_name
  
  # Change names of some columns
  names(temp)[1] <- "frame" 
  names(temp)[2] <- "rg" 
  names(temp)[3] <- "re"
  
  # Appen the files
  temp_all <- rbind(temp_all, temp)
  
  
}

# Compute the mean Rg by each sequence variant
df_lea45 <- temp_all %>% group_by(ID, treatment) %>% summarise(Mean_Rg = mean(rg), Mean_Re = mean(re))

# Convert to factor the treatmen column
df_lea45$treatment <- factor(df_lea45$treatment, levels = c("150", "200", "400", "600", "800", "1000", "1500"))





# LEA4-5 without tail (100 ns) ---- 
# Name of the sequence folder
seq_name <- "lea45_wo_tail"

# Generate a vector with the NaCl concentrations used
vec_conc <- c("150", "200", "400", "600", "800", "1000", "1500") 

# Create an empty data frame
temp_all <- data.frame(matrix(nrow = 0, ncol = 6))

# Change column names
names(temp_all) <- c("frame", "rg", "re", "Delta", "S", "ID")

# Read multiple files
for (i in vec_conc) {
  # Read file one by one
  temp <- read.csv(file = paste0("D:/MASTER_FILES/DATA/calvados_100ns/", seq_name, "/Rep1/", seq_name, "_NaCl_", i, 
                                 "_IDRLab/", seq_name, "_NaCl_", i, "_IDRLab/", 
                                 "time_series_Rg_Ree_Delta_S.csv"), header = TRUE)
  
  # Add an ID for the concentration
  temp$treatment <- i
  
  # Add an ID for the sequence
  temp$ID <- seq_name
  
  # Change names of some columns
  names(temp)[1] <- "frame" 
  names(temp)[2] <- "rg" 
  names(temp)[3] <- "re"
  
  # Appen the files
  temp_all <- rbind(temp_all, temp)
  
  
}

# Compute the mean Rg by each sequence variant
df_lea45_wot <- temp_all %>% group_by(ID, treatment) %>% summarise(Mean_Rg = mean(rg), Mean_Re = mean(re))

# Convert to factor the treatmen column
df_lea45_wot$treatment <- factor(df_lea45_wot$treatment, levels = c("150", "200", "400", "600", "800", "1000", "1500"))

# Merge data sets
df_lea45_all <- rbind(df_lea45, df_lea45_wot)

# Plot
ggplot(data = df_lea45_all, aes(x = treatment, y = Mean_Rg, group = ID, color = ID)) +
  geom_point(size = 1.8) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  labs(x = "[NaCl] (mM)", y = "Mean Rg (nm)") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(2, 4)) + 
  scale_y_continuous(breaks = seq(2, 4, 0.5)) +
  scale_color_manual(name = NULL, labels = c("LEA4-5", "LEA4-5 no tail"),
                     values = c( "#42B540", "#e37100"))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/lea45_tail_calvados_rg.pdf", device = "pdf", 
       width = 5, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/lea45_tail_alvados_rg.png", device = "png", 
       width = 5, height = 2, units = "in", dpi = 450)



# Plot
ggplot(data = df_lea45_all, aes(x = treatment, y = Mean_Re, group = ID, color = ID)) +
  geom_point(size = 1.8) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  labs(x = "[NaCl] (mM)", y = "Mean Ree (nm)") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(5.5, 9.5)) + 
  scale_y_continuous(breaks = seq(5.5, 9.5, 0.8)) +
  scale_color_manual(name = NULL, labels = c("LEA4-5", "LEA4-5 no tail"),
                     values = c( "#42B540", "#e37100"))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/lea45_tail_calvados_re.pdf", device = "pdf", 
       width = 5, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/lea45_tail_alvados_re.png", device = "png", 
       width = 5, height = 2, units = "in", dpi = 450)

