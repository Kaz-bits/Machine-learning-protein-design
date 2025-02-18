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
seq_name <- "seq1_FCR_0.4_high"

# List of concentration 
c1 <- "150"
c2 <- "200"
c3 <- "400"
c4 <- "600"
c5 <- "800"
c6 <- "1000"
c7 <- "1500"

# Combine vector of concentrations
vec_conc <- c(c1, c7)

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
ggsave(filename = "D:/MASTER_FILES/PLOTS/seq1_150_rg_distribution.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/seq1_150_rg_distribution.png", device = "png", 
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
ggsave(filename = "D:/MASTER_FILES/PLOTS/seq1_150_rg_time_serie.pdf", device = "pdf", 
       width = 4, height = 1.5, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/seq1_150_rg_time_serie.png", device = "png", 
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
ggsave(filename = "D:/MASTER_FILES/PLOTS/seq1_400_rg_distribution.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/seq1_400_rg_distribution.png", device = "png", 
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
ggsave(filename = "D:/MASTER_FILES/PLOTS/seq1_400_rg_time_serie.pdf", device = "pdf", 
       width = 4, height = 1.5, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/seq1_400_rg_time_serie.png", device = "png", 
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
ggsave(filename = "D:/MASTER_FILES/PLOTS/seq1_1500_rg_distribution.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/seq1_1500_rg_distribution.png", device = "png", 
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
ggsave(filename = "D:/MASTER_FILES/PLOTS/seq1_1500_rg_time_serie.pdf", device = "pdf", 
       width = 4, height = 1.5, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/seq1_1500_rg_time_serie.png", device = "png", 
       width = 4, height = 1.5, units = "in", dpi = 450)

