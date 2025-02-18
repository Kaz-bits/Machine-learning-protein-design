# Packages
library(ggplot2)
library(dplyr)
library(car)

# Load pLDDT data set from library
df_plddt <- read.csv("D:/MASTER_FILES/DATA/IDRBS_DN_SVM_goose_890_all_sequences_pLDDT.csv", header = FALSE)

# Load data set from SPARROW (188 IDRs)
df_sparrow <- read.csv(file = "D:/MASTER_FILES/DATA/IDRBS_DN_SVM_goose_890_all.csv", header = TRUE)

# Change the name of the categtory of FRET response
df_sparrow$Prediction <- ifelse(df_sparrow$Prediction == "Baja", "Low", "High")

# Create a new data frame
df_all_plddt <- data.frame(matrix(nrow = 0, ncol = 5))

# Change name sof columns of new data frame
names(df_all_plddt) <- c("construct", "very_low", "low", "high", "very_high")

# Iterate over each construct or IDR
for (i in 1:nrow(df_plddt)) {
  
  # Select row
  temp <- unname((df_plddt[i, -1]))
  
  # Generate vector of plddt values only
  temp <- as.numeric(temp[!is.na(temp)])
  
  # Data below or equal to 50 of plddt (very low)
  plddt_vlow <- ((length(temp[temp <= 50])) * 100) / length(temp)
  
  # Data between 50 and 70 of plddt (low)
  plddt_low <- ((length(temp[temp > 50 & temp < 70])) * 100) / length(temp)
  
  # Data between 70 and 90 of plddt (high)
  plddt_high <- ((length(temp[temp > 70 & temp < 90])) * 100) / length(temp)
  
  # Data above or equal to 90 of plddt (very high)
  plddt_vhigh <- ((length(temp[temp >= 90])) * 100) / length(temp)
  
  # Merge vector
  temp <- data.frame("construct" = i, 
                     "very_low" = round(plddt_vlow, 4), 
                     "low" = round(plddt_low, 4), 
                     "high" = round(plddt_high, 4), 
                     "very_high" = round(plddt_vhigh, 4))
  
  # Add to data frame
  df_all_plddt <- rbind(df_all_plddt, temp)
  
}

# Add data to df_sparrow
df_library_all <- cbind(df_sparrow, df_all_plddt[, -1])


# Determine probablity values for high confident values
a <- qnorm(p = 0.95, mean = mean(df_library_all$high), sd(df_library_all$high)) # 69.02
b <- qnorm(p = 0.95, mean = mean(df_library_all$very_high), sd(df_library_all$very_high)) # 9.97

# Extract data using the probability values
temp_library <- df_library_all[df_library_all$high >= a | df_library_all$very_high >= b, ]


temp_library[temp_library$Prediction == "High", ]$Construct
temp_library[temp_library$Prediction == "High", ]$very_low
temp_library[temp_library$Prediction == "High", ]$low
temp_library[temp_library$Prediction == "High", ]$high
temp_library[temp_library$Prediction == "High", ]$very_high

temp1 <- temp_library[temp_library$Prediction == "Low", ]


# Create a new data frame
df_all_plddt1 <- data.frame(matrix(nrow = 0, ncol = 2))

# Change name sof columns of new data frame
names(df_all_plddt1) <- c("construct", "pLDDT")

# Iterate over each construct or IDR
for (i in 1:nrow(df_plddt)) {
  
  # Select row
  temp <- unname((df_plddt[i, -1]))
  
  # Generate vector of plddt values only
  temp <- as.numeric(temp[!is.na(temp)])
  
  # Compute the mean of pLDDT 
  temp_plddt <- mean(temp)
  
  # Merge vector
  temp <- data.frame("construct" = i, 
                     "pLDDT" = temp_plddt)
  
  # Add to data frame
  df_all_plddt1 <- rbind(df_all_plddt1, temp)
  
}

# Compute densities
densidad <- ggplot2::ggplot_build(ggplot(df_all_plddt1, aes(x = pLDDT)) + 
                                    geom_density())$data[[1]]

# Normalize density
densidad$scaled_density <- densidad$density / max(densidad$density)

# Plot
ggplot(data = densidad, aes(x = x, y = scaled_density)) +
  # pLDDT > 90
  geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = 90, xmax = Inf), 
            fill = "#868dff", alpha = 1/8) +
  # 90 > pLDDT > 70
  geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = 70, xmax = 90),
            fill = "#95edff", alpha = 1/8) +
  # 50 > pLDDT < 70
  geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = 50, xmax = 70),
            fill = "#f8ff9e", alpha = 1/6) +
  # pLDDT < 50
  geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = -Inf, xmax = 50, alpha = 0.1),
            fill = "#ffbb6f", alpha = 1/6) +
  geom_line(linewidth = 0.8) +
  theme_bw() +
  labs(x = "pLDDT score (%)", y = "Normalized density") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(0, 1.1)) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = seq(0, 1.2, 0.2)) +
  scale_x_continuous(expand = c(0, 0))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/plddt_density_plot_all_DN.pdf", device = "pdf", 
       width = 4, height = 2.5, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/plddt_density_plot_all_DN.png", device = "png", 
       width = 4, height = 2.5, units = "in", dpi = 450)






