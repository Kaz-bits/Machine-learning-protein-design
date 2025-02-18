# Packages
library(ggplot2)

# Load data of FRET
df_fret <- read.csv(file = "D:/MASTER_FILES/DATA/all_fret_data_normalized.csv", header = TRUE)

# Load data from SPARROW
df_sparrow <- read.csv(file = "D:/MASTER_FILES/DATA/IDRBS_library_sparrow_188.csv", header = TRUE)

# Change the ID of the constructs from df_sparrow
df_sparrow$Construct <- as.integer(substr(df_sparrow$Construct, 7, 9))

# Fiter data
df_fret <- df_fret[df_fret$Construct %in% df_sparrow$Construct, ]

# Order by construct df_fret data set
df_fret <- df_fret[order(df_fret$Construct), ]

# Order by construct df_sparrow data set
df_sparrow <- df_sparrow[order(df_sparrow$Construct), ]

# Crate new data frame
df_fret_kappa <- data.frame(matrix(ncol = 7, nrow = 0))

# Change names of the columns of the new data frame
names(df_fret_kappa) <- c("DxAm.DxDm", "Construct", "kappa", "FCR", "rg", "re", "treatment")

# Iterate each IDR
for (i in df_fret$Construct) {
  
  # Extract FRET data from each biosensor
  temp <- as.data.frame(t(unname(df_fret[df_fret$Construct == i, ][2:ncol(df_fret)])))
  
  # Change name of the first column
  names(temp)[1] <- "DxAm.DxDm"
  
  # Add ID
  temp$Construct <- i
  
  # Add value of kappa
  temp$kappa <- df_sparrow[df_sparrow$Construct == i, ]$kappa
  
  # Add value of FCR
  temp$FCR <- df_sparrow[df_sparrow$Construct == i, ]$FCR
  
  # Add value of Rg
  temp$rg <- df_sparrow[df_sparrow$Construct == i, ]$rg
  
  # Add value of Ree
  temp$re <- df_sparrow[df_sparrow$Construct == i, ]$re
  
  # Add a column with NaCl concentrations
  temp$treatment <- c("150", "200", "400", "600", "800", "1000", "1500")
  
  # Append data frames
  df_fret_kappa <- rbind(df_fret_kappa, temp)
  
}

# Change column
df_fret_kappa <- df_fret_kappa[, c(2, 1, 3:ncol(df_fret_kappa))]

# Convert treatment to factor
df_fret_kappa$treatment <- factor(df_fret_kappa$treatment, levels = c("150", "200", "400", "600", "800", "1000", "1500"))

# Save IDRs to extract
temp <- c(151, 182, 034, 058, 122, 023, 088, 162, 041, 139, 084, 079, 109, 133, 
          036, 142, 039, 038, 110, 024, 064, 029, 096, 063, 171, 120, 076, 053, 
          015, 126, 022, 135, 099, 040, 047, 062, 166, 118, 014, 121, 005, 017, 
          173, 030, 040, 006, 094, 025, 125, 108, 128, 129, 141, 086, 097, 020)

# Extract data from sparrow
df_sparrow_test <- df_sparrow[df_sparrow$Construct %in% temp, ]

# Proof of concept ----
# Extract three IDRs
df_fret_test <- df_fret_kappa[df_fret_kappa$Construct == 025 |
                              df_fret_kappa$Construct == 099 |
                              df_fret_kappa$Construct == 135 |
                              df_fret_kappa$Construct == 096 |
                              df_fret_kappa$Construct == 039 |
                              df_fret_kappa$Construct == 006 |
                              df_fret_kappa$Construct == 023, ]


df_fret_test <- df_fret_kappa[df_fret_kappa$Construct == 025 |
                              df_fret_kappa$Construct == 099 |
                              df_fret_kappa$Construct == 135, ]

# Convert kappa values into rounded values
df_fret_test$kappa <- round(df_fret_test$kappa, digits = 2)

# Plot
ggplot(data = df_fret_test, aes(x = treatment, y = DxAm.DxDm, color = as.character(kappa), group = kappa)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(x = "[NaCl] (mM)", y = "DxAm/DxDm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_color_manual(name = "kappa", 
                     values = c("#ff7f0e", "#1f77b4", "#2ca02c"))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/kappa_low_high_FRET_response.pdf", device = "pdf", 
       width = 4.6, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/kappa_low_high_FRET_response.png", device = "png", 
       width = 4.6, height = 2, units = "in", dpi = 450)


# Convert FCR values into rounded values
df_fret_test$FCR <- round(df_fret_test$FCR, digits = 2)

# Plot
ggplot(data = df_fret_test, aes(x = treatment, y = DxAm.DxDm, color = as.character(FCR), group = FCR)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(x = "[NaCl] (mM)", y = "DxAm/DxDm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_color_manual(name = "FCR", 
                     values = c("#2ca02c", "#ff7f0e", "#1f77b4"))


# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/fcr_low_high_FRET_response.pdf", device = "pdf", 
       width = 4.6, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/fcr_low_high_FRET_response.png", device = "png", 
       width = 4.6, height = 2, units = "in", dpi = 450)


# Convert FCR values into rounded values
df_fret_test$rg <- round(df_fret_test$rg/10, digits = 2)

# Plot
ggplot(data = df_fret_test, aes(x = treatment, y = DxAm.DxDm, color = as.character(rg), group = rg)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(x = "[NaCl] (mM)", y = "DxAm/DxDm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_color_manual(name = "Rg (nm)", 
                     values = c("#ff7f0e", "#1f77b4", "#2ca02c"))
