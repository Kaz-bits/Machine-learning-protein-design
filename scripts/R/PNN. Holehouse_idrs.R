# Packages
library(ggplot2)

# Load data
df_idr <- read.csv(file = "D:/Documentos/MASTER_FILES/DATA/SVM_prediction/Holehouse/summary_GOOSE_mergedInfo.csv", 
                   header = TRUE)

# Load data
df_idr1 <- read.csv(file = "D:/Documentos/MASTER_FILES/DATA/SVM_prediction/Holehouse/svm_all_sequences_pSVM.csv", 
                    header = TRUE)

# List of IDRs to reomve
temp <- c(6, 7, 10, 19, 23, 24)

# Remove IDRs fro data set
df_idr <- df_idr[!df_idr$SequenceNumber %in% temp, ]
df_idr1 <- df_idr1[!df_idr1$fasta_header %in% temp, ]

# Order data
df_idr <- df_idr[order(df_idr$SequenceNumber), ]
df_idr1 <- df_idr1[order(df_idr1$fasta_header), ]

# Add classification data to the data set
df_idr$FRET_100mOsm <- df_idr1$FRET_100mOsm
df_idr$FRET_750mOsm <- df_idr1$FRET_750mOsm
df_idr$Prediction_2 <- df_idr1$Prediction_2

# Change numbers to factors
df_idr$SequenceNumber <- as.factor(df_idr$SequenceNumber)

# Plot
ggplot(data = df_idr, aes(x = deltaFRET_mean_750, y = SequenceNumber, color = Prediction_2)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(size = 4, shape = 21, fill = "white", stroke = 1.5) +
  geom_errorbar(aes(xmin = deltaFRET_mean_750 - deltaFRET_std_750, 
                    xmax = deltaFRET_mean_750 + deltaFRET_std_750)) +
  labs(x = expression(Delta * E[f]), 
       y = "IDR") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_color_manual(name = "Prediction", values = c("#FFA500", "#00D0D6"))
  

# Save plot
ggsave(filename = "D:/Documentos/MASTER_FILES/PLOTS/holehouse_idrs_svm.pdf", 
       device = "pdf", width = 5.4, height = 4.8, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/Documentos/MASTER_FILES/PLOTS/holehouse_idrs_svm.png", 
       device = "png", width = 5.4, height = 4.8, units = "in", dpi = 450)
