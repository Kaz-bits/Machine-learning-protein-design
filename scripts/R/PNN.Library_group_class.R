# Packages
library(ggplot2)
library(ggcorrplot)

# load data base of the library with the organism properties
df_library <- readxl::read_xlsx(path = "D:/MASTER_FILES/DATA/IDR_Properties.xlsx", col_names = TRUE)

# Load data set from SPARROW (188 IDRs)
df_sparrow <- read.csv(file = "D:/MASTER_FILES/DATA/IDRBS_library_sparrow_188.csv", header = TRUE)

# Order the data sets by construct 
df_library <- df_library[order(as.integer(substr(df_library$entry, 7, 9))), ]
df_sparrow <- df_sparrow[order(df_sparrow$Construct), ]

# Filter data sets to obtain the 188 IDRs from the df_library data set
df_library <- df_library[as.integer(substr(df_library$entry, 7, 9)) %in% df_sparrow$Construct, ]

# Merge data sets
df_library_prop <- cbind(df_library[, -c(2, 3)], df_sparrow[, -c(1, 6)])

# Rearrange data set
df_library_prop <- df_library_prop[, c(1:27, 49:65, 28:48)]

# Rename some columns
names(df_library_prop)[1] <- "Construct"
names(df_library_prop)[21] <- "Abs_0.1_res"
names(df_library_prop)[22] <- "Abs_0.1_val"

# Change negative delta fret values to zero
df_library_prop$delta <- ifelse(df_library_prop$delta < 0, yes = 0, no = df_library_prop$delta)

# Save data frame
write.csv(x = df_library_prop, file = "D:/MASTER_FILES/DATA/IDRBS_library_sparrow_organisms.csv", 
          quote = FALSE, row.names = FALSE)

# Create sub-data sets by kigdom
temp_animalia <- df_library_prop[df_library_prop$kingdom == "Animalia", ] 
temp_fungi <- df_library_prop[df_library_prop$kingdom == "Fungi", ] 
temp_monera <- df_library_prop[df_library_prop$kingdom == "Monera", ] 
temp_plantae <- df_library_prop[df_library_prop$kingdom == "Plantae", ] 
temp_protista <- df_library_prop[df_library_prop$kingdom == "Protista", ] 

# Plot animalia
ggplot(data = temp_animalia, aes(x = helix_percentage, y = delta)) +
  geom_point(color = "#42B540") +
  geom_hline(yintercept = summary(temp_animalia$delta)[[2]], 
             lty = 2, linewidth = 0.5, alpha = 0.6) +
  geom_hline(yintercept = summary(temp_animalia$delta)[[5]], 
             lty = 2, linewidth = 0.5, alpha = 0.6) +
  theme_bw() +
  labs(x = "Helix propensity (%)", y = expression(Delta*"FRET"), title = "Animalia") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 2)) +
  scale_x_continuous(breaks = seq(0, 100, 20))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/fret_vs_helix_propensity_animalia.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/fret_vs_helix_propensity_animalia.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)




# Plot fungi
ggplot(data = temp_fungi, aes(x = helix_percentage, y = delta)) +
  geom_point(color = "#42B540") +
  geom_hline(yintercept = summary(temp_fungi$delta)[[2]], 
             lty = 2, linewidth = 0.5, alpha = 0.6) +
  geom_hline(yintercept = summary(temp_fungi$delta)[[5]], 
             lty = 2, linewidth = 0.5, alpha = 0.6) +
  theme_bw() +
  labs(x = "Helix propensity (%)", y = expression(Delta*"FRET"), title = "Fungi") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 2)) +
  scale_x_continuous(breaks = seq(0, 100, 20))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/fret_vs_helix_propensity_fungi.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/fret_vs_helix_propensity_fungi.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)




# Plot monera
ggplot(data = temp_monera, aes(x = helix_percentage, y = delta)) +
  geom_point(color = "#42B540") +
  geom_hline(yintercept = summary(temp_monera$delta)[[2]], 
             lty = 2, linewidth = 0.5, alpha = 0.6) +
  geom_hline(yintercept = summary(temp_monera$delta)[[5]], 
             lty = 2, linewidth = 0.5, alpha = 0.6) +
  theme_bw() +
  labs(x = "Helix propensity (%)", y = expression(Delta*"FRET"), title = "Monera") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 2)) +
  scale_x_continuous(breaks = seq(0, 100, 20))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/fret_vs_helix_propensity_monera.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/fret_vs_helix_propensity_monera.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)




# Plot plantae
ggplot(data = temp_plantae, aes(x = helix_percentage, y = delta)) +
  geom_point(color = "#42B540") +
  geom_hline(yintercept = summary(temp_plantae$delta)[[2]], 
             lty = 2, linewidth = 0.5, alpha = 0.6) +
  geom_hline(yintercept = summary(temp_plantae$delta)[[5]], 
             lty = 2, linewidth = 0.5, alpha = 0.6) +
  theme_bw() +
  labs(x = "Helix propensity (%)", y = expression(Delta*"FRET"), title = "Plantae") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 2)) +
  scale_x_continuous(breaks = seq(0, 100, 20))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/fret_vs_helix_propensity_plantae.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/fret_vs_helix_propensity_plantae.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)




# Plot protista
ggplot(data = temp_protista, aes(x = helix_percentage, y = delta)) +
  geom_jitter(color = "#42B540") +
  geom_hline(yintercept = summary(temp_protista$delta)[[2]], 
             lty = 2, linewidth = 0.5, alpha = 0.6) +
  geom_hline(yintercept = summary(temp_protista$delta)[[5]], 
             lty = 2, linewidth = 0.5, alpha = 0.6) +
  theme_bw() +
  labs(x = "Helix propensity (%)", y = expression(Delta*"FRET"), title = "Protista") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 2)) +
  scale_x_continuous(breaks = seq(0, 100, 20))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/fret_vs_helix_propensity_protista.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/fret_vs_helix_propensity_protista.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)








# Correlation matrix ----
# Extract data to be computed as a correlation matrix
temp_animalia_matx <- temp_animalia[, c(26:29, 32:44)]

# Create a correlatiox matrix
corr <- round(cor(temp_animalia_matx), 1)

# Plot animalia
ggcorrplot(corr, hc.order = TRUE, type = "lower", outline.color = "white", 
           lab = TRUE, insig = "blank")
