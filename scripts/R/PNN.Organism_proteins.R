# Packages
library(ggplot2)

# Load data
df_prop <- read.csv(file = "D:/Documentos/MASTER_FILES/DATA/IDR_Properties.csv", header = TRUE)

# Compute quantity of IDRs per organism
table(df_prop$organism)[order(table(df_prop$organism), decreasing = TRUE)]

# Plot
ggplot(data = df_prop, aes(y = organism, fill = organism)) +
  geom_bar(show.legend = FALSE) +
  theme_classic() +
  labs(x = "Count", y = "Organism") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        axis.text.y = element_text(face = "italic")) +
  coord_cartesian(xlim = c(0, 62)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 60, 10))


# Save plot
ggsave(filename = "D:/Documentos/MASTER_FILES/PLOTS/PDF/organism_distribution.pdf", device = "pdf", 
       width = 6.5, height = 6.8, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/Documentos/MASTER_FILES/PLOTS/PNG/frganism_distribution.png", device = "png", 
       width = 6.5, height = 6.8, units = "in", dpi = 450)



# Compute quaqntity of proteins by name
table(df_prop$protein_name)[order(table(df_prop$protein_name), decreasing = TRUE)]

# Plot
ggplot(data = df_prop, aes(y = protein_name, fill = protein_name)) +
  geom_bar(show.legend = FALSE) +
  theme_classic() +
  labs(x = "Count", y = "Protein name") +
  theme(axis.title = element_text(size = 50), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 50)) +
  coord_cartesian(xlim = c(0, 18)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 20, 5))


# Save plot
ggsave(filename = "D:/Documentos/MASTER_FILES/PLOTS/PDF/protein_distribution.pdf", device = "pdf", 
       width = 20, height = 28, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/Documentos/MASTER_FILES/PLOTS/PNG/protein_distribution.png", device = "png", 
       width = 20, height = 28, units = "in", dpi = 450)
