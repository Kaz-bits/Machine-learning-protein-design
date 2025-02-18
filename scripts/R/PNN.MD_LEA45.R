# Packages
library(ggplot2)

# Load data
pnn_rg <- read.csv(file = "D:/MASTER_FILES/DATA/MD/cesar/rg.csv", header = TRUE)

# Plot
ggplot(data = pnn_rg, aes(x = Time, y = Ca/10)) +
  geom_line(color = "#e37100") +
  theme_bw() +
  labs(x = "Time (ns)", y = "Rg (nm)") +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12)) +
  coord_cartesian(xlim = c(0, 250)) +
  scale_x_continuous(expand = c(0, 0))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/rg_lea45_md_0.5nacl_250ns.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/rg_lea45_md_0.5nacl_250ns.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)





# Load data
pnn_ree <- read.csv(file = "D:/MASTER_FILES/DATA/MD/cesar/distance.csv", header = TRUE)

# Plot
ggplot(data = pnn_ree, aes(x = Time, y = distance/10)) +
  geom_line(color = "#e37100") +
  theme_bw() +
  labs(x = "Time (ns)", y = "Ree (nm)") +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12)) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 8)) +
  scale_x_continuous(expand = c(0, 0))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/ree_lea45_md_0.5nacl_250ns.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/ree_lea45_md_0.5nacl_250ns.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)





# Load data
pnn_ss <- read.csv(file = "D:/MASTER_FILES/DATA/MD/cesar/ss_per_frame.csv", header = TRUE)

# Filter data
pnn_ss1 <- pnn_ss[seq(0, nrow(pnn_ss), 10),]

# Plot
ggplot(data = pnn_ss1, aes(x = Frame/100)) +
  geom_line(aes(y = Disordered, color = "Disordered")) +
  geom_line(aes(y = Bend, color = "Bend")) +
  geom_line(aes(y = Turn, color = "Turn")) +
  geom_line(aes(y = a.helix, color = "α-Helix")) +
  theme_bw() +
  labs(x = "Time (ns)", y = "Percentage (%)") +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12)) +
  coord_cartesian(xlim = c(0, 250)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(name = "", values = c("#FF6E00", "#006666", "#490092", "#BB00BB"))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/ss_lea45_md_0.5nacl_250ns.pdf", device = "pdf", 
       width = 5, height = 3, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/ss_lea45_md_0.5nacl_250ns.png", device = "png", 
       width = 5, height = 3, units = "in", dpi = 450)





# Load data
pnn_rmsd <- read.csv(file = "D:/MASTER_FILES/DATA/MD/cesar/rmsd.csv", header = TRUE)


# Plot
ggplot(data = pnn_rmsd, aes(x = Time)) +
  geom_line(aes(y = Protein, color = "Protein")) +
  geom_line(aes(y = Backbone, color = "Backbone")) +
  geom_line(aes(y = C.alpha, color = "Cα")) +
  theme_bw() +
  labs(x = "Time (ns)", y = "RMSD") +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12)) +
  coord_cartesian(xlim = c(0, 250)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(name = "", values = c("#FF6E00", "#006666", "#490092"))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/rmsd_lea45_md_0.5nacl_250ns.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/rmsd_lea45_md_0.5nacl_250ns.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)





# Load data
pnn_rmsf <- read.csv(file = "D:/MASTER_FILES/DATA/MD/cesar/rmsf.csv", header = TRUE)


# Plot
ggplot(data = pnn_rmsf, aes(x = residue, rmsf)) +
  geom_line(color = "#490092") +
  theme_bw() +
  labs(x = "Residue (aa)", y = "RMSF") +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12)) +
  coord_cartesian(xlim = c(0, 158)) +
  scale_x_continuous(expand = c(0, 0))

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/rmsf_lea45_md_0.5nacl_250ns.pdf", device = "pdf", 
       width = 4, height = 2, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/rmsf_lea45_md_0.5nacl_250ns.png", device = "png", 
       width = 4, height = 2, units = "in", dpi = 450)
