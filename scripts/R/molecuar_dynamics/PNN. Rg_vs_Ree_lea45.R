# Libraries
library(ggplot2)
library(ggridges)
library(dplyr)

# Add directory for simulation at 200 mM of NaCl
dir200 <- "D:/MASTER_FILES/DATA/sed1_NaCl_200_IDRLab/sed1_NaCl_200_IDRLab/time_series_Rg_Ree_Delta_S.csv"

# Load data of CALVADOS at 200 mM
temp_md1 <- read.csv(file = dir200, header = TRUE)

# Add a new column with the NaCl concentration used
temp_md1$NaCl <- "200"

# Save in as a new data frame
temp_all <- temp_md1

# Create scatterplot
ggplot(data = temp_md1, aes(x = Rg..nm., y = Ree..nm.)) + 
  geom_point(color  = "#8B00DF", alpha = 0.8) +
  theme_bw() + 
  labs(x = "Rg (nm)", y = "Ree (nm)") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        panel.grid = element_blank()) +
  coord_cartesian(xlim = c(1.8, 5.3), ylim = c(2, 16))

# Save plots
ggsave(filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/IDR DESIGN/PDF/rg_vs_re_seq1_200.pdf", 
       device = "pdf", width = 3.5, height = 4, units = "in", dpi = 450)

ggsave(filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/IDR DESIGN/PNG/rg_vs_re_seq1_200.png", 
       device = "png", width = 3.5, height = 4, units = "in", dpi = 400)





# Add directory for simulation at 200 mM of NaCl
dir400 <- "D:/MASTER_FILES/DATA/sed1_NaCl_400_IDRLab/sed1_NaCl_400_IDRLab/time_series_Rg_Ree_Delta_S.csv"

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

# Save plots
ggsave(filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/IDR DESIGN/PDF/rg_vs_re_seq1_400.pdf", 
       device = "pdf", width = 3.5, height = 4, units = "in", dpi = 450)

ggsave(filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/IDR DESIGN/PNG/rg_vs_re_seq1_400.png", 
       device = "png", width = 3.5, height = 4, units = "in", dpi = 400)





# Add directory for simulation at 200 mM of NaCl
dir600 <- "D:/MASTER_FILES/DATA/sed1_NaCl_600_IDRLab/sed1_NaCl_600_IDRLab/time_series_Rg_Ree_Delta_S.csv"

# Load data of CALVADOS at 200 mM
temp_md3 <- read.csv(file = dir600, header = TRUE)

# Add a new column with the NaCl concentration used
temp_md3$NaCl <- "600"

# Add temp_md with temp_all
temp_all <- rbind(temp_all, temp_md3)

# Create scatterplot
ggplot(data = temp_md3, aes(x = Rg..nm., y = Ree..nm.)) + 
  geom_point(color  = "#8B00DF", alpha = 0.8) +
  theme_bw() + 
  labs(x = "Rg (nm)", y = "Ree (nm)") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        panel.grid = element_blank())

# Save plots
ggsave(filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/IDR DESIGN/PDF/rg_vs_re_seq1_600.pdf", 
       device = "pdf", width = 3.5, height = 4, units = "in", dpi = 450)

ggsave(filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/IDR DESIGN/PNG/rg_vs_re_seq1_600.png", 
       device = "png", width = 3.5, height = 4, units = "in", dpi = 400)






# Add directory for simulation at 200 mM of NaCl
dir800 <- "D:/MASTER_FILES/DATA/sed1_NaCl_800_IDRLab/sed1_NaCl_800_IDRLab/time_series_Rg_Ree_Delta_S.csv"

# Load data of CALVADOS at 200 mM
temp_md4 <- read.csv(file = dir800, header = TRUE)


# Add a new column with the NaCl concentration used
temp_md4$NaCl <- "800"

# Add temp_md with temp_all
temp_all <- rbind(temp_all, temp_md4)

# Create scatterplot
ggplot(data = temp_md4, aes(x = Rg..nm., y = Ree..nm.)) + 
  geom_point(color  = "#8B00DF", alpha = 0.8) +
  theme_bw() + 
  labs(x = "Rg (nm)", y = "Ree (nm)") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        panel.grid = element_blank())

# Save plots
ggsave(filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/IDR DESIGN/PDF/rg_vs_re_seq1_800.pdf", 
       device = "pdf", width = 3.5, height = 4, units = "in", dpi = 450)

ggsave(filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/IDR DESIGN/PNG/rg_vs_re_seq1_800.png", 
       device = "png", width = 3.5, height = 4, units = "in", dpi = 400)






# Add directory for simulation at 200 mM of NaCl
dir1000 <- "D:/MASTER_FILES/DATA/sed1_NaCl_1000_IDRLab/sed1_NaCl_1000_IDRLab/time_series_Rg_Ree_Delta_S.csv"

# Load data of CALVADOS at 200 mM
temp_md5 <- read.csv(file = dir1000, header = TRUE)


# Add a new column with the NaCl concentration used
temp_md5$NaCl <- "1000"

# Add temp_md with temp_all
temp_all <- rbind(temp_all, temp_md5)

# Create scatterplot
ggplot(data = temp_md5, aes(x = Rg..nm., y = Ree..nm.)) + 
  geom_point(color  = "#8B00DF", alpha = 0.8) +
  theme_bw() + 
  labs(x = "Rg (nm)", y = "Ree (nm)") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        panel.grid = element_blank())

# Save plots
ggsave(filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/IDR DESIGN/PDF/rg_vs_re_seq1_1000.pdf", 
       device = "pdf", width = 3.5, height = 4, units = "in", dpi = 450)

ggsave(filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/IDR DESIGN/PNG/rg_vs_re_seq1_1000.png", 
       device = "png", width = 3.5, height = 4, units = "in", dpi = 400)





# Add directory for simulation at 200 mM of NaCl
dir1500 <- "D:/MASTER_FILES/DATA/sed1_NaCl_1500_IDRLab/sed1_NaCl_1500_IDRLab/time_series_Rg_Ree_Delta_S.csv"

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

# Save plots
ggsave(filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/IDR DESIGN/PDF/rg_vs_re_seq1_1500.pdf", 
       device = "pdf", width = 3.5, height = 4, units = "in", dpi = 450)

ggsave(filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/IDR DESIGN/PNG/rg_vs_re_seq1_1500.png", 
       device = "png", width = 3.5, height = 4, units = "in", dpi = 400)


# Order data by NaCl concentration
temp_all$NaCl <- factor(temp_all$NaCl, levels = c("200", "400", "600", "800", "1000", "1500"))

# Generate violin plot 
ggplot(data = temp_all, aes(x = NaCl, y = Rg..nm., fill = NaCl)) +
  geom_violin(show.legend = FALSE) +
  theme_bw() +
  labs(x = "[NaCl] (mM)", y = "Rg (nm)") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(1.5, 6)) +
  scale_fill_manual(values = c("#e8c2ff", "#d99fff", "#ca7cff", "#bc5aff", "#ad37ff", "#9e14ff"))

# Save plots
ggsave(filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/IDR DESIGN/PDF/rg_violin_seq1.pdf", 
       device = "pdf", width = 4, height = 4, units = "in", dpi = 450)

ggsave(filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/IDR DESIGN/PNG/rg_violin_seq1.png", 
       device = "png", width = 4, height = 4, units = "in", dpi = 400)




# Generate violin plot 
ggplot(data = temp_all, aes(x = NaCl, y = Ree..nm., fill = NaCl)) +
  geom_violin(show.legend = FALSE) +
  theme_bw() +
  labs(x = "[NaCl] (mM)", y = "Ree (nm)") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(0, 20)) + 
  scale_fill_manual(values = c("#86cd6b", "#73b56a", "#5f9d6a", "#4c8669", "#386e69", "#255668"))

# Save plots
ggsave(filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/IDR DESIGN/PDF/ree_violin_seq1.pdf", 
       device = "pdf", width = 4, height = 4, units = "in", dpi = 450)

ggsave(filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/IDR DESIGN/PNG/ree_violin_seq1.png", 
       device = "png", width = 4, height = 4, units = "in", dpi = 400)


# Change name of columns of temp_all
names(temp_all) <- c("Frame", "rg", "re", "delta", "s", "Treatment")

# Load data set from all de novo IDRs with their Rg
df_novo <- read.csv(file = "D:/MASTER_FILES/DATA/all_de_novo_IDR_calvados.csv", header = TRUE)

# Obtain the mean of tge Rg from eahc condition
temp_sed1 <- temp_all %>% group_by(NaCl) %>% summarise(Mean_Rg = mean(Rg..nm.), Mean_Ree = mean(Ree..nm.))

# Add an ID for the sequence
temp_sed1$seq <- "LEA4-5"

# Change order of columns of temp_sed1
temp_sed1 <- temp_sed1[, c(4, 1, 2, 3)]

# Change name of the columns
names(temp_sed1) <- c("ID", "Treatment", "rg", "re")

# Merge data sets
df_novo <- rbind(df_novo, temp_sed1)

# Order ID column as factor
df_novo$Treatment <- factor(df_novo$Treatment, levels = c("0", "200", "400", "600", "800", "1000", "1500"))

# Remove NaCl at 0 mM
df_novo <- df_novo[!df_novo$Treatment == "0", ]

# Plot 
ggplot(data = df_novo, aes(x = Treatment, y = rg, color = ID, group = ID)) +
  geom_point(size = 1.8) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  labs(x = "NaCl (mM)", y = "Rg (nm)") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(1.5, 4.1)) + 
  scale_y_continuous(breaks = seq(1.5, 4.1, 0.5)) +
  scale_color_manual(name = NULL, 
                     values = c("#00468B", "#ED0000", "#42B540", "#0099B4", "#9A0BA9"))

# Save plots
ggsave(filename = "D:/MASTER_FILES/PLOTS/rg_calvados_idr_de_novo.pdf", 
       device = "pdf", width = 4.5, height = 3.5, units = "in", dpi = 450)

ggsave(filename = "D:/MASTER_FILES/PLOTS/rg_calvados_idr_de_novo.png", 
       device = "png", width = 4.5, height = 3.5, units = "in", dpi = 450)



# Plot Rg series for each condition
ggplot(data = temp_all[temp_all$Treatment == 200, ], aes(x = Frame, y = rg)) +
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
ggplot(data = temp_all[temp_all$Treatment == 200, ], aes(x = rg)) +
  geom_density(color = "#0011df", linewidth = 0.6) +
  theme_bw() +
  labs(x = "Rg (nm)", y = "") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank())



