# Packages
library(ggplot2)
library(dplyr)
library(car)

# Load pLDDT data set from library
df_plddt <- read.csv("D:/MASTER_FILES/DATA/library_plddt_META.csv", header = FALSE)

# load data base of the library with the organism properties
df_library <- read.csv(file = "D:/MASTER_FILES/DATA/IDR_Properties.csv", header = TRUE)

# Load data set from SPARROW (188 IDRs)
df_sparrow <- read.csv(file = "D:/MASTER_FILES/DATA/IDRBS_library_sparrow_188.csv", header = TRUE)

# Order the data sets by construct 
df_library <- df_library[order(as.integer(substr(df_library$entry, 7, 9))), ]
df_sparrow <- df_sparrow[order(df_sparrow$Construct), ]

# Filter data sets to obtain the 188 IDRs from the df_library data set
df_library <- df_library[df_library$entry %in% df_sparrow$Construct, ]

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
                     "very_low" = round(plddt_vlow, 2), 
                     "low" = round(plddt_low, 2), 
                     "high" = round(plddt_high, 2), 
                     "very_high" = round(plddt_vhigh, 2))
  
  # Add to data frame
  df_all_plddt <- rbind(df_all_plddt, temp)
  
}

# Filter data 
df_all_plddt <- df_all_plddt[df_all_plddt$construct %in% as.integer(substr(df_library_prop$Construct, 7, 9)), ]

# Add data to df_sparrow
df_library_all <- cbind(df_library_prop, df_all_plddt[, -1])


# Determine probablity values for high confident values
a <- qnorm(p = 0.95, mean = mean(df_library_all$high), sd(df_library_all$high)) # 50.94 (95) - 63.76 (99)
b <- qnorm(p = 0.95, mean = mean(df_library_all$very_high), sd(df_library_all$very_high)) # 46.69 (95) - 61.46 (99)

# Extract data using the probability values
temp_library <- df_library_all[df_library_all$high >= a | df_library_all$very_high >= b, ]

temp_library[temp_library$Response == "Low", ]$Construct


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
ggsave(filename = "D:/MASTER_FILES/PLOTS/plddt_density_plot_all.pdf", device = "pdf", 
       width = 4, height = 2.5, units = "in", dpi = 450)

# Save plot
ggsave(filename = "D:/MASTER_FILES/PLOTS/plddt_density_plot_all.png", device = "png", 
       width = 4, height = 2.5, units = "in", dpi = 450)









# Physicochemical properties of selected IDRs ----
# Filter data with a high FRET response 
df_sparrow_a <- temp_library[temp_library$Response == "High", ]

# Filter data with a low FRET response
df_sparrow_b <- temp_library[temp_library$Response == "Low", ]

# Transformar datos de hidropatía a una escala menor
df_sparrow_a$hydropathy <- (df_sparrow_a$hydropathy) / 10
df_sparrow_b$hydropathy <- (df_sparrow_b$hydropathy) / 10

# Transformar datos de Rg a una escala menor
df_sparrow_a$rg <- (df_sparrow_a$rg) / 100
df_sparrow_b$rg <- (df_sparrow_b$rg) / 100

# Transformar datos de SCD a una escala menor
df_sparrow_a$SCD <- (df_sparrow_a$SCD) / 10
df_sparrow_b$SCD <- (df_sparrow_b$SCD) / 10


# Construir data frames con los datos de sparrow_a
temp_df1a <- df_sparrow_a[, c(1, 32)]
temp_df2a <- df_sparrow_a[, c(1, 33)]
temp_df3a <- df_sparrow_a[, c(1, 34)]
temp_df4a <- df_sparrow_a[, c(1, 35)]
temp_df5a <- df_sparrow_a[, c(1, 44)]
temp_df6a <- df_sparrow_a[, c(1, 29)]
temp_df7a <- df_sparrow_a[, c(1, 39)]
temp_df8a <- df_sparrow_a[, c(1, 40)]
temp_df9a <- df_sparrow_a[, c(1, 41)]
temp_df10a <- df_sparrow_a[, c(1, 36)]


# Construir data frames con los datos de sparrow_b
temp_df1b <- df_sparrow_b[, c(1, 32)]
temp_df2b <- df_sparrow_b[, c(1, 33)]
temp_df3b <- df_sparrow_b[, c(1, 34)]
temp_df4b <- df_sparrow_b[, c(1, 35)]
temp_df5b <- df_sparrow_b[, c(1, 44)]
temp_df6b <- df_sparrow_b[, c(1, 29)]
temp_df7b <- df_sparrow_b[, c(1, 39)]
temp_df8b <- df_sparrow_b[, c(1, 40)]
temp_df9b <- df_sparrow_b[, c(1, 41)]
temp_df10b <- df_sparrow_b[, c(1, 36)]

# Cambiar nombres de las columnas por uno común
names(temp_df1a)[2] <- "Value"
names(temp_df2a)[2] <- "Value"
names(temp_df3a)[2] <- "Value"
names(temp_df4a)[2] <- "Value"
names(temp_df5a)[2] <- "Value"
names(temp_df6a)[2] <- "Value"
names(temp_df7a)[2] <- "Value"
names(temp_df8a)[2] <- "Value"
names(temp_df9a)[2] <- "Value"
names(temp_df10a)[2] <- "Value"

# Cambiar nombres de las columnas por uno común
names(temp_df1b)[2] <- "Value"
names(temp_df2b)[2] <- "Value"
names(temp_df3b)[2] <- "Value"
names(temp_df4b)[2] <- "Value"
names(temp_df5b)[2] <- "Value"
names(temp_df6b)[2] <- "Value"
names(temp_df7b)[2] <- "Value"
names(temp_df8b)[2] <- "Value"
names(temp_df9b)[2] <- "Value"
names(temp_df10b)[2] <- "Value"

# Agregar columna con variable categórica a cada temp_df
temp_df1a$Parameter <- "Kappa"
temp_df2a$Parameter <- "FCR"
temp_df3a$Parameter <- "NCPR"
temp_df4a$Parameter <- "Hydropathy"
temp_df5a$Parameter <- "Rg"
temp_df6a$Parameter <- "Re"
temp_df7a$Parameter <- "Aliphatic"
temp_df8a$Parameter <- "Polar"
temp_df9a$Parameter <- "Proline"
temp_df10a$Parameter <- "SCD"

# Agregar columna con variable categórica a cada temp_df
temp_df1b$Parameter <- "Kappa"
temp_df2b$Parameter <- "FCR"
temp_df3b$Parameter <- "NCPR"
temp_df4b$Parameter <- "Hydropathy"
temp_df5b$Parameter <- "Rg"
temp_df6b$Parameter <- "Re"
temp_df7b$Parameter <- "Aliphatic"
temp_df8b$Parameter <- "Polar"
temp_df9b$Parameter <- "Proline"
temp_df10b$Parameter <- "SCD"

# Unir data frames
df_sparrow1 <- rbind(temp_df1a, temp_df2a, temp_df3a, temp_df4a, temp_df5a, temp_df6a, temp_df7a, temp_df8a, temp_df9a, temp_df10a)
df_sparrow2 <- rbind(temp_df1b, temp_df2b, temp_df3b, temp_df4b, temp_df5b, temp_df6b, temp_df7b, temp_df8b, temp_df9b, temp_df10b)

# Convertir columna de parámetros en un factor
df_sparrow1$Parameter <- factor(df_sparrow1$Parameter)
df_sparrow2$Parameter <- factor(df_sparrow2$Parameter)

# Añadir columna de clasificación
df_sparrow1$Group <- "High"
df_sparrow2$Group <- "Low"

# Unir data frames df_sparrow 1 y 2
df_sparrow1 <- rbind(df_sparrow1, df_sparrow2)

# Convertir columna Group en factor
df_sparrow1$Group <- factor(df_sparrow1$Group)

# Convertir columna del grupo a factor
df_sparrow1$Group <- factor(df_sparrow1$Group)

# Anova de una vía entre los parámetros vs los grupos de alta y baja
res.aov <- aov(data = df_sparrow1, formula = Value ~ Group * Parameter)

# Verificar la estadística del anova
summary(res.aov)

# Realizar prueba post-hoc de Tukey
TukeyHSD(res.aov)

# Realizar prueba de homogeneidad de varianzas
leveneTest(data = df_sparrow1, Value ~ Group * Parameter)

# Realizar prueba de Welch (no hay homogeneidad de varianzas)
res.aov1 <- oneway.test(data = df_sparrow1, formula = Value ~ Group * Parameter)

# Encontrar diferencias entre los grupos
df_sparrow1 %>% pairwise_t_test(Value ~ Parameter)


# Constuir boxplot con los parámetros asignados
ggplot(data = df_sparrow1, aes(x = Parameter, y = Value, fill = Group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_classic() +
  labs(x = NULL, y = "CIDER value") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = "top") +
  scale_x_discrete(labels = c("Aliphatic","FCR", "Hydro", "Kappa", "NCPR", "Polar", "Proline", "Re", "Rg / 100", "SCD")) +
  scale_fill_manual(values = c("#FFA500", "#00D0D6")) +
  annotate(geom = "text", x = 1, y = 0.42, label = "ns", size = 5) +
  annotate(geom = "text", x = 2, y = 0.65, label = "ns", size = 5) +
  annotate(geom = "text", x = 3, y = 0.50, label = "ns", size = 5) +
  annotate(geom = "text", x = 4, y = 0.46, label = "ns", size = 5) +
  annotate(geom = "text", x = 5, y = 0.34, label = "ns", size = 5) +
  annotate(geom = "text", x = 6, y = 0.68, label = "ns", size = 5) +
  annotate(geom = "text", x = 7, y = 0.25, label = "ns", size = 5) +
  annotate(geom = "text", x = 8, y = 0.57, label = "ns", size = 5) +
  annotate(geom = "text", x = 9, y = 0.54, label = "ns", size = 5) +
  annotate(geom = "text", x = 10, y = 0.64, label = "ns", size = 5)


# Create a vector with the comparisons to extract
temp_vec <- c(
  "Low:Kappa-High:Kappa",
  "Low:FCR-High:FCR",
  "Low:NCPR-High:NCPR",
  "Low:Hydropathy-High:Hydropathy",
  "Low:Rg-High:Rg",
  "Low:Re-High:Re",
  "Low:Aliphatic-High:Aliphatic",
  "Low:Polar-High:Polar",
  "Low:Proline-High:Proline",
  "Low:SCD-High:SCD"
)

# Save all the comparison's name
temp_tukey <- row.names(TukeyHSD(res.aov)[[3]])

# Extract all the p-values
temp_p_values <- unname(TukeyHSD(res.aov)[[3]])[, 4]

# Filter data from Tukey test
temp_p_values[temp_tukey %in% temp_vec]


