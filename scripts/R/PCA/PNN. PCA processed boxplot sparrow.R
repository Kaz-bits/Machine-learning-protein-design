# Paqueterías
library(ggplot2)
library(ggfortify)
library(cluster)
library(factoextra)
library(parameters)

# Cargar archivo
df_sparrow <- read.csv(file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/FASTA/IDRBS_library_sparrow_186.csv", 
                       header = TRUE)

# Cargar datos de CIDER (186 IDRs)
cider <- read.csv(file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/df_fret_hml_delta.csv", 
                  header = TRUE)

# Cargar biblioteca de biosensores
bios <- readxl::read_xlsx(path = "D:/ME/ALL/Project_IDP_D2P2/DATA/DATA FOR R/DATA FOR ANALYSIS/200 IDRBS Library.xlsx", 
                          sheet = 1, col_names = TRUE)


# Añadir secuencias a los data frames de parrot
temp <- bios[as.numeric(substr(bios$Entry, start = 7, 9)) %in% cider$construct,]

# Ordenar datos por constructo
temp <- temp[order(temp$Entry), ]

# Ordenar datos de cider por constructo
cider <- cider[order(cider$construct), ]

# Añadir secuencias a cider_parrot_1
cider$sequence <- temp$`IDR sequence to order`

# Juntar columnas de cider con df_sparrow
df_sparrow <- cbind(df_sparrow, cider[, c(2, 10, 21)])

# Arreglar columnas
df_sparrow <- df_sparrow[, c(15, 1:14, 16)]

# Guardar archivo 
write.table(x = df_sparrow, file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/FASTA/IDRBS_library_sparrow_186_FRET.csv", 
            quote = FALSE, sep = ",", row.names = FALSE)


# Extraer datos de cluster
# Cambiar nombres de los datos
df_sparrow[df_sparrow$Response == "High", ]$Response <- "Alta"
df_sparrow[df_sparrow$Response == "Medium", ]$Response <- "Media"
df_sparrow[df_sparrow$Response == "Low", ]$Response <- "Baja"

# Ordenar datos
df_sparrow <- df_sparrow[order(df_sparrow$Entry), ]

# Añadir datos del clustering (PAM)
# Escalar los datos de sparrow
sp_sc <- scale(df_sparrow[, c(3:15)])

# Determinar el número de clusters
fviz_nbclust(sp_sc, pam, method = "silhouette")

# Realizar análisis de PAM
pam.res <- pam(sp_sc, 2)

# Añadir cluster a los datos de sparrow
df_sparrow <- cbind(df_sparrow, cluster=pam.res$cluster)


# Gráficos de distribución ----

# Extraer datos de kappa del cluster 1
i <- df_sparrow[df_sparrow$cluster == 1, ]$kappa
h <- df_sparrow[df_sparrow$cluster == 2, ]$kappa

# Realizar prueba de t para los datos
a <- t.test(i, h)[[3]]


# Nivel de significancia
signif_level <- c(ns = 1, "*" = 0.05, "**" = 0.01, 
                  "***" = 0.001, "****" = 0.0001)

# Evaluar el valor de p para determinar significancia
# Para 0 vs 200
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

#Extraer los simbolos de significancia del valor de probabilidad anterior
p_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1

# Gráfico 1. kappa vs FRET
ggplot(data = df_sparrow, aes(x = as.factor(cluster), y = kappa, fill = as.factor(cluster))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = expression(kappa)) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD")) +
  annotate(geom = "text", x = 2, y = max(df_sparrow$kappa) + 0.1, size = 4, label = p_box)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/kappa_vs_fret_box_sparrow.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/kappa_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del cluster 1
i <- df_sparrow[df_sparrow$cluster == 1, ]$FCR
h <- df_sparrow[df_sparrow$cluster == 2, ]$FCR

# Realizar prueba de t para los datos
a <- t.test(i, h)[[3]]


# Nivel de significancia
signif_level <- c(ns = 1, "*" = 0.05, "**" = 0.01, 
                  "***" = 0.001, "****" = 0.0001)

# Evaluar el valor de p para determinar significancia
# Para 0 vs 200
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

#Extraer los simbolos de significancia del valor de probabilidad anterior
p_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1

# Gráfico 2. FCR vs FRET
ggplot(data = df_sparrow, aes(x = as.factor(cluster), y = FCR, fill = as.factor(cluster))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "FCR") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD")) +
  annotate(geom = "text", x = 2, y = max(df_sparrow$FCR) + 0.15, size = 4, label = p_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/FCR_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/FCR_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)



# Extraer datos de kappa del cluster 1
i <- df_sparrow[df_sparrow$cluster == 1, ]$NCPR
h <- df_sparrow[df_sparrow$cluster == 2, ]$NCPR

# Realizar prueba de t para los datos
a <- t.test(i, h)[[3]]


# Nivel de significancia
signif_level <- c(ns = 1, "*" = 0.05, "**" = 0.01, 
                  "***" = 0.001, "****" = 0.0001)

# Evaluar el valor de p para determinar significancia
# Para 0 vs 200
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

#Extraer los simbolos de significancia del valor de probabilidad anterior
p_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1


# Gráfico 3. NCPR vs FRET
ggplot(data = df_sparrow, aes(x = as.factor(cluster), y = NCPR, fill = as.factor(cluster))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "NCPR") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD")) +
  annotate(geom = "text", x = 2, y = max(df_sparrow$NCPR) + 0.15, size = 4, label = p_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/NCPR_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/NCPR_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del cluster 1
i <- df_sparrow[df_sparrow$cluster == 1, ]$hydropathy
h <- df_sparrow[df_sparrow$cluster == 2, ]$hydropathy

# Realizar prueba de t para los datos
a <- t.test(i, h)[[3]]


# Nivel de significancia
signif_level <- c(ns = 1, "*" = 0.05, "**" = 0.01, 
                  "***" = 0.001, "****" = 0.0001)

# Evaluar el valor de p para determinar significancia
# Para 0 vs 200
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

#Extraer los simbolos de significancia del valor de probabilidad anterior
p_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1


# Gráfico 4. NCPR vs FRET
ggplot(data = df_sparrow, aes(x = as.factor(cluster), y = hydropathy, fill = as.factor(cluster))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "Hidropatia") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD")) +
  annotate(geom = "text", x = 2, y = max(df_sparrow$hydropathy) + 0.75, size = 4, label = p_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/Hidropatia_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/Hidropatia_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del cluster 1
i <- df_sparrow[df_sparrow$cluster == 1, ]$SCD
h <- df_sparrow[df_sparrow$cluster == 2, ]$SCD

# Realizar prueba de t para los datos
a <- t.test(i, h)[[3]]


# Nivel de significancia
signif_level <- c(ns = 1, "*" = 0.05, "**" = 0.01, 
                  "***" = 0.001, "****" = 0.0001)

# Evaluar el valor de p para determinar significancia
# Para 0 vs 200
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

#Extraer los simbolos de significancia del valor de probabilidad anterior
p_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1


# Gráfico 4. NCPR vs FRET
ggplot(data = df_sparrow, aes(x = as.factor(cluster), y = SCD, fill = as.factor(cluster))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "SCD") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD")) +
  annotate(geom = "text", x = 2, y = max(df_sparrow$SCD)+ 5, size = 4, label = p_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/SCD_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/SCD_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del cluster 1
i <- df_sparrow[df_sparrow$cluster == 1, ]$SHD
h <- df_sparrow[df_sparrow$cluster == 2, ]$SHD

# Realizar prueba de t para los datos
a <- t.test(i, h)[[3]]


# Nivel de significancia
signif_level <- c(ns = 1, "*" = 0.05, "**" = 0.01, 
                  "***" = 0.001, "****" = 0.0001)

# Evaluar el valor de p para determinar significancia
# Para 0 vs 200
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

#Extraer los simbolos de significancia del valor de probabilidad anterior
p_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1


# Gráfico 4. NCPR vs FRET
ggplot(data = df_sparrow, aes(x = as.factor(cluster), y = SHD, fill = as.factor(cluster))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "SHD") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD")) +
  annotate(geom = "text", x = 2, y = max(df_sparrow$SHD) + 1, size = 4, label = p_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/SHD_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/SHD_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del cluster 1
i <- df_sparrow[df_sparrow$cluster == 1, ]$aromatic
h <- df_sparrow[df_sparrow$cluster == 2, ]$aromatic

# Realizar prueba de t para los datos
a <- t.test(i, h)[[3]]


# Nivel de significancia
signif_level <- c(ns = 1, "*" = 0.05, "**" = 0.01, 
                  "***" = 0.001, "****" = 0.0001)

# Evaluar el valor de p para determinar significancia
# Para 0 vs 200
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

#Extraer los simbolos de significancia del valor de probabilidad anterior
p_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1


# Gráfico 4. NCPR vs FRET
ggplot(data = df_sparrow, aes(x = as.factor(cluster), y = aromatic, fill = as.factor(cluster))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "Fracción de residuos aromáticos") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD")) +
  annotate(geom = "text", x = 2, y = max(df_sparrow$aromatic) + 0.02, size = 4, label = p_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/aromatic_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/aromatic_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del cluster 1
i <- df_sparrow[df_sparrow$cluster == 1, ]$aliphatic
h <- df_sparrow[df_sparrow$cluster == 2, ]$aliphatic

# Realizar prueba de t para los datos
a <- t.test(i, h)[[3]]


# Nivel de significancia
signif_level <- c(ns = 1, "*" = 0.05, "**" = 0.01, 
                  "***" = 0.001, "****" = 0.0001)

# Evaluar el valor de p para determinar significancia
# Para 0 vs 200
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

#Extraer los simbolos de significancia del valor de probabilidad anterior
p_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1


# Gráfico 4. NCPR vs FRET
ggplot(data = df_sparrow, aes(x = as.factor(cluster), y = aliphatic, fill = as.factor(cluster))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "Fracción de residuos alifáticos") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD")) +
  annotate(geom = "text", x = 2, y = max(df_sparrow$aliphatic) + 0.04, size = 4, label = p_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/aliphatic_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/aliphatic_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del cluster 1
i <- df_sparrow[df_sparrow$cluster == 1, ]$polar
h <- df_sparrow[df_sparrow$cluster == 2, ]$polar

# Realizar prueba de t para los datos
a <- t.test(i, h)[[3]]


# Nivel de significancia
signif_level <- c(ns = 1, "*" = 0.05, "**" = 0.01, 
                  "***" = 0.001, "****" = 0.0001)

# Evaluar el valor de p para determinar significancia
# Para 0 vs 200
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

#Extraer los simbolos de significancia del valor de probabilidad anterior
p_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1


# Gráfico 4. NCPR vs FRET
ggplot(data = df_sparrow, aes(x = as.factor(cluster), y = polar, fill = as.factor(cluster))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "Fracción de residuos polares") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD")) +
  annotate(geom = "text", x = 2, y = max(df_sparrow$polar) + 0.04, size = 4, label = p_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/polar_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/polar_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del cluster 1
i <- df_sparrow[df_sparrow$cluster == 1, ]$proline
h <- df_sparrow[df_sparrow$cluster == 2, ]$proline

# Realizar prueba de t para los datos
a <- t.test(i, h)[[3]]


# Nivel de significancia
signif_level <- c(ns = 1, "*" = 0.05, "**" = 0.01, 
                  "***" = 0.001, "****" = 0.0001)

# Evaluar el valor de p para determinar significancia
# Para 0 vs 200
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

#Extraer los simbolos de significancia del valor de probabilidad anterior
p_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1


# Gráfico 4. NCPR vs FRET
ggplot(data = df_sparrow, aes(x = as.factor(cluster), y = proline, fill = as.factor(cluster))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "Fracción de prolina") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD")) +
  annotate(geom = "text", x = 2, y = max(df_sparrow$proline) + 0.04, size = 4, label = p_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/proline_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/proline_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)



# Extraer datos de kappa del cluster 1
i <- df_sparrow[df_sparrow$cluster == 1, ]$complex
h <- df_sparrow[df_sparrow$cluster == 2, ]$complex

# Realizar prueba de t para los datos
a <- t.test(i, h)[[3]]


# Nivel de significancia
signif_level <- c(ns = 1, "*" = 0.05, "**" = 0.01, 
                  "***" = 0.001, "****" = 0.0001)

# Evaluar el valor de p para determinar significancia
# Para 0 vs 200
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

#Extraer los simbolos de significancia del valor de probabilidad anterior
p_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1


# Gráfico 4. NCPR vs FRET
ggplot(data = df_sparrow, aes(x = as.factor(cluster), y = complex, fill = as.factor(cluster))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "Complejidad") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD")) +
  annotate(geom = "text", x = 2, y = max(df_sparrow$complex) + 0.04, size = 4, label = p_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/complex_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/complex_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)





# Extraer datos de kappa del cluster 1
i <- df_sparrow[df_sparrow$cluster == 1, ]$MW
h <- df_sparrow[df_sparrow$cluster == 2, ]$MW

# Realizar prueba de t para los datos
a <- t.test(i, h)[[3]]


# Nivel de significancia
signif_level <- c(ns = 1, "*" = 0.05, "**" = 0.01, 
                  "***" = 0.001, "****" = 0.0001)

# Evaluar el valor de p para determinar significancia
# Para 0 vs 200
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

#Extraer los simbolos de significancia del valor de probabilidad anterior
p_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1


# Gráfico 4. NCPR vs FRET
ggplot(data = df_sparrow, aes(x = as.factor(cluster), y = MW, fill = as.factor(cluster))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "Peso molecular") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD")) +
  annotate(geom = "text", x = 2, y = max(df_sparrow$MW) + 0.04, size = 4, label = p_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/MW_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/MW_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del cluster 1
i <- df_sparrow[df_sparrow$cluster == 1, ]$mean_delta
h <- df_sparrow[df_sparrow$cluster == 2, ]$mean_delta

# Realizar prueba de t para los datos
a <- t.test(i, h)[[3]]


# Nivel de significancia
signif_level <- c(ns = 1, "*" = 0.05, "**" = 0.01, 
                  "***" = 0.001, "****" = 0.0001)

# Evaluar el valor de p para determinar significancia
# Para 0 vs 200
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

#Extraer los simbolos de significancia del valor de probabilidad anterior
p_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1


# Gráfico 4. NCPR vs FRET
ggplot(data = df_sparrow, aes(x = as.factor(cluster), y = mean_delta, fill = as.factor(cluster))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "Peso molecular") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD")) +
  annotate(geom = "text", x = 2, y = max(df_sparrow$mean_delta) + 0.04, size = 4, label = p_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/delta_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/delta_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)

