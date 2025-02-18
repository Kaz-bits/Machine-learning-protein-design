# Paqueterías
library(ggplot2)
library(ggfortify)
library(cluster)
library(factoextra)
library(parameters)

# Cargar datos de CIDER (186 IDRs)
cider <- read.csv(file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/df_fret_hml_delta.csv", 
                  header = TRUE)

# Cambiar nombres de los datos
cider[cider$Response == "High", ]$Response <- "Alta"
cider[cider$Response == "Medium", ]$Response <- "Media"
cider[cider$Response == "Low", ]$Response <- "Baja"

# Cambiar nombre de columnas específicas
names(cider)[18] <- "Hidropatia"
names(cider)[19] <- "Desorden"

# Ordenar datos
cider <- cider[order(cider$construct), ]

# Extraer datos de clustering
temp_df <- unname(pam(cider[15:19], k = 2)[3])[[1]]

# Crear un data frame con los datos
temp_df <- data.frame("rownames" = as.numeric(names(temp_df)), 
                      "cluster_pam" = unname(temp_df))

# Ordenar datos
temp_df <- temp_df[order(temp_df$rownames), ]

# Añadir renglones a data frame de cider
cider$rownames <- as.numeric(rownames(cider))

# Ordenar datos de cider por renglones
cider <- cider[order(cider$rownames), ]

# Añadir datos de clustering a data frame de cider
cider$clustering <- temp_df$cluster_pam

# Remover columna de renglones
cider$rownames <- NULL


# Gráficos de distribución ----

# Extraer datos de kappa del cluster 1
i <- cider[cider$clustering == 1, ]$kappa
h <- cider[cider$clustering == 2, ]$kappa

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
ggplot(data = cider, aes(x = as.factor(clustering), y = kappa, fill = as.factor(clustering))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = expression(kappa)) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD")) +
  annotate(geom = "text", x = 2, y = max(cider$kappa) + 0.1, size = 4, label = p_box)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/kappa_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/kappa_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del cluster 1
i <- cider[cider$clustering == 1, ]$FCR
h <- cider[cider$clustering == 2, ]$FCR

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
ggplot(data = cider, aes(x = as.factor(clustering), y = FCR, fill = as.factor(clustering))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "FCR") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD")) +
  annotate(geom = "text", x = 2, y = max(cider$FCR) + 0.15, size = 4, label = p_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/FCR_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/FCR_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)



# Extraer datos de kappa del cluster 1
i <- cider[cider$clustering == 1, ]$NCPR
h <- cider[cider$clustering == 2, ]$NCPR

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
ggplot(data = cider, aes(x = as.factor(clustering), y = NCPR, fill = as.factor(clustering))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "NCPR") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD")) +
  annotate(geom = "text", x = 2, y = max(cider$NCPR) + 0.15, size = 4, label = p_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/NCPR_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/NCPR_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del cluster 1
i <- cider[cider$clustering == 1, ]$Hidropatia
h <- cider[cider$clustering == 2, ]$Hidropatia

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
ggplot(data = cider, aes(x = as.factor(clustering), y = Hidropatia, fill = as.factor(clustering))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "Hidropatia") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD")) +
  annotate(geom = "text", x = 2, y = max(cider$Hidropatia) - 0.75, size = 4, label = p_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/Hidropatia_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/Hidropatia_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del cluster 1
i <- cider[cider$clustering == 1, ]$Desorden
h <- cider[cider$clustering == 2, ]$Desorden

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

# Gráfico 5. NCPR vs FRET
ggplot(data = cider, aes(x = as.factor(clustering), y = Desorden, fill = as.factor(clustering))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "Desorden") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD")) +
  annotate(geom = "text", x = 2, y = max(cider$Desorden) + 0.05, size = 4, label = p_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/disorder_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/disorder_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del cluster 1
i <- cider[cider$clustering == 1, ]$IDR.length
h <- cider[cider$clustering == 2, ]$IDR.length

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


# Gráfico 6. Longitud vs FRET
ggplot(data = cider, aes(x = as.factor(clustering), y = IDR.length, fill = as.factor(clustering))) +
  geom_boxplot(show.legend = FALSE) +
  theme_classic() +
  labs(x = "Cluster", y = "Longitud de la IDR (aa)") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD")) +
  annotate(geom = "text", x = 2, y = max(cider$IDR.length) + 15, size = 4, label = p_box)



# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/longitud_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/longitud_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Gráficos de correlación ----
ggplot(data = cider, aes(x = FCR, y = kappa, color = Response)) +
  geom_point() +
  geom_vline(xintercept = 0.2, lty = 2, linewidth = 0.6) +
  geom_hline(yintercept = 0.25, lty = 2, linewidth = 0.6) +
  geom_hline(yintercept = 0.12, lty = 2, linewidth = 0.6) +
  theme_bw() +
  labs(y = expression(kappa)) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_color_manual(values = c("#3288FF", "#99E5FF", "#65BFFF"))


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/longitud_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/longitud_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




















