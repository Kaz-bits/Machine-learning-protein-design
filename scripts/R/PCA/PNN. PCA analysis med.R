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

# Extraer datos de alta
cider <- cider[cider$Response == "Media", ]


# Gráficos de correlación ----

# Valor de correlación
a <- round(unname(cor.test(cider$mean_delta, cider$kappa, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(cider$mean_delta, cider$kappa, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = cider, aes(x = kappa, y = mean_delta)) +
  geom_point() +
  theme_bw() +
  labs(x = expression(kappa), y = "DxAm/DxDm\nnormalizado") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(-0.1, 2)) +
  annotate(geom = "text", y = 2, x = 0.10, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 0.21, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/MED/kappa_vs_fret.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/MED/kappa_vs_fret.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)



# Valor de correlación
a <- round(unname(cor.test(cider$mean_delta, cider$FCR, method = "pearson")$estimate), 2)

# Valor de p
b <- round(cor.test(cider$mean_delta, cider$FCR, method = "pearson")$p.value, digits = 5)


# Gráfico 1. FCR vs FRET
ggplot(data = cider, aes(x = FCR, y = mean_delta)) +
  geom_point() +
  theme_bw() +
  labs(x = "FCR", y = "DxAm/DxDm\nnormalizado") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(-0.1, 2)) +
  annotate(geom = "text", y = 2, x = 0.08, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 0.25, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/MED/FCR_vs_fret.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/MED/FCR_vs_fret.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)



# Valor de correlación
a <- round(unname(cor.test(cider$mean_delta, cider$NCPR, method = "pearson")$estimate), 2)

# Valor de 
b <- round(cor.test(cider$mean_delta, cider$NCPR, method = "pearson")$p.value, digits = 2)


# Gráfico 1. NCPR vs FRET
ggplot(data = cider, aes(x = NCPR, y = mean_delta)) +
  geom_point() +
  theme_bw() +
  labs(x = "NCPR", y = "DxAm/DxDm\nnormalizado") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(-0.1, 2)) +
  annotate(geom = "text", y = 2, x = -0.33, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = -0.20, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/MED/NCPR_vs_fret.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/MED/NCPR_vs_fret.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)




# Valor de correlación
a <- round(unname(cor.test(cider$mean_delta, cider$Hidropatia, method = "pearson")$estimate), 2)

# Valor de p
b <- round(cor.test(cider$mean_delta, cider$Hidropatia, method = "pearson")$p.value, digits = 3)


# Gráfico 1. hidropatia vs FRET
ggplot(data = cider, aes(x = Hidropatia, y = mean_delta)) +
  geom_point() +
  theme_bw() +
  labs(x = "Hidropatía", y = "DxAm/DxDm\nnormalizado") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(-0.1, 2)) +
  annotate(geom = "text", y = 2, x = 2.35, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 2.90, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/MED/Hidropatia_vs_fret.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/MED/Hidropatia_vs_fret.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)



# Valor de correlación
a <- round(unname(cor.test(cider$mean_delta, cider$Desorden, method = "pearson")$estimate), 2)

# Valor de p
b <- round(cor.test(cider$mean_delta, cider$Desorden, method = "pearson")$p.value, digits = 2)


# Gráfico 1. disorder vs FRET
ggplot(data = cider, aes(x = Desorden, y = mean_delta)) +
  geom_point() +
  theme_bw() +
  labs(x = "Desorden promovido", y = "DxAm/DxDm\nnormalizado") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(-0.1, 2)) +
  annotate(geom = "text", y = 2, x = 0.64, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 0.71, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/MED/disorder_vs_fret.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/MED/disorder_vs_fret.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)




# Valor de correlación
a <- round(unname(cor.test(cider$mean_delta, cider$IDR.length, method = "pearson")$estimate), 3)

# Valor de p
b <- format(cor.test(cider$mean_delta, cider$IDR.length, method = "pearson")$p.value, digits = 2)


# Gráfico 1. longitud vs FRET
ggplot(data = cider, aes(x = IDR.length, y = mean_delta)) +
  geom_point() +
  theme_bw() +
  labs(x = "Longitud de la IDR (aa)", y = "DxAm/DxDm\nnormalizado") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(-0.1, 2)) +
  annotate(geom = "text", y = 2, x = 68, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 108, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/MED/longitud_vs_fret.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/MED/longitud_vs_fret.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)



# Gráficos de PCA ----
# Preparar los datos para el PCA
prcomp(cider[, 15:19], scale. = TRUE)

# Realizar análisis de PCA
principal_components(cider[, 15:19], n = "all")


# Graficar PCA con cinco variables
autoplot(prcomp(cider[, c(15:19)], scale. = TRUE), 
         data = cider, 
         color = "Response", 
         loadings = TRUE, 
         loadings.label = TRUE, 
         loadings.color = "black",
         loadings.label.size = 3,
         loadings.label.color = "red") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none") +
  coord_cartesian(xlim = c(-0.3, 0.3)) +
  scale_color_manual(name = NULL, 
                     values = c("#3288FF", "#65BFFF", "#99E5FF"))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/PCA/MED/PCA_var_5_high.pdf",
       device = "pdf", width = 4.5, height = 3.5, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/PCA/MED/PCA_var_5_high.png",
       device = "png", width = 4.5, height = 3.5, units = "in", dpi = 400)


# Análisis de clusters
# Determinar el valor de grupos (kmeans)
fviz_nbclust(cider[, c(15:19)], kmeans, method = "wss")

# Graficar los cluster de k-means
autoplot(kmeans(cider[15:19], 3), data = cider, size = 1.5) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_color_manual(name = NULL, 
                     values = c("#2FA3EE", "#4BCAAD", "#86C157", 
                                "#D99C3F", "#CE6633", "#A35DD1"), 
                     labels = c("Grupo 1", "Grupo 2", "Grupo 3", 
                                "Grupo 4", "Grupo 5", "Grupo 6"))


# Análisis de clusters (pam)
# Determinar el valor de grupos (pam)
fviz_nbclust(cider[, c(15:19)], pam, method = "silhouette")

# Extraer datos de clustering
temp_df <- unname(pam(cider[15:19], k = 3)[3])[[1]]

# Crear un data frame con los datos
temp_df <- data.frame("construct" = as.numeric(names(temp_df)), 
                      "cluster_pam" = unname(temp_df))

# Ordenar datos
temp_df <- temp_df[order(temp_df$construct), ]

# Añadir datos a data frame de cider
cider$clustering <- temp_df$cluster

# Gráfico de clusters (pam)
autoplot(pam(cider[15:19], k = 2), frame = TRUE, frame.type = "norm") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank(), 
        legend.position = "none") +
  scale_color_manual(values = c("#BCBD22", "#9467BD")) +
  scale_fill_manual(values = c("#BCBD22", "#9467BD"))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/PCA/MED/clustering_pam_group.pdf",
       device = "pdf", width = 6, height = 4.5, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/PCA/MED/clustering_pam_group.png",
       device = "png", width = 6, height = 4.5, units = "in", dpi = 400)
