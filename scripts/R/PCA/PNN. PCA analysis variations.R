# Paqueterías
library(ggplot2)
library(ggfortify)
library(cluster)
library(factoextra)
library(parameters)
library(hopkins)


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
temp_df <- unname(pam(cider[c(2, 15:19)], k = 2)[3])[[1]]

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

# Reemplazar NAs con (No)
cider[is.na(cider)] <- "-"

# Eliminar columnas
cider_f <- cider[-c(4:14, 20)]

# Guardar archivo
write.table(x = cider_f, file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/cider_library_PCA.csv", 
            sep = ",", quote = FALSE, row.names = FALSE)


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
  annotate(geom = "text", y = 2, x = 0.1, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 0.25, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/kappa_vs_fret.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/kappa_vs_fret.png",
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
  annotate(geom = "text", y = 2, x = 0.28, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/FCR_vs_fret.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/FCR_vs_fret.png",
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
  annotate(geom = "text", y = 2, x = -0.17, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/NCPR_vs_fret.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/NCPR_vs_fret.png",
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
  annotate(geom = "text", y = 2, x = 2.10, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 2.80, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/Hidropatia_vs_fret.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/Hidropatia_vs_fret.png",
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
  annotate(geom = "text", y = 2, x = 0.74, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/disorder_vs_fret.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/disorder_vs_fret.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)




# Valor de correlación
a <- round(unname(cor.test(cider$mean_delta, cider$IDR.length, method = "pearson")$estimate), 2)

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
  annotate(geom = "text", y = 2, x = 65, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 105, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/longitud_vs_fret.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/longitud_vs_fret.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)



# Gráficos de PCA ----
# Preparar los datos para el PCA
prcomp(cider[, c(15:19)], scale. = TRUE)

# Realizar análisis de PCA
principal_components(cider[, c(15:19)], n = "auto")


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
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  coord_cartesian(xlim = c(-0.3, 0.3)) +
  scale_color_manual(name = NULL, 
                     values = c("#3288FF", "#65BFFF", "#99E5FF"))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/PCA/PCA_var5.pdf",
       device = "pdf", width = 4.5, height = 3.5, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/PCA/PCA_var5.png",
       device = "png", width = 4.5, height = 3.5, units = "in", dpi = 400)


# Preparar los datos para el PCA
prcomp(cider[, c(2, 15:19)], scale. = TRUE)

# Realizar análisis de PCA
principal_components(cider[, c(2, 15:19)], n = "all")


# Graficar PCA con cinco variables
autoplot(prcomp(cider[, c(2, 15:19)], scale. = TRUE), 
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
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  coord_cartesian(xlim = c(-0.3, 0.3)) +
  scale_color_manual(name = NULL, 
                     values = c("#3288FF", "#65BFFF", "#99E5FF"))


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/PCA/PCA_var5_WTD.pdf",
       device = "pdf", width = 4.5, height = 3.5, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/PCA/PCA_var5_WTD.png",
       device = "png", width = 4.5, height = 3.5, units = "in", dpi = 400)



# Análisis de clusters
# Determinar el valor de grupos (kmeans)
fviz_nbclust(cider[, c(2, 15:19)], kmeans, method = "silhouette")

# Graficar los cluster de k-means
autoplot(kmeans(cider[c(2, 15:19)], 2), data = cider, size = 1.5) +
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
fviz_nbclust(cider[, c(2, 15:19)], pam, method = "silhouette")

# Gráfico de clusters (pam)
autoplot(pam(cider[c(2, 15:19)], k = 2, ), frame = TRUE, frame.type = "norm") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank(), 
        legend.position = "none") +
  scale_color_manual(values = c("#2FA3EE", "#4BCAAD")) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD"))


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/PCA/clustering_pam_group.pdf",
       device = "pdf", width = 6, height = 4.5, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/PCA/clustering_pam_group.png",
       device = "png", width = 6, height = 4.5, units = "in", dpi = 400)


# Ordenar los datos por su respuesta
cider$Response <- factor(cider$Response, levels = c("Alta", "Media", "Baja"))

# Realizar stack barplot
ggplot(data = cider) +
  geom_bar(aes(x = as.factor(clustering), y = ..count../sum(..count..), 
               fill = Response)) +
  theme_classic() +
  labs(x = "Cluster", y = "Cantidad") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) + 
  coord_cartesian(ylim = c(0, 0.65)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.65, 0.1)) +
  scale_fill_manual(name = "Respuesta", 
                     values = c("#3288FF", "#65BFFF", "#99E5FF"))


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/PCA/ALL/barplot_clustering_pam_group.pdf",
       device = "pdf", width = 4, height = 4.5, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/PCA/ALL/barplot_clustering_pam_group.png",
       device = "png", width = 4, height = 4.5, units = "in", dpi = 400)


# Datos de clustering analisis ----
pam(cider[c(2, 15:19)], k = 2)[7]

# Validar que los datos contienen grupos (assesing cluster tendency)
fviz_pca_ind(prcomp(cider[c(2, 15:19)]), habillage = cider$Response,
             geom = "point", 
             palette = "jco", 
             title = "PCA - Datos FRET")

# Generar un clustering
km_cider <- kmeans(x = cider[c(2, 15:19)], centers = 3)
fviz_cluster(object = km_cider, data = cider[c(2, 15:19)])
fviz_dend(x = hclust(dist(cider[c(2, 15:19)])), k = 2, k_colors = "jco")

set.seed(1)
hopkins(X = cider[c(2, 15:19)])
df_dit <- dist(cider[c(2, 15:19)], method = "euclidean")
fviz_dist(dist.obj = df_dit, show_labels = FALSE)


# Gráficos de distribución ----

# Nivel de significancia
signif_level <- c(ns = 1, "*" = 0.05, "**" = 0.01, 
                  "***" = 0.001, "****" = 0.0001)

# Extraer datos de kappa del cluster 1 por respuesta
p1 <- cider[cider$clustering == 1 & cider$Response == "Alta", ]$kappa
p2 <- cider[cider$clustering == 1 & cider$Response == "Media", ]$kappa
p3 <- cider[cider$clustering == 1 & cider$Response == "Baja", ]$kappa

# Realizar prueba de t para los datos
a <- t.test(p1, p2)[[3]]
b <- t.test(p1, p3)[[3]]
c <- t.test(p2, p3)[[3]]

# Juntar valores de p
a <- c(a, b, c)

# Evaluar el valor de p para determinar significancia
# Para alta vs media
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para alta vs baja
j <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para media vs baja
k <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

#Extraer los simbolos de significancia del valor de probabilidad anterior
p1_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p2_box <- (names(signif_level[j]))[length((names(signif_level[j])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p3_box <- (names(signif_level[k]))[length((names(signif_level[k])))] #p1


# Extraer datos de kappa del cluster 2 por respuesta
p4 <- cider[cider$clustering == 2 & cider$Response == "Alta", ]$kappa
p5 <- cider[cider$clustering == 2 & cider$Response == "Media", ]$kappa
p6 <- cider[cider$clustering == 2 & cider$Response == "Baja", ]$kappa

# Realizar prueba de t para los datos
d <- t.test(p4, p5)[[3]]
e <- t.test(p4, p6)[[3]]
f <- t.test(p5, p6)[[3]]

# Juntar valores de p
a <- c(d, e, f)

# Evaluar el valor de p para determinar significancia
# Para alta vs media
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para alta vs baja
j <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para media vs baja
k <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)


#Extraer los simbolos de significancia del valor de probabilidad anterior
p4_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p5_box <- (names(signif_level[j]))[length((names(signif_level[j])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p6_box <- (names(signif_level[k]))[length((names(signif_level[k])))] #p1


# Gráfico 1. kappa vs FRET
ggplot(data = cider, aes(x = as.factor(clustering), y = kappa, 
                         fill = as.factor(Response))) +
  geom_boxplot(width = 0.5) +
  # Comparación alta vs media (c1)
  geom_segment(aes(x = 0.8, xend = 1, y = 0.8, yend = 0.8), linewidth = 0.8) +
  geom_segment(aes(x = 0.8, xend = 0.8, y = 0.8, yend = 0.78), linewidth = 0.8) +
  geom_segment(aes(x = 1, xend = 1, y = 0.8, yend = 0.78), linewidth = 0.8) +
  
  # Comparación alta vs baja (c1)
  geom_segment(aes(x = 0.8, xend = 1.15, y = 0.87, yend = 0.87), linewidth = 0.8) +
  geom_segment(aes(x = 0.8, xend = 0.8, y = 0.87, yend = 0.84), linewidth = 0.8) +
  geom_segment(aes(x = 1.15, xend = 1.15, y = 0.87, yend = 0.84), linewidth = 0.8) +
  
  # Comparación media vs baja (c1)
  geom_segment(aes(x = 0.97, xend = 1.15, y = 0.95, yend = 0.95), linewidth = 0.8) +
  geom_segment(aes(x = 0.97, xend = 0.97, y = 0.95, yend = 0.92), linewidth = 0.8) +
  geom_segment(aes(x = 1.15, xend = 1.15, y = 0.95, yend = 0.92), linewidth = 0.8) +
  
  # Comparación alta vs mediac (c2)
  geom_segment(aes(x = 1.8, xend = 2, y = 0.8, yend = 0.8), linewidth = 0.8) +
  geom_segment(aes(x = 1.8, xend = 1.8, y = 0.8, yend = 0.78), linewidth = 0.8) +
  geom_segment(aes(x = 2, xend = 2, y = 0.8, yend = 0.78), linewidth = 0.8) +
  
  # Comparación alta vs baja (c2)
  geom_segment(aes(x = 1.8, xend = 2.2, y = 0.87, yend = 0.87), linewidth = 0.8) +
  geom_segment(aes(x = 1.8, xend = 1.8, y = 0.87, yend = 0.84), linewidth = 0.8) +
  geom_segment(aes(x = 2.2, xend = 2.2, y = 0.87, yend = 0.84), linewidth = 0.8) +
  
  # Comparación media vs baja (c2)
  geom_segment(aes(x = 2, xend = 2.2, y = 0.95, yend = 0.95), linewidth = 0.8) +
  geom_segment(aes(x = 2, xend = 2, y = 0.95, yend = 0.92), linewidth = 0.8) +
  geom_segment(aes(x = 2.2, xend = 2.2, y = 0.95, yend = 0.92), linewidth = 0.8) +
  
  # Añadir estéticas
  theme_classic() +
  labs(x = "Cluster", y = expression(kappa)) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(name = "Respuesta", 
                    values = c("#3288FF", "#65BFFF", "#99E5FF")) +
  # Datos de cluster 1
  annotate(geom = "text", x = 0.9, y = 0.72 + 0.1, size = 4, 
           label = p1_box) +
  annotate(geom = "text", x = 0.97, y = 0.79 + 0.1, size = 4, 
           label = p2_box) +
  annotate(geom = "text", x = 1.07, y = 0.87 + 0.1, size = 4, 
           label = p3_box) +
  # Datos de cluster 2
  annotate(geom = "text", x = 1.9, y = 0.72 + 0.1, size = 4, 
           label = p4_box) +
  annotate(geom = "text", x = 1.99, y = 0.79 + 0.1, size = 4, 
           label = p5_box) +
  annotate(geom = "text", x = 2.10, y = 0.87 + 0.1, size = 4, 
           label = p6_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/kappa_vs_fret_r.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/kappa_vs_fret_r.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# Nivel de significancia
signif_level <- c(ns = 1, "*" = 0.05, "**" = 0.01, 
                  "***" = 0.001, "****" = 0.0001)

# Extraer datos de kappa del cluster 1 por respuesta
p1 <- cider[cider$clustering == 1 & cider$Response == "Alta", ]$FCR
p2 <- cider[cider$clustering == 1 & cider$Response == "Media", ]$FCR
p3 <- cider[cider$clustering == 1 & cider$Response == "Baja", ]$FCR

# Realizar prueba de t para los datos
a <- t.test(p1, p2)[[3]]
b <- t.test(p1, p3)[[3]]
c <- t.test(p2, p3)[[3]]

# Juntar valores de p
a <- c(a, b, c)

# Evaluar el valor de p para determinar significancia
# Para alta vs media
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para alta vs baja
j <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para media vs baja
k <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

#Extraer los simbolos de significancia del valor de probabilidad anterior
p1_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p2_box <- (names(signif_level[j]))[length((names(signif_level[j])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p3_box <- (names(signif_level[k]))[length((names(signif_level[k])))] #p1


# Extraer datos de kappa del cluster 2 por respuesta
p4 <- cider[cider$clustering == 2 & cider$Response == "Alta", ]$FCR
p5 <- cider[cider$clustering == 2 & cider$Response == "Media", ]$FCR
p6 <- cider[cider$clustering == 2 & cider$Response == "Baja", ]$FCR

# Realizar prueba de t para los datos
d <- t.test(p4, p5)[[3]]
e <- t.test(p4, p6)[[3]]
f <- t.test(p5, p6)[[3]]

# Juntar valores de p
a <- c(d, e, f)

# Evaluar el valor de p para determinar significancia
# Para alta vs media
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para alta vs baja
j <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para media vs baja
k <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)


#Extraer los simbolos de significancia del valor de probabilidad anterior
p4_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p5_box <- (names(signif_level[j]))[length((names(signif_level[j])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p6_box <- (names(signif_level[k]))[length((names(signif_level[k])))] #p1

# Gráfico 2. FCR vs FRET
ggplot(data = cider, aes(x = as.factor(clustering), y = FCR, 
                         fill = as.factor(Response))) +
  geom_boxplot(width = 0.5) +
  
  # Comparación alta vs media (c1)
  geom_segment(aes(x = 0.8, xend = 1, y = 0.8, yend = 0.8), linewidth = 0.8) +
  geom_segment(aes(x = 0.8, xend = 0.8, y = 0.8, yend = 0.78), linewidth = 0.8) +
  geom_segment(aes(x = 1, xend = 1, y = 0.8, yend = 0.78), linewidth = 0.8) +
  
  # Comparación alta vs baja (c1)
  geom_segment(aes(x = 0.8, xend = 1.15, y = 0.87, yend = 0.87), linewidth = 0.8) +
  geom_segment(aes(x = 0.8, xend = 0.8, y = 0.87, yend = 0.84), linewidth = 0.8) +
  geom_segment(aes(x = 1.15, xend = 1.15, y = 0.87, yend = 0.84), linewidth = 0.8) +
  
  # Comparación media vs baja (c1)
  geom_segment(aes(x = 0.97, xend = 1.15, y = 0.95, yend = 0.95), linewidth = 0.8) +
  geom_segment(aes(x = 0.97, xend = 0.97, y = 0.95, yend = 0.92), linewidth = 0.8) +
  geom_segment(aes(x = 1.15, xend = 1.15, y = 0.95, yend = 0.92), linewidth = 0.8) +
  
  # Comparación alta vs mediac (c2)
  geom_segment(aes(x = 1.8, xend = 2, y = 0.8, yend = 0.8), linewidth = 0.8) +
  geom_segment(aes(x = 1.8, xend = 1.8, y = 0.8, yend = 0.78), linewidth = 0.8) +
  geom_segment(aes(x = 2, xend = 2, y = 0.8, yend = 0.78), linewidth = 0.8) +
  
  # Comparación alta vs baja (c2)
  geom_segment(aes(x = 1.8, xend = 2.2, y = 0.87, yend = 0.87), linewidth = 0.8) +
  geom_segment(aes(x = 1.8, xend = 1.8, y = 0.87, yend = 0.84), linewidth = 0.8) +
  geom_segment(aes(x = 2.2, xend = 2.2, y = 0.87, yend = 0.84), linewidth = 0.8) +
  
  # Comparación media vs baja (c2)
  geom_segment(aes(x = 2, xend = 2.2, y = 0.95, yend = 0.95), linewidth = 0.8) +
  geom_segment(aes(x = 2, xend = 2, y = 0.95, yend = 0.92), linewidth = 0.8) +
  geom_segment(aes(x = 2.2, xend = 2.2, y = 0.95, yend = 0.92), linewidth = 0.8) +
  
  # Añadir estéticas
  theme_classic() +
  labs(x = "Cluster", y = "FCR") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(name = "Respuesta", 
                    values = c("#3288FF", "#65BFFF", "#99E5FF")) +
  # Datos de cluster 1
  annotate(geom = "text", x = 0.9, y = 0.72 + 0.1, size = 4, 
           label = p1_box) +
  annotate(geom = "text", x = 0.97, y = 0.79 + 0.1, size = 4, 
           label = p2_box) +
  annotate(geom = "text", x = 1.07, y = 0.87 + 0.1, size = 4, 
           label = p3_box) +
  # Datos de cluster 2
  annotate(geom = "text", x = 1.9, y = 0.72 + 0.1, size = 4, 
           label = p4_box) +
  annotate(geom = "text", x = 1.99, y = 0.79 + 0.1, size = 4, 
           label = p5_box) +
  annotate(geom = "text", x = 2.10, y = 0.87 + 0.1, size = 4, 
           label = p6_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/FCR_vs_fret_r.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/FCR_vs_fret_r.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)




# Nivel de significancia
signif_level <- c(ns = 1, "*" = 0.05, "**" = 0.01, 
                  "***" = 0.001, "****" = 0.0001)

# Extraer datos de kappa del cluster 1 por respuesta
p1 <- cider[cider$clustering == 1 & cider$Response == "Alta", ]$NCPR
p2 <- cider[cider$clustering == 1 & cider$Response == "Media", ]$NCPR
p3 <- cider[cider$clustering == 1 & cider$Response == "Baja", ]$NCPR

# Realizar prueba de t para los datos
a <- t.test(p1, p2)[[3]]
b <- t.test(p1, p3)[[3]]
c <- t.test(p2, p3)[[3]]

# Juntar valores de p
a <- c(a, b, c)

# Evaluar el valor de p para determinar significancia
# Para alta vs media
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para alta vs baja
j <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para media vs baja
k <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

#Extraer los simbolos de significancia del valor de probabilidad anterior
p1_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p2_box <- (names(signif_level[j]))[length((names(signif_level[j])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p3_box <- (names(signif_level[k]))[length((names(signif_level[k])))] #p1


# Extraer datos de kappa del cluster 2 por respuesta
p4 <- cider[cider$clustering == 2 & cider$Response == "Alta", ]$NCPR
p5 <- cider[cider$clustering == 2 & cider$Response == "Media", ]$NCPR
p6 <- cider[cider$clustering == 2 & cider$Response == "Baja", ]$NCPR

# Realizar prueba de t para los datos
d <- t.test(p4, p5)[[3]]
e <- t.test(p4, p6)[[3]]
f <- t.test(p5, p6)[[3]]

# Juntar valores de p
a <- c(d, e, f)

# Evaluar el valor de p para determinar significancia
# Para alta vs media
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para alta vs baja
j <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para media vs baja
k <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)


#Extraer los simbolos de significancia del valor de probabilidad anterior
p4_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p5_box <- (names(signif_level[j]))[length((names(signif_level[j])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p6_box <- (names(signif_level[k]))[length((names(signif_level[k])))] #p1

# Gráfico 3. NCPR vs FRET
ggplot(data = cider, aes(x = as.factor(clustering), y = NCPR, 
                         fill = as.factor(Response))) +
  geom_boxplot(width = 0.5) +
  
  #Comparación alta vs media (c1)
  geom_segment(aes(x = 0.8, xend = 1, y = 0.8, yend = 0.8), linewidth = 0.8) +
  geom_segment(aes(x = 0.8, xend = 0.8, y = 0.8, yend = 0.78), linewidth = 0.8) +
  geom_segment(aes(x = 1, xend = 1, y = 0.8, yend = 0.78), linewidth = 0.8) +
  
  # Comparación alta vs baja (c1)
  geom_segment(aes(x = 0.8, xend = 1.15, y = 0.87, yend = 0.87), linewidth = 0.8) +
  geom_segment(aes(x = 0.8, xend = 0.8, y = 0.87, yend = 0.84), linewidth = 0.8) +
  geom_segment(aes(x = 1.15, xend = 1.15, y = 0.87, yend = 0.84), linewidth = 0.8) +
  
  # Comparación media vs baja (c1)
  geom_segment(aes(x = 0.97, xend = 1.15, y = 0.95, yend = 0.95), linewidth = 0.8) +
  geom_segment(aes(x = 0.97, xend = 0.97, y = 0.95, yend = 0.92), linewidth = 0.8) +
  geom_segment(aes(x = 1.15, xend = 1.15, y = 0.95, yend = 0.92), linewidth = 0.8) +
  
  # Comparación alta vs mediac (c2)
  geom_segment(aes(x = 1.8, xend = 2, y = 0.8, yend = 0.8), linewidth = 0.8) +
  geom_segment(aes(x = 1.8, xend = 1.8, y = 0.8, yend = 0.78), linewidth = 0.8) +
  geom_segment(aes(x = 2, xend = 2, y = 0.8, yend = 0.78), linewidth = 0.8) +
  
  # Comparación alta vs baja (c2)
  geom_segment(aes(x = 1.8, xend = 2.2, y = 0.87, yend = 0.87), linewidth = 0.8) +
  geom_segment(aes(x = 1.8, xend = 1.8, y = 0.87, yend = 0.84), linewidth = 0.8) +
  geom_segment(aes(x = 2.2, xend = 2.2, y = 0.87, yend = 0.84), linewidth = 0.8) +
  
  # Comparación media vs baja (c2)
  geom_segment(aes(x = 2, xend = 2.2, y = 0.95, yend = 0.95), linewidth = 0.8) +
  geom_segment(aes(x = 2, xend = 2, y = 0.95, yend = 0.92), linewidth = 0.8) +
  geom_segment(aes(x = 2.2, xend = 2.2, y = 0.95, yend = 0.92), linewidth = 0.8) +
  
  # Añadir estéticas
  theme_classic() +
  labs(x = "Cluster", y = "NCPR") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(name = "Respuesta", 
                    values = c("#3288FF", "#65BFFF", "#99E5FF")) +
  # Datos de cluster 1
  annotate(geom = "text", x = 0.9, y = 0.72 + 0.1, size = 4, 
           label = p1_box) +
  annotate(geom = "text", x = 0.97, y = 0.79 + 0.1, size = 4, 
           label = p2_box) +
  annotate(geom = "text", x = 1.07, y = 0.87 + 0.1, size = 4, 
           label = p3_box) +
  # Datos de cluster 2
  annotate(geom = "text", x = 1.9, y = 0.72 + 0.1, size = 4, 
           label = p4_box) +
  annotate(geom = "text", x = 1.99, y = 0.79 + 0.1, size = 4, 
           label = p5_box) +
  annotate(geom = "text", x = 2.10, y = 0.87 + 0.1, size = 4, 
           label = p6_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/NCPR_vs_fret_r.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/NCPR_vs_fret_r.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# Nivel de significancia
signif_level <- c(ns = 1, "*" = 0.05, "**" = 0.01, 
                  "***" = 0.001, "****" = 0.0001)

# Extraer datos de kappa del cluster 1 por respuesta
p1 <- cider[cider$clustering == 1 & cider$Response == "Alta", ]$Hidropatia
p2 <- cider[cider$clustering == 1 & cider$Response == "Media", ]$Hidropatia
p3 <- cider[cider$clustering == 1 & cider$Response == "Baja", ]$Hidropatia

# Realizar prueba de t para los datos
a <- t.test(p1, p2)[[3]]
b <- t.test(p1, p3)[[3]]
c <- t.test(p2, p3)[[3]]

# Juntar valores de p
a <- c(a, b, c)

# Evaluar el valor de p para determinar significancia
# Para alta vs media
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para alta vs baja
j <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para media vs baja
k <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

#Extraer los simbolos de significancia del valor de probabilidad anterior
p1_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p2_box <- (names(signif_level[j]))[length((names(signif_level[j])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p3_box <- (names(signif_level[k]))[length((names(signif_level[k])))] #p1


# Extraer datos de kappa del cluster 2 por respuesta
p4 <- cider[cider$clustering == 2 & cider$Response == "Alta", ]$Hidropatia
p5 <- cider[cider$clustering == 2 & cider$Response == "Media", ]$Hidropatia
p6 <- cider[cider$clustering == 2 & cider$Response == "Baja", ]$Hidropatia

# Realizar prueba de t para los datos
d <- t.test(p4, p5)[[3]]
e <- t.test(p4, p6)[[3]]
f <- t.test(p5, p6)[[3]]

# Juntar valores de p
a <- c(d, e, f)

# Evaluar el valor de p para determinar significancia
# Para alta vs media
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para alta vs baja
j <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para media vs baja
k <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)


#Extraer los simbolos de significancia del valor de probabilidad anterior
p4_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p5_box <- (names(signif_level[j]))[length((names(signif_level[j])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p6_box <- (names(signif_level[k]))[length((names(signif_level[k])))] #p1


# Gráfico 4. NCPR vs FRET
ggplot(data = cider, aes(x = as.factor(clustering), y = Hidropatia, 
                         fill = as.factor(Response))) +
  geom_boxplot(width = 0.5) +
  
  #Comparación alta vs media (c1)
  geom_segment(aes(x = 0.8, xend = 1.0, y = 5, yend = 5.0), linewidth = 0.8) +
  geom_segment(aes(x = 0.8, xend = 0.8, y = 5, yend = 4.9), linewidth = 0.8) +
  geom_segment(aes(x = 1.0, xend = 1.0, y = 5, yend = 4.9), linewidth = 0.8) +
  
  # Comparación alta vs baja (c1)
  geom_segment(aes(x = 0.80, xend = 1.15, y = 5.3, yend = 5.3), linewidth = 0.8) +
  geom_segment(aes(x = 0.80, xend = 0.80, y = 5.3, yend = 5.2), linewidth = 0.8) +
  geom_segment(aes(x = 1.15, xend = 1.15, y = 5.3, yend = 5.2), linewidth = 0.8) +
  
  # Comparación media vs baja (c1)
  geom_segment(aes(x = 0.97, xend = 1.15, y = 5.6, yend = 5.6), linewidth = 0.8) +
  geom_segment(aes(x = 0.97, xend = 0.97, y = 5.6, yend = 5.5), linewidth = 0.8) +
  geom_segment(aes(x = 1.15, xend = 1.15, y = 5.6, yend = 5.5), linewidth = 0.8) +
  
  # Comparación alta vs mediac (c2)
  geom_segment(aes(x = 1.8, xend = 2.0, y = 5.0, yend = 5.0), linewidth = 0.8) +
  geom_segment(aes(x = 1.8, xend = 1.8, y = 5.0, yend = 4.9), linewidth = 0.8) +
  geom_segment(aes(x = 2.0, xend = 2.0, y = 5.0, yend = 4.9), linewidth = 0.8) +
  
  # Comparación alta vs baja (c2)
  geom_segment(aes(x = 1.8, xend = 2.2, y = 5.3, yend = 5.3), linewidth = 0.8) +
  geom_segment(aes(x = 1.8, xend = 1.8, y = 5.3, yend = 5.2), linewidth = 0.8) +
  geom_segment(aes(x = 2.2, xend = 2.2, y = 5.3, yend = 5.2), linewidth = 0.8) +
  
  # Comparación media vs baja (c2)
  geom_segment(aes(x = 2.0, xend = 2.2, y = 5.6, yend = 5.6), linewidth = 0.8) +
  geom_segment(aes(x = 2.0, xend = 2.0, y = 5.6, yend = 5.5), linewidth = 0.8) +
  geom_segment(aes(x = 2.2, xend = 2.2, y = 5.6, yend = 5.5), linewidth = 0.8) +
  
  # Añadir estéticas
  theme_classic() +
  labs(x = "Cluster", y = "Hidropatia") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(name = "Respuesta", 
                    values = c("#3288FF", "#65BFFF", "#99E5FF")) +
  # Datos de cluster 1
  annotate(geom = "text", x = 0.9, y = 5 + 0.1, size = 4, 
           label = p1_box) +
  annotate(geom = "text", x = 0.97, y = 5.3 + 0.1, size = 4, 
           label = p2_box) +
  annotate(geom = "text", x = 1.07, y = 5.6 + 0.1, size = 4, 
           label = p3_box) +
  # Datos de cluster 2
  annotate(geom = "text", x = 1.9, y = 5 + 0.1, size = 4, 
           label = p4_box) +
  annotate(geom = "text", x = 1.99, y = 5.3 + 0.1, size = 4, 
           label = p5_box) +
  annotate(geom = "text", x = 2.10, y = 5.6 + 0.1, size = 4, 
           label = p6_box)



# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/Hidropatia_vs_fret_r.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/Hidropatia_vs_fret_r.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)






# Nivel de significancia
signif_level <- c(ns = 1, "*" = 0.05, "**" = 0.01, 
                  "***" = 0.001, "****" = 0.0001)

# Extraer datos de kappa del cluster 1 por respuesta
p1 <- cider[cider$clustering == 1 & cider$Response == "Alta", ]$Desorden
p2 <- cider[cider$clustering == 1 & cider$Response == "Media", ]$Desorden
p3 <- cider[cider$clustering == 1 & cider$Response == "Baja", ]$Desorden

# Realizar prueba de t para los datos
a <- t.test(p1, p2)[[3]]
b <- t.test(p1, p3)[[3]]
c <- t.test(p2, p3)[[3]]

# Juntar valores de p
a <- c(a, b, c)

# Evaluar el valor de p para determinar significancia
# Para alta vs media
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para alta vs baja
j <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para media vs baja
k <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

#Extraer los simbolos de significancia del valor de probabilidad anterior
p1_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p2_box <- (names(signif_level[j]))[length((names(signif_level[j])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p3_box <- (names(signif_level[k]))[length((names(signif_level[k])))] #p1


# Extraer datos de kappa del cluster 2 por respuesta
p4 <- cider[cider$clustering == 2 & cider$Response == "Alta", ]$Desorden
p5 <- cider[cider$clustering == 2 & cider$Response == "Media", ]$Desorden
p6 <- cider[cider$clustering == 2 & cider$Response == "Baja", ]$Desorden

# Realizar prueba de t para los datos
d <- t.test(p4, p5)[[3]]
e <- t.test(p4, p6)[[3]]
f <- t.test(p5, p6)[[3]]

# Juntar valores de p
a <- c(d, e, f)

# Evaluar el valor de p para determinar significancia
# Para alta vs media
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para alta vs baja
j <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para media vs baja
k <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)


#Extraer los simbolos de significancia del valor de probabilidad anterior
p4_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p5_box <- (names(signif_level[j]))[length((names(signif_level[j])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p6_box <- (names(signif_level[k]))[length((names(signif_level[k])))] #p1


# Gráfico 5. NCPR vs FRET
ggplot(data = cider, aes(x = as.factor(clustering), y = Desorden, 
                         fill = as.factor(Response))) +
  geom_boxplot(width = 0.5) +
  
  #Comparación alta vs media (c1)
  geom_segment(aes(x = 0.8, xend = 1.0, y = 1.3, yend = 1.3), linewidth = 0.8) +
  geom_segment(aes(x = 0.8, xend = 0.8, y = 1.3, yend = 1.28), linewidth = 0.8) +
  geom_segment(aes(x = 1.0, xend = 1.0, y = 1.3, yend = 1.28), linewidth = 0.8) +
  
  # Comparación alta vs baja (c1)
  geom_segment(aes(x = 0.80, xend = 1.15, y = 1.4, yend = 1.4), linewidth = 0.8) +
  geom_segment(aes(x = 0.80, xend = 0.80, y = 1.4, yend = 1.38), linewidth = 0.8) +
  geom_segment(aes(x = 1.15, xend = 1.15, y = 1.4, yend = 1.38), linewidth = 0.8) +
  
  # Comparación media vs baja (c1)
  geom_segment(aes(x = 0.97, xend = 1.15, y = 1.5, yend = 1.5), linewidth = 0.8) +
  geom_segment(aes(x = 0.97, xend = 0.97, y = 1.5, yend = 1.48), linewidth = 0.8) +
  geom_segment(aes(x = 1.15, xend = 1.15, y = 1.5, yend = 1.48), linewidth = 0.8) +
  
  # Comparación alta vs mediac (c2)
  geom_segment(aes(x = 1.8, xend = 2.0, y = 1.3, yend = 1.3), linewidth = 0.8) +
  geom_segment(aes(x = 1.8, xend = 1.8, y = 1.3, yend = 1.28), linewidth = 0.8) +
  geom_segment(aes(x = 2.0, xend = 2.0, y = 1.3, yend = 1.28), linewidth = 0.8) +
  
  # Comparación alta vs baja (c2)
  geom_segment(aes(x = 1.8, xend = 2.2, y = 1.4, yend = 1.4), linewidth = 0.8) +
  geom_segment(aes(x = 1.8, xend = 1.8, y = 1.4, yend = 1.38), linewidth = 0.8) +
  geom_segment(aes(x = 2.2, xend = 2.2, y = 1.4, yend = 1.38), linewidth = 0.8) +
  
  # Comparación media vs baja (c2)
  geom_segment(aes(x = 2.0, xend = 2.2, y = 1.5, yend = 1.5), linewidth = 0.8) +
  geom_segment(aes(x = 2.0, xend = 2.0, y = 1.5, yend = 1.48), linewidth = 0.8) +
  geom_segment(aes(x = 2.2, xend = 2.2, y = 1.5, yend = 1.48), linewidth = 0.8) +
  
  # Añadir estéticas
  theme_classic() +
  labs(x = "Cluster", y = "Desorden promovido") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(name = "Respuesta", 
                    values = c("#3288FF", "#65BFFF", "#99E5FF")) +
  # Datos de cluster 1
  annotate(geom = "text", x = 0.9, y = 1.3 + 0.01, size = 4, 
           label = p1_box) +
  annotate(geom = "text", x = 0.97, y = 1.4 + 0.01, size = 4, 
           label = p2_box) +
  annotate(geom = "text", x = 1.07, y = 1.5 + 0.01, size = 4, 
           label = p3_box) +
  # Datos de cluster 2
  annotate(geom = "text", x = 1.9, y = 1.3 + 0.02, size = 4, 
           label = p4_box) +
  annotate(geom = "text", x = 1.99, y = 1.4 + 0.02, size = 4, 
           label = p5_box) +
  annotate(geom = "text", x = 2.10, y = 1.5 + 0.02, size = 4, 
           label = p6_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/disorder_vs_fret_r.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/disorder_vs_fret_r.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# Nivel de significancia
signif_level <- c(ns = 1, "*" = 0.05, "**" = 0.01, 
                  "***" = 0.001, "****" = 0.0001)

# Extraer datos de kappa del cluster 1 por respuesta
p1 <- cider[cider$clustering == 1 & cider$Response == "Alta", ]$IDR.length
p2 <- cider[cider$clustering == 1 & cider$Response == "Media", ]$IDR.length
p3 <- cider[cider$clustering == 1 & cider$Response == "Baja", ]$IDR.length

# Realizar prueba de t para los datos
a <- t.test(p1, p2)[[3]]
b <- t.test(p1, p3)[[3]]
c <- t.test(p2, p3)[[3]]

# Juntar valores de p
a <- c(a, b, c)

# Evaluar el valor de p para determinar significancia
# Para alta vs media
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para alta vs baja
j <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para media vs baja
k <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

#Extraer los simbolos de significancia del valor de probabilidad anterior
p1_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p2_box <- (names(signif_level[j]))[length((names(signif_level[j])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p3_box <- (names(signif_level[k]))[length((names(signif_level[k])))] #p1


# Extraer datos de kappa del cluster 2 por respuesta
p4 <- cider[cider$clustering == 2 & cider$Response == "Alta", ]$IDR.length
p5 <- cider[cider$clustering == 2 & cider$Response == "Media", ]$IDR.length
p6 <- cider[cider$clustering == 2 & cider$Response == "Baja", ]$IDR.length

# Realizar prueba de t para los datos
d <- t.test(p4, p5)[[3]]
e <- t.test(p4, p6)[[3]]
f <- t.test(p5, p6)[[3]]

# Juntar valores de p
a <- c(d, e, f)

# Evaluar el valor de p para determinar significancia
# Para alta vs media
i <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para alta vs baja
j <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)

# Para media vs baja
k <- c(a[1] > 0.05, a[1] <= 0.05, a[1] <= 0.01, a[1] <= 0.001, a[1] <= 0.0001)


#Extraer los simbolos de significancia del valor de probabilidad anterior
p4_box <- (names(signif_level[i]))[length((names(signif_level[i])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p5_box <- (names(signif_level[j]))[length((names(signif_level[j])))] #p1

#Extraer los simbolos de significancia del valor de probabilidad anterior
p6_box <- (names(signif_level[k]))[length((names(signif_level[k])))] #p1


# Gráfico 6. Longitud vs FRET
ggplot(data = cider, aes(x = as.factor(clustering), y = IDR.length, 
                         fill = as.factor(Response))) +
  geom_boxplot(width = 0.5) +
  
  #Comparación alta vs media (c1)
  geom_segment(aes(x = 0.8, xend = 1.0, y = 215, yend = 215), linewidth = 0.8) +
  geom_segment(aes(x = 0.8, xend = 0.8, y = 215, yend = 213), linewidth = 0.8) +
  geom_segment(aes(x = 1.0, xend = 1.0, y = 215, yend = 213), linewidth = 0.8) +
  
  # Comparación alta vs baja (c1)
  geom_segment(aes(x = 0.80, xend = 1.15, y = 225, yend = 225), linewidth = 0.8) +
  geom_segment(aes(x = 0.80, xend = 0.80, y = 225, yend = 223), linewidth = 0.8) +
  geom_segment(aes(x = 1.15, xend = 1.15, y = 225, yend = 223), linewidth = 0.8) +
  
  # Comparación media vs baja (c1)
  geom_segment(aes(x = 0.97, xend = 1.15, y = 235, yend = 235), linewidth = 0.8) +
  geom_segment(aes(x = 0.97, xend = 0.97, y = 235, yend = 233), linewidth = 0.8) +
  geom_segment(aes(x = 1.15, xend = 1.15, y = 235, yend = 233), linewidth = 0.8) +
  
  # Comparación alta vs mediac (c2)
  geom_segment(aes(x = 1.8, xend = 2.0, y = 215, yend = 215), linewidth = 0.8) +
  geom_segment(aes(x = 1.8, xend = 1.8, y = 215, yend = 213), linewidth = 0.8) +
  geom_segment(aes(x = 2.0, xend = 2.0, y = 215, yend = 213), linewidth = 0.8) +
  
  # Comparación alta vs baja (c2)
  geom_segment(aes(x = 1.8, xend = 2.2, y = 225, yend = 225), linewidth = 0.8) +
  geom_segment(aes(x = 1.8, xend = 1.8, y = 225, yend = 223), linewidth = 0.8) +
  geom_segment(aes(x = 2.2, xend = 2.2, y = 225, yend = 223), linewidth = 0.8) +
  
  # Comparación media vs baja (c2)
  geom_segment(aes(x = 2.0, xend = 2.2, y = 235, yend = 235), linewidth = 0.8) +
  geom_segment(aes(x = 2.0, xend = 2.0, y = 235, yend = 233), linewidth = 0.8) +
  geom_segment(aes(x = 2.2, xend = 2.2, y = 235, yend = 233), linewidth = 0.8) +
  
  # Añadir estéticas
  theme_classic() +
  labs(x = "Cluster", y = "Longitud de la IDR (aa)") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(name = "Respuesta", 
                    values = c("#3288FF", "#65BFFF", "#99E5FF")) +
  # Datos de cluster 1
  annotate(geom = "text", x = 0.9, y = 215 + 1, size = 4, 
           label = p1_box) +
  annotate(geom = "text", x = 0.97, y = 225 + 1, size = 4, 
           label = p2_box) +
  annotate(geom = "text", x = 1.07, y = 235 + 1, size = 4, 
           label = p3_box) +
  # Datos de cluster 2
  annotate(geom = "text", x = 1.9, y = 215 + 4, size = 4, 
           label = p4_box) +
  annotate(geom = "text", x = 2.0, y = 225 + 4, size = 4, 
           label = p5_box) +
  annotate(geom = "text", x = 2.10, y = 235 + 4, size = 4, 
           label = p6_box)


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/longitud_vs_fret_r.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/longitud_vs_fret_r.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)

