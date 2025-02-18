# Paqueterías
library(ggplot2)
library(ggfortify)
library(cluster)
library(factoextra)
library(parameters)
library(hopkins)
library(ggpubr)

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
df_sparrow <- df_sparrow[, c(20, 1:19, 21)]

# Guardar archivo 
write.table(x = df_sparrow, file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/FASTA/IDRBS_library_sparrow_186_FRET.csv", 
            quote = FALSE, sep = ",", row.names = FALSE)


# Extraer datos de cluster
# Cambiar nombres de los datos
df_sparrow[df_sparrow$Response == "High", ]$Response <- "Alta"
df_sparrow[df_sparrow$Response == "Medium", ]$Response <- "Media"
df_sparrow[df_sparrow$Response == "Low", ]$Response <- "Baja"

# Ordenar datos
df_sparrow <- df_sparrow[order(df_sparrow$construct), ]

principal_components(df_sparrow[, c(3:17, 20)], n = "all")

# Graficar PCA con cinco variables
autoplot(prcomp(df_sparrow[, c(3:17, 20)], scale. = TRUE), 
         data = df_sparrow, 
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
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/PCA/PCA_sparrow.pdf",
       device = "pdf", width = 4.5, height = 3.5, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/PCA/PCA_sparrow.png",
       device = "png", width = 4.5, height = 3.5, units = "in", dpi = 400)


# Análisis de clusters (pam) ----
# Escalar los datos de sparrow
sp_sc <- scale(df_sparrow[, c(3:17, 20)])

# Determinar el número de clusters
fviz_nbclust(sp_sc, pam, method = "silhouette")

# Realizar análisis de PAM
pam.res <- pam(sp_sc, 2)

# Añadir cluster a los datos de sparrow
df_sparrow <- cbind(df_sparrow, cluster=pam.res$cluster)

# Verificar el cluster (repel)
fviz_cluster(pam.res, palette = c("#2FA3EE", "#4BCAAD"), #Concentration ellipse
             repel = FALSE, #Avoid label overplotting
             ggtheme = theme_classic())


# Gráfico de clusters (pam)
autoplot(pam.res, frame = TRUE, frame.type = "norm") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank(), 
        legend.position = "none") +
  scale_color_manual(values = c("#2FA3EE", "#4BCAAD")) +
  scale_fill_manual(values = c("#2FA3EE", "#4BCAAD"))


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/PCA/cluster_pam_group_sparrow.pdf",
       device = "pdf", width = 6, height = 4.5, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/PCA/cluster_pam_group_sparrow.png",
       device = "png", width = 6, height = 4.5, units = "in", dpi = 400)



# Ordenar los datos por su respuesta
df_sparrow$Response <- factor(df_sparrow$Response, levels = c("Alta", "Media", "Baja"))

# Realizar stack barplot
ggplot(data = df_sparrow) +
  geom_bar(aes(x = as.factor(cluster), y = ..count../sum(..count..), 
               fill = Response)) +
  theme_classic() +
  labs(x = "Cluster", y = "Cantidad") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) + 
  coord_cartesian(ylim = c(0, 0.65)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.65, 0.1)) +
  scale_fill_manual(name = "Respuesta", 
                    values = c("#3288ff", "#54a7ff", "#77c6ff"))


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/PCA/ALL/barplot_cluster_pam_group_sparrow.pdf",
       device = "pdf", width = 4, height = 4.5, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/PCA/ALL/barplot_cluster_pam_group_sparrow.png",
       device = "png", width = 4, height = 4.5, units = "in", dpi = 400)


# Determinar la cantidad de datos en la barras
df_sparrow[df_sparrow$cluster == 1 & df_sparrow$Response == "Alta", ]$Entry
df_sparrow[df_sparrow$cluster == 2 & df_sparrow$Response == "Alta", ]$Entry

df_sparrow[df_sparrow$cluster == 1 & df_sparrow$Response == "Baja", ]$Entry
df_sparrow[df_sparrow$cluster == 2 & df_sparrow$Response == "Baja", ]$Entry

# Extraer los datos de alta del cluster 1
df_sw_alta <- df_sparrow[df_sparrow$cluster == 1 & df_sparrow$Response == "Alta", ] # 38
df_sw_baja <- df_sparrow[df_sparrow$cluster == 2 & df_sparrow$Response == "Baja", ] # 41

# Juntar data frames por renglones
cider <- rbind(df_sw_alta, df_sw_baja)




# Gráficos de distribución ----

# Extraer datos de kappa del Response 1
i <- cider[cider$Response == "Baja", ]$kappa
h <- cider[cider$Response == "Alta", ]$kappa

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
a <- ggplot(data = cider, aes(x = as.factor(Response), y = kappa, fill = as.factor(Response))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Respuesta", y = expression(kappa)) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#3288FF", "#99E5FF")) +
  annotate(geom = "text", x = 2, y = max(cider$kappa) + 0.1, size = 4, label = p_box)

# Guardar gráfico
ggsave(plot = a, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/kappa_vs_fret_box_sparrow.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(plot = a, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/kappa_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del Response 1
i <- cider[cider$Response == "Baja", ]$FCR
h <- cider[cider$Response == "Alta", ]$FCR

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
b <- ggplot(data = cider, aes(x = as.factor(Response), y = FCR, fill = as.factor(Response))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Respuesta", y = "FCR") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#3288FF", "#99E5FF")) +
  annotate(geom = "text", x = 2, y = max(cider$FCR) + 0.15, size = 4, label = p_box)


# Guardar gráfico
ggsave(plot = b, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/FCR_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(plot = b, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/FCR_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)



# Extraer datos de kappa del Response 1
i <- cider[cider$Response == "Baja", ]$NCPR
h <- cider[cider$Response == "Alta", ]$NCPR

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
c <- ggplot(data = cider, aes(x = as.factor(Response), y = NCPR, fill = as.factor(Response))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Respuesta", y = "NCPR") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#3288FF", "#99E5FF")) +
  annotate(geom = "text", x = 2, y = max(cider$NCPR) + 0.15, size = 4, label = p_box)


# Guardar gráfico
ggsave(plot = c, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/NCPR_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(plot = c, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/NCPR_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del Response 1
i <- cider[cider$Response == "Baja", ]$hydropathy
h <- cider[cider$Response == "Alta", ]$hydropathy

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
d <- ggplot(data = cider, aes(x = as.factor(Response), y = hydropathy, fill = as.factor(Response))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Respuesta", y = "Hidrofobicidad") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#3288FF", "#99E5FF")) +
  annotate(geom = "text", x = 2, y = max(cider$hydropathy) + 0.75, size = 4, label = p_box)


# Guardar gráfico
ggsave(plot = d, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/Hidropatia_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(plot = d, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/Hidropatia_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del Response 1
i <- cider[cider$Response == "Baja", ]$SCD
h <- cider[cider$Response == "Alta", ]$SCD

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
e <- ggplot(data = cider, aes(x = as.factor(Response), y = SCD, fill = as.factor(Response))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Respuesta", y = "SCD") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#3288FF", "#99E5FF")) +
  annotate(geom = "text", x = 2, y = max(cider$SCD)+ 5, size = 4, label = p_box)


# Guardar gráfico
ggsave(plot = e, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/SCD_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(plot = e, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/SCD_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del Response 1
i <- cider[cider$Response == "Baja", ]$SHD
h <- cider[cider$Response == "Alta", ]$SHD

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
f <- ggplot(data = cider, aes(x = as.factor(Response), y = SHD, fill = as.factor(Response))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Respuesta", y = "SHD") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#3288FF", "#99E5FF")) +
  annotate(geom = "text", x = 2, y = max(cider$SHD) + 1, size = 4, label = p_box)


# Guardar gráfico
ggsave(plot = f, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/SHD_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(plot = f, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/SHD_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del Response 1
i <- cider[cider$Response == "Baja", ]$aromatic
h <- cider[cider$Response == "Alta", ]$aromatic

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
g <- ggplot(data = cider, aes(x = as.factor(Response), y = aromatic, fill = as.factor(Response))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Respuesta", y = "Fracc. de residuos aromáticos") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#3288FF", "#99E5FF")) +
  annotate(geom = "text", x = 2, y = max(cider$aromatic) + 0.02, size = 4, label = p_box)


# Guardar gráfico
ggsave(plot = g, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/aromatic_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(plot = g, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/aromatic_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del Response 1
i <- cider[cider$Response == "Baja", ]$aliphatic
h <- cider[cider$Response == "Alta", ]$aliphatic

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
h <- ggplot(data = cider, aes(x = as.factor(Response), y = aliphatic, fill = as.factor(Response))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Respuesta", y = "Fracc. de residuos alifáticos") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#3288FF", "#99E5FF")) +
  annotate(geom = "text", x = 2, y = max(cider$aliphatic) + 0.04, size = 4, label = p_box)


# Guardar gráfico
ggsave(plot = h, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/aliphatic_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(plot = h, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/aliphatic_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del Response 1
i <- cider[cider$Response == "Baja", ]$polar
h <- cider[cider$Response == "Alta", ]$polar

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
i <- ggplot(data = cider, aes(x = as.factor(Response), y = polar, fill = as.factor(Response))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Respuesta", y = "Fracc. de residuos polares") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#3288FF", "#99E5FF")) +
  annotate(geom = "text", x = 2, y = max(cider$polar) + 0.04, size = 4, label = p_box)


# Guardar gráfico
ggsave(plot = i, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/polar_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(plot = i, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/polar_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del Response 1
i <- cider[cider$Response == "Baja", ]$proline
h <- cider[cider$Response == "Alta", ]$proline

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
j <- ggplot(data = cider, aes(x = as.factor(Response), y = proline, fill = as.factor(Response))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Respuesta", y = "Fracc. de prolina") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#3288FF", "#99E5FF")) +
  annotate(geom = "text", x = 2, y = max(cider$proline) + 0.04, size = 4, label = p_box)


# Guardar gráfico
ggsave(plot = j, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/proline_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(plot = j, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/proline_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)



# Extraer datos de kappa del Response 1
i <- cider[cider$Response == "Baja", ]$complex
h <- cider[cider$Response == "Alta", ]$complex

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
k <- ggplot(data = cider, aes(x = as.factor(Response), y = complex, fill = as.factor(Response))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Respuesta", y = "Complejidad") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#3288FF", "#99E5FF")) +
  annotate(geom = "text", x = 2, y = max(cider$complex) + 0.04, size = 4, label = p_box)


# Guardar gráfico
ggsave(plot = k, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/complex_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(plot = k, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/complex_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)





# Extraer datos de kappa del Response 1
i <- cider[cider$Response == "Baja", ]$MW
h <- cider[cider$Response == "Alta", ]$MW

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
l <- ggplot(data = cider, aes(x = as.factor(Response), y = MW, fill = as.factor(Response))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Respuesta", y = "Peso molecular") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#3288FF", "#99E5FF")) +
  annotate(geom = "text", x = 2, y = max(cider$MW) + 0.04, size = 4, label = p_box)


# Guardar gráfico
ggsave(plot = l, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/MW_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(plot = l, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/MW_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Extraer datos de kappa del Response 1
i <- cider[cider$Response == "Baja", ]$mean_delta
h <- cider[cider$Response == "Alta", ]$mean_delta

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
m <- ggplot(data = cider, aes(x = as.factor(Response), y = mean_delta, fill = as.factor(Response))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Respuesta", y = expression(Delta*FRET)) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#3288FF", "#99E5FF")) +
  annotate(geom = "text", x = 2, y = max(cider$mean_delta) + 0.04, size = 4, label = p_box)


# Guardar gráfico
ggsave(plot = m, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/delta_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(plot = m, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/delta_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)






# Extraer datos de kappa del Response 1
i <- cider[cider$Response == "Baja", ]$radius_of_gyration
h <- cider[cider$Response == "Alta", ]$radius_of_gyration

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
n <- ggplot(data = cider, aes(x = as.factor(Response), y = radius_of_gyration, fill = as.factor(Response))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Respuesta", y = "Radio de giro (Rg)") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#3288FF", "#99E5FF")) +
  annotate(geom = "text", x = 2, y = max(cider$radius_of_gyration) + 0.04, size = 4, label = p_box)


# Guardar gráfico
ggsave(plot = n, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/rg_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(plot = n, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/rg_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)







# Extraer datos de kappa del Response 1
i <- cider[cider$Response == "Baja", ]$end_to_end_distance
h <- cider[cider$Response == "Alta", ]$end_to_end_distance

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
o <- ggplot(data = cider, aes(x = as.factor(Response), y = end_to_end_distance, fill = as.factor(Response))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Respuesta", y = "Longitud de desplazamiento (Re)") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#3288FF", "#99E5FF")) +
  annotate(geom = "text", x = 2, y = max(cider$end_to_end_distance) + 0.04, size = 4, label = p_box)


# Guardar gráfico
ggsave(plot = o, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/re_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(plot = o, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/re_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)







# Extraer datos de kappa del Response 1
i <- cider[cider$Response == "Baja", ]$asphericity
h <- cider[cider$Response == "Alta", ]$asphericity

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
p <- ggplot(data = cider, aes(x = as.factor(Response), y = asphericity, fill = as.factor(Response))) +
  geom_boxplot(show.legend = FALSE, width = 0.5) +
  theme_classic() +
  labs(x = "Respuesta", y = "Asfericidad") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#3288FF", "#99E5FF")) +
  annotate(geom = "text", x = 2, y = max(cider$asphericity) + 0.04, size = 4, label = p_box)


# Guardar gráfico
ggsave(plot = p, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/asphericity_vs_fret.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(plot = p, filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/CORR/ALL/PCA/SPARROW/GROUPS/asphericity_vs_fret.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)



# Juntar gráficos de boxplot
ggarrange(p, k, m, b, d, a, l, i, o, n, f, 
          ncol = 4, nrow = 4)


# Valores promedio de las variables ----
mean(cider[cider$Response == "Alta", ]$kappa) # 0.15
mean(cider[cider$Response == "Baja", ]$kappa) # 0.25

mean(cider[cider$Response == "Alta", ]$FCR) # 0.23
mean(cider[cider$Response == "Baja", ]$FCR) # 0.37

mean(cider[cider$Response == "Alta", ]$hydropathy) # 3.61
mean(cider[cider$Response == "Baja", ]$hydropathy) # 3.250

mean(cider[cider$Response == "Alta", ]$asphericity) # 0.42
mean(cider[cider$Response == "Baja", ]$asphericity) # 0.43

mean(cider[cider$Response == "Alta", ]$SHD) # 4.60
mean(cider[cider$Response == "Baja", ]$SHD) # 3.38

mean(cider[cider$Response == "Alta", ]$radius_of_gyration) # 0.15
mean(cider[cider$Response == "Baja", ]$radius_of_gyration) # 0.25

mean(cider[cider$Response == "Alta", ]$polar) # 0.42
mean(cider[cider$Response == "Baja", ]$polar) # 0.30

mean(cider[cider$Response == "Alta", ]$complex) # 0.86
mean(cider[cider$Response == "Baja", ]$complex) # 0.81

mean(cider[cider$Response == "Alta", ]$MW) # 16309.21
mean(cider[cider$Response == "Baja", ]$MW) # 8382.66

mean(cider[cider$Response == "Alta", ]$end_to_end_distance) # 90.83
mean(cider[cider$Response == "Baja", ]$end_to_end_distance) # 59.89

mean(cider[cider$Response == "Alta", ]$mean_delta) # 0.74
mean(cider[cider$Response == "Baja", ]$mean_delta) # 0.13




# Rango de valores ----
# Valores máximo y mínimo de alta
# Rango de valores ----
# Valores máximo y mínimo de alta
min(cider[cider$Response == "Alta", ]$kappa) # 0.15
max(cider[cider$Response == "Alta", ]$kappa) # 0.25

min(cider[cider$Response == "Alta", ]$FCR) # 0.23
max(cider[cider$Response == "Alta", ]$FCR) # 0.37

min(cider[cider$Response == "Alta", ]$hydropathy) # 3.61
max(cider[cider$Response == "Alta", ]$hydropathy) # 3.250

min(cider[cider$Response == "Alta", ]$asphericity) # 0.42
max(cider[cider$Response == "Alta", ]$asphericity) # 0.43

min(cider[cider$Response == "Alta", ]$SHD) # 4.60
max(cider[cider$Response == "Alta", ]$SHD) # 3.38

min(cider[cider$Response == "Alta", ]$radius_of_gyration) # 0.15
max(cider[cider$Response == "Alta", ]$radius_of_gyration) # 0.25

min(cider[cider$Response == "Alta", ]$polar) # 0.42
max(cider[cider$Response == "Alta", ]$polar) # 0.30

min(cider[cider$Response == "Alta", ]$complex) # 0.86
max(cider[cider$Response == "Alta", ]$complex) # 0.81

min(cider[cider$Response == "Alta", ]$MW) # 16309.21
max(cider[cider$Response == "Alta", ]$MW) # 8382.66

min(cider[cider$Response == "Alta", ]$end_to_end_distance) # 90.83
max(cider[cider$Response == "Alta", ]$end_to_end_distance) # 59.89

min(cider[cider$Response == "Alta", ]$mean_delta) # 0.74
max(cider[cider$Response == "Alta", ]$mean_delta) # 0.13



# Valores máximo y mínimo de baja
min(cider[cider$Response == "Baja", ]$kappa) # 0.15
max(cider[cider$Response == "Baja", ]$kappa) # 0.25

min(cider[cider$Response == "Baja", ]$FCR) # 0.23
max(cider[cider$Response == "Baja", ]$FCR) # 0.37

min(cider[cider$Response == "Baja", ]$hydropathy) # 3.61
max(cider[cider$Response == "Baja", ]$hydropathy) # 3.250

min(cider[cider$Response == "Baja", ]$asphericity) # 0.42
max(cider[cider$Response == "Baja", ]$asphericity) # 0.43

min(cider[cider$Response == "Baja", ]$SHD) # 4.60
max(cider[cider$Response == "Baja", ]$SHD) # 3.38

min(cider[cider$Response == "Baja", ]$radius_of_gyration) # 0.15
max(cider[cider$Response == "Baja", ]$radius_of_gyration) # 0.25

min(cider[cider$Response == "Baja", ]$polar) # 0.42
max(cider[cider$Response == "Baja", ]$polar) # 0.30

min(cider[cider$Response == "Baja", ]$complex) # 0.86
max(cider[cider$Response == "Baja", ]$complex) # 0.81

min(cider[cider$Response == "Baja", ]$MW) # 16309.21
max(cider[cider$Response == "Baja", ]$MW) # 8382.66

min(cider[cider$Response == "Baja", ]$end_to_end_distance) # 90.83
max(cider[cider$Response == "Baja", ]$end_to_end_distance) # 59.89

min(cider[cider$Response == "Baja", ]$mean_delta) # 0.74
max(cider[cider$Response == "Baja", ]$mean_delta) # 0.13



# Calcular los cuartiles de cada grupo de respuesta
summary(cider[cider$Response == "Alta", ])
summary(cider[cider$Response == "Baja", ])





