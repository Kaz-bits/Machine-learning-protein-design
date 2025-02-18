# Packages
library(ggplot2)
library(scatterplot3d)
library(magick)
library(rgl)

# Load data
df_sparrow <- read.csv(file = "D:/MASTER_FILES/DATA/IDRBS_library_sparrow_188.csv", 
                       header = TRUE)

# Load data
df_goose <- read.csv(file = "D:/MASTER_FILES/DATA/IDRBS_DN_SVM_goose_kappa_890.csv", 
                     header = TRUE)

# Remove data with medium response
df_sparrow <- df_sparrow[!df_sparrow$Response == "Medium", ]

# Convert rg data to nanometers
df_sparrow$rg <- df_sparrow$rg / 10
df_goose$rg <- df_goose$rg / 10

# Vector fo colors
temp_colors <- c("#FFA500", "#00D0D6")

# Assign each color to each point
temp_colors <- temp_colors[as.numeric(factor(df_sparrow$Response))]

# Plot
scatterplot3d(df_sparrow[, c(6, 7, 8)], 
              angle = 45, pch = 20, 
              color = temp_colors, 
              grid = TRUE, box = TRUE,
              xlab = expression(kappa),
              ylab = "FCR",
              zlab = "Rg (nm)")

#Assign each color to each point
temp_colors1 <- temp_colors[as.numeric(factor(df_goose$Prediction))]
                           
# Plot
scatterplot3d(df_goose[, c(1, 2, 3)], 
              angle = 15, pch = 20, 
              color = temp_colors1, 
              grid = TRUE, box = TRUE,
              xlab = expression(kappa),
              ylab = "FCR",
              zlab = "Rg (nm)")



# Build plot with magick and rgl
# Plot
plot3d(df_sparrow[, 6], df_sparrow[, 7], df_sparrow[, 8], 
       col = temp_colors, type = "s", radius = 0.6, 
       xlab = expression(kappa),
       ylab = "FCR",
       zlab = "Rg (nm)")

# Plot
plot3d(df_goose[, 1], 
       df_goose[, 2], 
       df_goose[, 3], 
       col = temp_colors1, type = "p",
       size = 5, 
       xlab = expression(kappa),
       ylab = "FCR",
       zlab = "Rg (nm)")

# Indicate axis and rotation velocity
play3d(spin3d(axis = c(0, 0, 1), rpm = 5), duration = 100)

# Save plot
movie3d(
  movie = "rg_fcr_kappa_3d_animation", 
  spin3d(axis = c(0, 0, 1), rpm = 5),
  duration = 10, 
  dir = "D:/MASTER_FILES/PLOTS",
  type = "gif", 
  clean = TRUE
)
