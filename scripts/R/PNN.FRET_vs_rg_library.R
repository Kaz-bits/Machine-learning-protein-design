# Packages
library(ggplot2)

# Load data
df_sparrow <- read.csv(file = "D:/MASTER_FILES/DATA/IDRBS_library_sparrow_188.csv", 
                       header = TRUE)

# Remove IDRs with medium response
df_sparrow <- df_sparrow[!df_sparrow$Response == "Medium", ]

# Plot
ggplot(data = df_sparrow, aes(x = FCR, y = kappa, color = Response)) +
  geom_point() +
  scale_color_viridis_d()

