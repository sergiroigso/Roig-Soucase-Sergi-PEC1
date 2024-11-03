## ACCESO A LOS DATOS

# utilizamos el paquete metabolomicsWorkbenchR
library(metabolomicsWorkbenchR)
# descargamos la clase SummarizedExperiment
# del estudio ST002892
SE = do_query(
  context = 'study',
  input_item = 'study_id',
  input_value = 'ST002892',
  output_item = 'SummarizedExperiment'
)
SE

## FORMATO Y TIPO DE LOS DATOS

# dimensión de la clase SE
dim(SE)
# podemos ver también la metadata de SE
# para ello usamos el paquete SummarizedExperiment
library(SummarizedExperiment)
# función metadata
metadata(SE)
# la información de las filas
rowData(SE)
# también podemos acceder a las columnas
colData(SE)
# podemos acceder a los datos con la función assays
# podríamos acceder a un solo ensayo con la función assay
assays(SE)[[1]]

## REPRESENTACIÓN DE LOS DATOS

# paquetes necesarios
library(ggplot2)
library(tidyr)
library(dplyr)

# Tendencia de la concentración de los metabolitos a lo largo del tiempo
# queremos un dataframe con el tipo de metabolito, día, réplica y concentración
# generamos un data frame con los datos
datos <- assays(SE)[[1]]
# renombramos la columna de metabolitos
datos$Metabolito <- rowData(SE)[,3]
# transformamos los datos en formato:
# tipo de metabolito, día, réplica y concentración
datos_long <- datos %>%
  pivot_longer(
    cols = -Metabolito,  # seleccionamos todas las columnas excepto "Metabolito"
    names_to = c("Día", "Réplica"),  # Nombres para las nuevas columnas
    names_pattern = ".*Day(\\d+)_Folate_(\\d+)",  # extraemos el día y réplica del nombre original
    values_to = "Concentración")
# convertimos "Día" y "Réplica" a variables numéricas
datos_long <- datos_long %>%
  mutate(Día = as.numeric(Día),
         Réplica = as.numeric(Réplica))
# representamos la tendencia de concentración
ggplot(datos_long, 
       aes(x = Día, y = Concentración, color = Metabolito, group = interaction(Metabolito, Réplica))) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1) +
  labs(title = "Tendencias de Concentración de Metabolitos a lo Largo del Tiempo",
       x = "Día",
       y = "Concentración") +
  theme_minimal()

# Análisis de los componentes princiales (PCA)
# queremos realizar un PCA de los diferentes metabolitos a diferentes días
# por si pudiera existir alguna relación
# para ello necesitamso los datos en formato columna, donde cada columna sea una condición, y hay dos primeras columnas con el día y la réplica
datos_wide <- datos_long %>%
  pivot_wider(names_from = Metabolito, values_from = Concentración) %>%
  arrange(Día, Réplica)  # ordenamos por día y réplica
# seleccionamos únicamente las columnas de metabolitos para el PCA
datos_pca <- datos_wide %>%
  select(-Día, -Réplica)  # excluimos las columnas "Día" y "Réplica"
# estandarizamos los datos antes del PCA 
datos_pca <- scale(datos_pca)
# realizamos el PCA
pca <- prcomp(datos_pca, center = TRUE, scale. = TRUE)
# vemos el del PCA para ver la varianza explicada
summary(pca)
# convertimos resultados del PCA en un data frame para representarlo mediante ggplot2
pca_data <- as.data.frame(pca$x)
pca_data$Día <- datos_wide$Día  # agregamos la información de día
pca_data$Réplica <- datos_wide$Réplica  # agregamos la información de réplica

# por último, representamos los puntos
ggplot(pca_data, aes(x = PC1, y = PC2, color = factor(Día))) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA de Concentraciones de Metabolitos",
       x = paste0("PC1 (", round(pca$sdev[1]^2 / sum(pca$sdev^2) * 100, 1), "% varianza)"),
       y = paste0("PC2 (", round(pca$sdev[2]^2 / sum(pca$sdev^2) * 100, 1), "% varianza)"),
       color = "Día") +
  theme_minimal()

# Boxplot de concentración de metabolitos por día 

# generamos un boxplot para cada metabolito
# donde mostremos la concentración a lo largo del tiempo
ggplot(datos_long, aes(x = factor(Día), y = Concentración, fill = Metabolito)) +
  geom_boxplot() +
  facet_wrap(~ Metabolito, scales = "free_y") +
  labs(title = "Distribución de Concentraciones de Metabolitos por Día",
       x = "Día",
       y = "Concentración") +
  theme_minimal()

# Heatmap  de Concentraciones Promedio de Metabolitos a lo largo de los días
# paquete necesario
library(reshape2)
# cambiamos formato de los datos para el heatmap
# calculamos la media de concentración de cada metabolito
datos_wide <- dcast(datos_long, Metabolito ~ Día, value.var = "Concentración", fun.aggregate = mean)
# creamos el heatmap
ggplot(melt(datos_wide), aes(x = variable, y = Metabolito, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Heatmap de Concentraciones Promedio de Metabolitos",
       x = "Día",
       y = "Metabolito") +
  theme_minimal()