# Distancia de inicio de vuelo (FID) de la gaviota patiamarilla

Código R y JAGS del informe técnico *Influencia de la presencia humana y la densidad colonial en la respuesta a las perturbaciones de la gaviota patiamarilla* (Máster en Bioestadística, Rocío Calderón Madrid, 2026).

## Contenido

| Archivo | Descripción |
|---|---|
| `informe_gaviotas_FID.pdf` | Informe técnico completo. |
| `analisis_gaviotas_FID.R` | Script con todos los análisis, tablas, figuras y cálculos del informe, en el mismo orden que el documento. |
| `FID.dat` | Datos: FID media, número de visitantes y densidad de parejas de las 18 unidades de estudio. |
| `figuras/` | Figuras del informe en PDF (las genera el script). |

## Cómo reproducir el análisis

1. Instalar [JAGS](https://mcmc-jags.sourceforge.io/).
2. Instalar los paquetes de R:
   ```r
   install.packages(c("tidyverse", "patchwork", "viridis", "GGally", "e1071",
                      "lawstat", "boot", "mgcv", "jagsUI"))
   ```
3. Situar el directorio de trabajo en la carpeta que contiene `FID.dat` y ejecutar `analisis_gaviotas_FID.R`.

El script crea el archivo del modelo bayesiano (`ModeloGaviotas`) y la carpeta `figuras/`. Los resultados son reproducibles: las cadenas de JAGS y la prueba `k.check` usan semillas fijas.

## Datos

Los datos proceden de Martínez-Abraín, A., Oro, D., Conesa, D. y Jiménez, J. (2008). Compromise between seabird enjoyment and disturbance: the role of observed and observers. *Environmental Conservation*, 35(2), 104–108. La licencia del repositorio se aplica únicamente al código.
