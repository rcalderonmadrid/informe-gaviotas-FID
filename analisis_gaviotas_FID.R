# ══════════════════════════════════════════════════════════════════════════════
#  Influencia de la presencia humana y la densidad colonial en la respuesta a
#  las perturbaciones de la gaviota patiamarilla (Larus michahellis)
#
#  Código R del informe técnico (Máster en Bioestadística).
#  Autora: Rocío Calderón Madrid
#
#  El script reproduce, en el mismo orden que el informe, todos los análisis,
#  tablas, figuras y cálculos que aparecen en él:
#
#    0. Preparación: paquetes, datos y estética de las figuras
#    1. Análisis exploratorio (apartado 2.2; Figuras A1-A3)
#    2. Selección de la distribución y de la codificación de los visitantes
#       (apartado 3.1; Tabla 1, Figuras A4-A5, Tabla B1)
#    3. Evaluación de la linealidad mediante GAM
#       (apartado 3.2; Tabla 2, Tablas B2-B3)
#    4. Selección del modelo final: inferencia frecuentista y bayesiana
#       (apartado 3.3; Tabla 3, Figuras A6-A7, Tabla B4, Figura 1, Figura A8)
#    5. Cálculos de la Discusión (distancias de seguridad)
#    6. Información de la sesión
#
#  Requisitos: el fichero FID.dat en el directorio de trabajo y JAGS
#  instalado en el sistema (https://mcmc-jags.sourceforge.io/).
#  Las figuras se guardan en PDF en la carpeta "figuras/".
# ══════════════════════════════════════════════════════════════════════════════


# ══════════════════════════════════════════════════════════════════════════════
# 0. PREPARACIÓN
# ══════════════════════════════════════════════════════════════════════════════

# ── Paquetes ──────────────────────────────────────────────────────────────────
library(tidyverse)   # manipulación de datos y ggplot2
library(patchwork)   # composición de figuras
library(viridis)     # paletas de color
library(GGally)      # gráficos de pares (ggpairs)
library(e1071)       # coeficiente de asimetría
library(lawstat)     # prueba de Levene
library(boot)        # validación cruzada (cv.glm)
library(mgcv)        # modelos aditivos generalizados (GAM)
library(jagsUI)      # inferencia bayesiana con JAGS

# ── Datos ─────────────────────────────────────────────────────────────────────
# FID.dat contiene, para cada una de las 18 unidades de estudio (colonias):
#   colony   : identificador de la colonia
#   visitors : número medio de visitantes (2 = sin visitas regulares)
#   density  : densidad de parejas reproductoras (parejas / superficie)
#   meanFID  : distancia media de inicio de vuelo (FID, m)
fid <- read.table("FID.dat", header = TRUE)

# Variables transformadas y factor "visitada" (Sí si V > 2, No si V = 2)
data <- fid |>
  mutate(
    sqrt_fid  = sqrt(meanFID),
    log_fid   = log(meanFID),
    log_vis   = log(visitors),
    sqrt_dens = sqrt(density),
    visitada  = factor(ifelse(visitors > 2, "Sí", "No"), levels = c("No", "Sí"))
  )

n <- nrow(data)
table(data$visitada)   # 11 colonias no visitadas y 7 visitadas

# ── Estética común de las figuras ─────────────────────────────────────────────
tema <- theme_bw(base_size = 16) +
  theme(aspect.ratio = 1,
        legend.position = "none",
        panel.grid.minor = element_blank())

col_vis   <- c("No" = viridis(1, option = "mako",   begin = 0.45),
               "Sí" = viridis(1, option = "plasma", begin = 0.75))
col_loess <- viridis(1, option = "plasma", begin = 0.55)

dir.create("figuras", showWarnings = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# 1. ANÁLISIS EXPLORATORIO (apartado 2.2)
# ══════════════════════════════════════════════════════════════════════════════

# ── 1.1. Distribución de la variable respuesta (Figura A1) ────────────────────

# Coeficiente de asimetría de la FID media (texto: 0,78)
skewness(data$meanFID)

# Histograma + densidad kernel + posición de cada colonia (rug)
hist_dens <- function(df, var, xlab, col, bins = 8) {
  ggplot(df, aes(x = {{ var }})) +
    geom_histogram(aes(y = after_stat(density)), bins = bins,
                   fill = col, colour = "white", alpha = 0.7) +
    geom_density(colour = col, linewidth = 1.1, adjust = 1.2) +
    geom_rug(colour = "grey30", alpha = 0.7, length = unit(0.03, "npc")) +
    labs(x = xlab, y = "Densidad") +
    tema
}

col_fid <- viridis(3, option = "viridis", begin = 0.15, end = 0.75)

# Figura A1: FID en escala original (a), raíz cuadrada (b) y logarítmica (c)
fig_A1 <- hist_dens(data, meanFID,  "FID (m)",             col_fid[1]) +
  hist_dens(data, sqrt_fid, expression(sqrt(FID)), col_fid[2]) +
  hist_dens(data, log_fid,  expression(log(FID)),  col_fid[3]) +
  plot_annotation(tag_levels = "a")
fig_A1
ggsave("figuras/expl_E1.pdf", fig_A1, width = 15, height = 5.5)

# ── 1.2. Relación entre la FID y las covariables (Figuras A2 y A3) ────────────
# Gráficos de pares (GGally::ggpairs) con la respuesta en la última columna:
#   - triángulo inferior: puntos coloreados según "visitada", recta de mínimos
#     cuadrados con IC 95 % (gris) y suavizado LOESS (discontinuo);
#   - diagonal: histograma y densidad kernel;
#   - triángulo superior: correlación de Pearson sobre fondo coloreado.
pal_cor <- colorRampPalette(hcl.colors(11, palette = "Blue-Red 3"))(201)

panel_inf <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "grey20",
                fill = "grey70", linewidth = 0.8, alpha = 0.25) +
    geom_smooth(method = "loess", formula = y ~ x, span = 0.9, se = FALSE,
                colour = col_loess, linetype = "dashed", linewidth = 0.9) +
    geom_point(aes(fill = visitada), shape = 21, size = 2.6,
               colour = "grey20", alpha = 0.9) +
    scale_fill_manual(values = col_vis, name = "Visitada")
}

panel_diag <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_histogram(aes(y = after_stat(density)), bins = 8,
                   fill = viridis(1, option = "viridis", begin = 0.45),
                   colour = "white", alpha = 0.7) +
    geom_density(colour = viridis(1, option = "viridis", begin = 0.25),
                 linewidth = 1, adjust = 1.2)
}

panel_sup <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  r <- cor(x, y)
  ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = pal_cor[round((r + 1) * 100) + 1], alpha = 0.9) +
    annotate("text", x = 0.5, y = 0.5, label = sprintf("r = %.2f", r),
             size = 6, fontface = "bold",
             colour = ifelse(abs(r) > 0.6, "white", "grey15")) +
    xlim(0, 1) + ylim(0, 1)
}

pares <- function(cols, etiquetas) {
  ggpairs(data, columns = cols,
          columnLabels = etiquetas, labeller = label_parsed,
          lower = list(continuous = panel_inf),
          diag  = list(continuous = panel_diag),
          upper = list(continuous = panel_sup),
          legend = c(2, 1)) +
    theme_bw(base_size = 16) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92"),
          strip.text = element_text(size = 15),
          legend.position = "bottom",
          legend.text  = element_text(size = 13),
          legend.title = element_text(size = 13))
}

# Figura A2: enfoque LM, respuesta sqrt(FID)
#   (a) covariables originales; (b) covariables transformadas
fig_A2a <- pares(c("visitors", "density", "sqrt_fid"),
                 c("Visitantes", "Densidad", "sqrt(FID)"))
fig_A2b <- pares(c("log_vis", "sqrt_dens", "sqrt_fid"),
                 c("log(Visitantes)", "sqrt(Densidad)", "sqrt(FID)"))
fig_A2a; fig_A2b
ggsave("figuras/expl_L1.pdf", ggmatrix_gtable(fig_A2a), width = 10, height = 10.5)
ggsave("figuras/expl_L2.pdf", ggmatrix_gtable(fig_A2b), width = 10, height = 10.5)

# Figura A3: enfoque GLM Gamma (enlace log), escala del predictor log(FID)
#   (a) covariables originales; (b) covariables transformadas
fig_A3a <- pares(c("visitors", "density", "log_fid"),
                 c("Visitantes", "Densidad", "log(FID)"))
fig_A3b <- pares(c("log_vis", "sqrt_dens", "log_fid"),
                 c("log(Visitantes)", "sqrt(Densidad)", "log(FID)"))
fig_A3a; fig_A3b
ggsave("figuras/expl_G1.pdf", ggmatrix_gtable(fig_A3a), width = 10, height = 10.5)
ggsave("figuras/expl_G2.pdf", ggmatrix_gtable(fig_A3b), width = 10, height = 10.5)

# Correlaciones citadas en el texto
with(data, c(
  r_logV_sqrtY = cor(log_vis,   sqrt_fid),    # -0,61
  r_logV_logY  = cor(log_vis,   log_fid),     # -0,64
  r_sqrtD_sqrtY = cor(sqrt_dens, sqrt_fid),   # -0,23
  r_sqrtD_logY  = cor(sqrt_dens, log_fid),    # -0,18
  r_logV_sqrtD  = cor(log_vis,   sqrt_dens))) # -0,40 (entre covariables)

# Hueco en el número de visitantes: ningún dato entre 154 y 2643
sort(unique(data$visitors))


# ══════════════════════════════════════════════════════════════════════════════
# 2. SELECCIÓN DE LA DISTRIBUCIÓN Y DE LA CODIFICACIÓN DE LOS VISITANTES
#    (apartado 3.1)
# ══════════════════════════════════════════════════════════════════════════════

# ── 2.1. Ajuste de los cuatro modelos ─────────────────────────────────────────
# LM sobre sqrt(FID) (modelo de Martínez-Abraín et al., 2008)
lm_c <- lm(sqrt_fid ~ log_vis  + sqrt_dens, data = data)   # visitantes continua
lm_d <- lm(sqrt_fid ~ visitada + sqrt_dens, data = data)   # factor visitada

# GLM Gamma con enlace logarítmico sobre la FID
glm_c <- glm(meanFID ~ log_vis  + sqrt_dens, family = Gamma(link = "log"),
             data = data)
glm_d <- glm(meanFID ~ visitada + sqrt_dens, family = Gamma(link = "log"),
             data = data)

summary(lm_c); summary(lm_d)
summary(glm_c); summary(glm_d)

# ── 2.2. AIC en una escala común (corrección del jacobiano) ───────────────────
# Si Z = sqrt(Y), f_Y(y) = f_Z(sqrt(y)) / (2 sqrt(y)), por lo que el AIC del LM
# en la escala de la FID es AIC(LM) + J, con J = 2 * sum(log(2 * sqrt(Y))).
J <- 2 * sum(log(2 * sqrt(data$meanFID)))
J

aic_fid <- c(lm_c  = AIC(lm_c) + J, lm_d  = AIC(lm_d) + J,
             glm_c = AIC(glm_c),    glm_d = AIC(glm_d))

# ── 2.3. Validación cruzada LOOCV con cv.glm (MSEP ajustado, delta[2]) ────────
# Los LM se reajustan con glm(family = gaussian), que es equivalente a lm().
# En los LM, respuesta y predicción se elevan al cuadrado para expresar el
# error en la escala de la FID (m^2).
lm_c_glm <- glm(sqrt_fid ~ log_vis  + sqrt_dens, family = gaussian, data = data)
lm_d_glm <- glm(sqrt_fid ~ visitada + sqrt_dens, family = gaussian, data = data)

coste_fid <- function(y, yhat) mean((y^2 - yhat^2)^2)

cv_lm_c  <- cv.glm(data, lm_c_glm, cost = coste_fid, K = n)
cv_lm_d  <- cv.glm(data, lm_d_glm, cost = coste_fid, K = n)
cv_glm_c <- cv.glm(data, glm_c, K = n)
cv_glm_d <- cv.glm(data, glm_d, K = n)

msep <- c(lm_c  = cv_lm_c$delta[2],  lm_d  = cv_lm_d$delta[2],
          glm_c = cv_glm_c$delta[2], glm_d = cv_glm_d$delta[2])

# ── 2.4. Tabla 1 ──────────────────────────────────────────────────────────────
# R2 de los LM y devianza explicada de los GLM (%)
dev_expl <- function(m) 100 * (1 - m$deviance / m$null.deviance)

tabla1 <- data.frame(
  Modelo    = c("LM", "LM", "GLM Gamma", "GLM Gamma"),
  Predictor = c("log V + sqrt(D)", "W + sqrt(D)", "log V + sqrt(D)", "W + sqrt(D)"),
  R2_DevExpl = c(100 * summary(lm_c)$r.squared, 100 * summary(lm_d)$r.squared,
                 dev_expl(glm_c), dev_expl(glm_d)),
  AIC   = aic_fid,
  dAIC  = aic_fid - min(aic_fid),
  MSEP  = msep)
print(tabla1, digits = 5, row.names = FALSE)

# ── 2.5. Diagnóstico de los residuos (Figuras A4 y A5) ────────────────────────
# Figura A4: LM con visitantes continua
pdf("figuras/diag_lm_c.pdf", width = 9, height = 8)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1))
plot(lm_c)
dev.off()

# Figura A5: GLM Gamma con visitantes continua
pdf("figuras/diag_glm_c.pdf", width = 9, height = 8)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1))
plot(glm_c)
dev.off()

# Distancia de Cook máxima (observación más influyente)
max(cooks.distance(lm_c)); max(cooks.distance(glm_c))

# ── 2.6. Tabla B1: Kolmogorov-Smirnov y Levene en los LM ──────────────────────
# KS: residuos estudentizados frente a una t de Student con n - p - 1 = 14 g.l.
# Levene: residuos estandarizados agrupados por cuartiles de los ajustados.
grupos_ajuste <- function(m) {
  cut(fitted(m), quantile(fitted(m), (0:4) / 4), include.lowest = TRUE)
}

validacion_lm <- function(m) {
  ks  <- ks.test(rstudent(m), "pt", df = df.residual(m) - 1)
  lev <- levene.test(rstandard(m), grupos_ajuste(m))
  c(KS_D = unname(ks$statistic), KS_p = ks$p.value,
    Levene_F = unname(lev$statistic), Levene_p = lev$p.value)
}

tablaB1 <- rbind(LM_logV = validacion_lm(lm_c), LM_W = validacion_lm(lm_d))
round(tablaB1, 3)


# ══════════════════════════════════════════════════════════════════════════════
# 3. EVALUACIÓN DE LA LINEALIDAD MEDIANTE GAM (apartado 3.2)
# ══════════════════════════════════════════════════════════════════════════════

# ── 3.1. Ajuste de los tres GAM (Gamma, enlace log, covariables originales) ───
# Splines de placas delgadas (bs = "tp") y suavizado estimado por REML.
# Se parte de k = 5 y se reduce la base de la densidad (k = 3) hasta
# incluirla como término lineal, ya que sus edf son próximos a 1.
gam_k5  <- gam(meanFID ~ s(visitors, bs = "tp", k = 5) + s(density, bs = "tp", k = 5),
               family = Gamma(link = "log"), data = data, method = "REML")
gam_k3  <- gam(meanFID ~ s(visitors, bs = "tp", k = 5) + s(density, bs = "tp", k = 3),
               family = Gamma(link = "log"), data = data, method = "REML")
gam_lin <- gam(meanFID ~ s(visitors, bs = "tp", k = 5) + density,
               family = Gamma(link = "log"), data = data, method = "REML")

gams <- list(gam_k5 = gam_k5, gam_k3 = gam_k3, gam_lin = gam_lin)
lapply(gams, summary)

# ── 3.2. Validación cruzada LOOCV de los GAM ──────────────────────────────────
# cv.glm no admite objetos gam, así que se reproduce el mismo cálculo
# (LOOCV con la corrección de sesgo de delta[2]) en una función propia.
cv_loo <- function(ajustar) {
  y  <- data$meanFID
  m0 <- ajustar(data)
  c0 <- mean((y - predict(m0, newdata = data, type = "response"))^2)
  e  <- ci <- numeric(n)
  for (i in seq_len(n)) {
    m     <- ajustar(data[-i, ])
    pred  <- predict(m, newdata = data, type = "response")
    e[i]  <- y[i] - pred[i]
    ci[i] <- mean((y - pred)^2)
  }
  d1 <- mean(e^2)
  c(delta1 = d1, delta2 = d1 + c0 - mean(ci))
}

# Comprobación: aplicada al GLM Gamma da exactamente lo mismo que cv.glm
cv_loo(function(d) glm(meanFID ~ log_vis + sqrt_dens,
                       family = Gamma(link = "log"), data = d))
cv_glm_c$delta

cv_gams <- sapply(gams, function(g)
  cv_loo(function(d) gam(formula(g), family = Gamma(link = "log"),
                         data = d, method = "REML")))

# ── 3.3. Tabla 2 ──────────────────────────────────────────────────────────────
edf_term <- function(m, v) {
  e <- m$edf[grep(paste0("s\\(", v), names(m$edf))]
  if (length(e) == 0) NA else sum(e)
}

tabla2 <- data.frame(
  Predictor = c("f1(V) + f2(D)", "f1(V) + f2(D)", "f1(V) + b2 D"),
  kV = c(5, 5, 5), kD = c("5", "3", "-"),
  edf_V = sapply(gams, edf_term, v = "visitors"),
  edf_D = sapply(gams, edf_term, v = "density"),   # NA = término lineal
  DevExpl = sapply(gams, function(m) 100 * summary(m)$dev.expl),
  AIC  = sapply(gams, AIC),
  MSEP = cv_gams["delta2", ])
print(tabla2, digits = 5, row.names = FALSE)

# ── 3.4. Tabla B2: dimensión de las bases (k.check) ───────────────────────────
# El p-valor del k-index se obtiene por aleatorización: se fija la semilla.
lapply(gams, function(m) { set.seed(2026); k.check(m) })

# ── 3.5. Tabla B3: concurvidad de los splines ─────────────────────────────────
# Se muestran las tres medidas (worst, observed, estimate) de cada término
# suave. En gam_lin, la columna "para" (parte paramétrica) depende de que la
# densidad no esté centrada y no indica concurvidad real, por lo que no se
# incluye en la tabla.
lapply(gams, function(m) round(concurvity(m, full = TRUE), 2))


# ══════════════════════════════════════════════════════════════════════════════
# 4. SELECCIÓN DEL MODELO FINAL: INFERENCIA FRECUENTISTA Y BAYESIANA
#    (apartado 3.3)
# ══════════════════════════════════════════════════════════════════════════════

# ── 4.1. GAM frente al GLM Gamma ──────────────────────────────────────────────
# El GAM tiene menor AIC, pero mayor error de predicción (sobreajuste)
c(AIC_GAM = AIC(gam_lin), AIC_GLM = AIC(glm_c))
c(MSEP_GAM = cv_gams["delta2", "gam_lin"], MSEP_GLM = msep[["glm_c"]])
100 * (cv_gams["delta2", "gam_lin"] / msep[["glm_c"]] - 1)   # % mayor (28 %)

# Modelo final: GLM Gamma con log(V) y sqrt(D)
dev_expl(glm_c)   # devianza explicada (60,2 %)

# ── 4.2. Inferencia por máxima verosimilitud ──────────────────────────────────
# Estimaciones, p-valores e IC 95 % por verosimilitud perfil
summary(glm_c)$coefficients
ic_mv <- suppressMessages(confint(glm_c))
mv <- data.frame(Estimacion = coef(glm_c), Inf95 = ic_mv[, 1], Sup95 = ic_mv[, 2])

# ── 4.3. Inferencia bayesiana con JAGS ────────────────────────────────────────
# Verosimilitud Gamma parametrizada por la media mu y la forma (shape):
#   FID_i ~ Gamma(shape, shape / mu_i),  log(mu_i) = b0 + b1 log V + b2 sqrt(D)
# A priori vagas: b_j ~ N(0, precisión 0,001) y shape ~ Gamma(0,01; 0,01).
# La dispersión del GLM es phi = 1 / shape.
cat(file = "ModeloGaviotas", "model{
  # Verosimilitud
  for (i in 1:n) {
    fid[i] ~ dgamma(shape, rate[i])
    rate[i] <- shape / mu[i]
    log(mu[i]) <- b0 + b1 * log_vis[i] + b2 * sqrt_dens[i]
  }
  # Distribuciones a priori vagas
  b0 ~ dnorm(0, 0.001)
  b1 ~ dnorm(0, 0.001)
  b2 ~ dnorm(0, 0.001)
  shape  ~ dgamma(0.01, 0.01)
  phi <- 1 / shape
}")

datos_jags <- list(fid = data$meanFID, log_vis = data$log_vis,
                   sqrt_dens = data$sqrt_dens, n = n)

# Valores iniciales dispersos y semilla fija por cadena (resultados
# reproducibles, también con parallel = TRUE)
iniciales <- list(
  list(b0 = 3, b1 =  0.0, b2 =  0.00, shape = 2,
       .RNG.name = "base::Mersenne-Twister", .RNG.seed = 101),
  list(b0 = 4, b1 = -0.2, b2 = -0.05, shape = 5,
       .RNG.name = "base::Mersenne-Twister", .RNG.seed = 202),
  list(b0 = 5, b1 = -0.4, b2 = -0.10, shape = 10,
       .RNG.name = "base::Mersenne-Twister", .RNG.seed = 303))

# 3 cadenas de 12000 iteraciones, 2000 de calentamiento y adelgazamiento 10
# (3000 muestras a posteriori en total)
ajuste_jags <- jags(data = datos_jags, inits = iniciales,
                    parameters.to.save = c("b0", "b1", "b2", "phi"),
                    model.file = "ModeloGaviotas",
                    n.chains = 3, n.iter = 12000, n.burnin = 2000, n.thin = 10,
                    parallel = TRUE)
ajuste_jags

# Convergencia: R-hat de Gelman-Rubin (< 1,1) y tamaño muestral efectivo (> 100)
unlist(ajuste_jags$Rhat[c("b0", "b1", "b2", "phi")])
unlist(ajuste_jags$n.eff[c("b0", "b1", "b2", "phi")])

# Figura A6: trazas de las tres cadenas
pdf("figuras/jags_trazas.pdf", width = 9, height = 10)
traceplot(ajuste_jags)
dev.off()

# Muestras a posteriori (3000 x 4)
mc <- do.call(cbind, ajuste_jags$sims.list[c("b0", "b1", "b2", "phi")])

# Media a posteriori e intervalo de credibilidad del 95 %
bayes <- data.frame(
  Estimacion = colMeans(mc[, c("b0", "b1", "b2")]),
  Inf95 = apply(mc[, c("b0", "b1", "b2")], 2, quantile, 0.025),
  Sup95 = apply(mc[, c("b0", "b1", "b2")], 2, quantile, 0.975))

# ── 4.4. Tabla 3: máxima verosimilitud frente a inferencia bayesiana ──────────
tabla3 <- rbind(mv, bayes)
tabla3 <- data.frame(Metodo = rep(c("Máxima verosimilitud", "Bayesiana"), each = 3),
                     Parametro = rep(c("Intercepto", "log V", "sqrt(D)"), 2),
                     round(tabla3, 3), row.names = NULL)
tabla3

# ── 4.5. Figura A7: distribuciones a posteriori ───────────────────────────────
etiq_par <- c(b0 = "beta[0]~(Intercepto)", b1 = "beta[1]~(log~Visitantes)",
              b2 = "beta[2]~(sqrt(Densidad))", phi = "phi~(Dispersión)")
post_df <- data.frame(valor = c(mc), par = rep(colnames(mc), each = nrow(mc))) |>
  mutate(par = factor(par, levels = names(etiq_par), labels = etiq_par))
res_post <- post_df |>
  group_by(par) |>
  summarise(media = mean(valor), inf = quantile(valor, 0.025),
            sup = quantile(valor, 0.975), .groups = "drop") |>
  mutate(mv = c(coef(glm_c), summary(glm_c)$dispersion))

fig_A7 <- ggplot(post_df, aes(valor)) +
  geom_rect(data = res_post, aes(xmin = inf, xmax = sup, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = viridis(1, begin = 0.45), alpha = 0.15) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = viridis(1, begin = 0.45), colour = "white", alpha = 0.6) +
  geom_density(colour = viridis(1, begin = 0.25), linewidth = 1) +
  geom_vline(data = res_post, aes(xintercept = media), linewidth = 0.8) +
  geom_vline(data = res_post, aes(xintercept = mv), colour = "#E74C3C",
             linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ par, scales = "free", labeller = label_parsed) +
  labs(x = "Valor del parámetro", y = "Densidad a posteriori") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"))
fig_A7
ggsave("figuras/jags_posterior.pdf", fig_A7, width = 10, height = 7)

# ── 4.6. Interpretación de los efectos ────────────────────────────────────────
# Probabilidad a posteriori de que cada efecto sea negativo
colMeans(mc[, c("b1", "b2")] < 0)                 # 0,999 y 0,993

# Duplicar el número de visitantes multiplica la FID media por 2^b1
dup <- 2^mc[, "b1"]
c(factor = 2^mean(mc[, "b1"]), quantile(dup, c(0.025, 0.975)))   # 0,896 (0,854-0,942)

# Pasar de 100 a 200 parejas por unidad de superficie: reducción (%)
red_dens <- 100 * (1 - exp(mc[, "b2"] * (sqrt(200) - sqrt(100))))
c(media = mean(red_dens), quantile(red_dens, c(0.025, 0.975)))   # 27 % (7-42 %)

# ── 4.7. Tabla B4: importancia relativa (coeficientes estandarizados) ─────────
# Estandarizar las covariables transformadas (media 0, DT 1) da el mismo
# modelo con beta*_j = beta_j * s_j. Por máxima verosimilitud se reajusta
# el GLM; en el enfoque bayesiano se reescalan las muestras a posteriori.
s_vis  <- sd(data$log_vis)
s_dens <- sd(data$sqrt_dens)

data$z_log_vis   <- as.numeric(scale(data$log_vis))
data$z_sqrt_dens <- as.numeric(scale(data$sqrt_dens))
glm_z <- glm(meanFID ~ z_log_vis + z_sqrt_dens, family = Gamma(link = "log"),
             data = data)
ic_z  <- suppressMessages(confint(glm_z))

b1_z <- mc[, "b1"] * s_vis
b2_z <- mc[, "b2"] * s_dens

tablaB4 <- data.frame(
  Metodo   = rep(c("Máxima verosimilitud", "Bayesiana"), each = 2),
  Variable = rep(c("log V", "sqrt(D)"), 2),
  beta_z = c(coef(glm_z)[-1], mean(b1_z), mean(b2_z)),
  Inf95  = c(ic_z[-1, 1], quantile(b1_z, 0.025), quantile(b2_z, 0.025)),
  Sup95  = c(ic_z[-1, 2], quantile(b1_z, 0.975), quantile(b2_z, 0.975)))
tablaB4$exp_beta_z <- exp(tablaB4$beta_z)
print(tablaB4, digits = 3, row.names = FALSE)

# Reducción de la FID media por cada DT (%): 43 % (log V) y 30 % (sqrt D)
100 * (1 - exp(coef(glm_z)[-1]))

# Cociente de importancia |b1*| / |b2*|: máxima verosimilitud y bayesiano
abs(coef(glm_z)["z_log_vis"] / coef(glm_z)["z_sqrt_dens"])        # 1,6
ratio_z <- abs(b1_z / b2_z)
quantile(ratio_z, c(0.025, 0.5, 0.975))                           # 0,8-5,6
mean(abs(b1_z) > abs(b2_z))                                       # 0,93

# ── 4.8. FID ajustada (Figura 1) ──────────────────────────────────────────────
# Rejilla regular en log(visitantes) y densidad que cubre el rango observado.
# Para cada muestra a posteriori se calcula mu = exp(b0 + b1 log V + b2 sqrt D)
# y se resume con la media a posteriori.
rejilla <- expand.grid(
  visitors = exp(seq(log(2), log(7500), length.out = 80)),
  density  = seq(11, 445, length.out = 80))
X_rej  <- cbind(1, log(rejilla$visitors), sqrt(rejilla$density))
mu_rej <- exp(X_rej %*% t(mc[, c("b0", "b1", "b2")]))   # puntos x muestras
rejilla$fid_media <- rowMeans(mu_rej)

# FID ajustada en cada colonia observada (media a posteriori)
X_obs  <- cbind(1, data$log_vis, data$sqrt_dens)
mu_obs <- exp(X_obs %*% t(mc[, c("b0", "b1", "b2")]))
data$fid_bayes <- rowMeans(mu_obs)

lim_fid    <- c(0, 80)
escala_fid <- scale_fill_viridis_c(option = "viridis", begin = 0.1, end = 0.95,
                                   limits = lim_fid, name = "FID (m)")
eje_x <- scale_x_continuous(breaks = 0:9)

# (a) FID observada
p_real <- ggplot(data, aes(log(visitors), density)) +
  geom_point(aes(size = meanFID, fill = meanFID), shape = 21, colour = "grey20") +
  escala_fid + eje_x + scale_size_continuous(range = c(2, 9), guide = "none") +
  labs(x = "log(Nº de visitantes)", y = "Densidad (parejas/superficie)",
       title = "FID real") + tema + theme(legend.position = "right")

# (b) FID ajustada (media a posteriori) con isolíneas cada 10 m
p_ajust <- ggplot(rejilla, aes(log(visitors), density)) +
  geom_raster(aes(fill = fid_media), interpolate = TRUE) +
  geom_contour(aes(z = fid_media), colour = "white", alpha = 0.7, binwidth = 10) +
  geom_point(data = data, aes(log(visitors), density), shape = 21, fill = "white",
             colour = "grey20", size = 1.8) +
  escala_fid + eje_x +
  labs(x = "log(Nº de visitantes)", y = "Densidad (parejas/superficie)",
       title = "FID ajustada (media a posteriori)") + tema +
  theme(legend.position = "right")

fig_1 <- p_real + p_ajust + plot_annotation(tag_levels = "a")
fig_1
ggsave("figuras/fig_superficie.pdf", fig_1, width = 12, height = 5.5)

# Valores citados en el texto: FID estimada en las colonias observadas.
# En el texto se calculan con la media a posteriori de los coeficientes
# (exp(media(b0) + media(b1) log V + media(b2) sqrt D)). Promediando la FID
# sobre las muestras a posteriori (fid_bayes, la escala de la Figura 1) los
# valores son algo mayores (63-17 m y 12-21 m), sin cambiar las conclusiones.
b_media <- colMeans(mc[, c("b0", "b1", "b2")])
data$fid_texto <- as.numeric(exp(X_obs %*% b_media))
data |>
  select(colony, visitors, density, meanFID, fid_texto, fid_bayes) |>
  arrange(visitors, density)
# Sin visitas: de ~61 m (menor densidad) a ~16 m (mayor densidad)
range(data$fid_texto[data$visitors == 2])
# Más de 2600 visitantes: entre ~11 y ~20 m
range(data$fid_texto[data$visitors > 2600])

# ── 4.9. Figura A8: representación tridimensional ─────────────────────────────
x_lv <- seq(log(2), log(7500), length.out = 40)
y_d  <- seq(11, 445, length.out = 40)
z_fid <- outer(x_lv, y_d, function(lv, d)
  rowMeans(exp(cbind(1, lv, sqrt(d)) %*% t(mc[, c("b0", "b1", "b2")]))))

# Colores de la superficie y de los puntos según la FID
z_fac   <- (z_fid[-1, -1] + z_fid[-1, -40] + z_fid[-40, -1] + z_fid[-40, -40]) / 4
# (pmin limita la FID a 79,9 m para que ninguna faceta quede fuera de la escala)
col_sup <- viridis(100, begin = 0.1, end = 0.95)[cut(pmin(z_fac, 79.9), seq(0, 80, length.out = 101))]
col_pts <- viridis(100, begin = 0.1, end = 0.95)[cut(data$meanFID, seq(0, 80, length.out = 101))]

pdf("figuras/fig_3d.pdf", width = 12, height = 6)
par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))

# (a) FID observada en cada colonia
pm <- persp(x_lv, y_d, matrix(NA, 40, 40), zlim = c(0, 80), theta = 40, phi = 20,
            xlab = "log(Visitantes)", ylab = "Densidad", zlab = "FID (m)",
            ticktype = "detailed", cex.axis = 0.7, main = "(a) FID real")
base <- trans3d(data$log_vis, data$density, 0, pm)
top  <- trans3d(data$log_vis, data$density, data$meanFID, pm)
segments(base$x, base$y, top$x, top$y, col = "grey50")
points(top, pch = 21, bg = col_pts, cex = 1.6)

# (b) Superficie ajustada con las colonias observadas
pm <- persp(x_lv, y_d, z_fid, zlim = c(0, 80), theta = 40, phi = 20,
            col = col_sup, border = adjustcolor("grey30", 0.3),
            xlab = "log(Visitantes)", ylab = "Densidad", zlab = "FID (m)",
            ticktype = "detailed", cex.axis = 0.7,
            main = "(b) FID ajustada (media a posteriori)")
sup <- trans3d(data$log_vis, data$density, data$fid_bayes, pm)
obs <- trans3d(data$log_vis, data$density, data$meanFID, pm)
segments(sup$x, sup$y, obs$x, obs$y, col = "grey20")
points(obs, pch = 21, bg = "white", cex = 1.3)
dev.off()


# ══════════════════════════════════════════════════════════════════════════════
# 5. CÁLCULOS DE LA DISCUSIÓN
# ══════════════════════════════════════════════════════════════════════════════

# Distancias de seguridad iguales al doble de la FID media (Martínez-Abraín et
# al., 2008), para una densidad de 80 parejas y distinto número de visitantes:
# media a posteriori e intervalo de credibilidad del 95 %
dist_seg <- t(sapply(c(2, 100, 1000, 5000), function(V) {
  x <- 2 * exp(mc[, "b0"] + mc[, "b1"] * log(V) + mc[, "b2"] * sqrt(80))
  c(Visitantes = V, Media = mean(x), quantile(x, c(0.025, 0.975)))
}))
round(dist_seg)
# Sin visitas: ~84 m (63-112 m); 5000 visitantes: ~25 m (16-38 m)


# ══════════════════════════════════════════════════════════════════════════════
# 6. INFORMACIÓN DE LA SESIÓN
# ══════════════════════════════════════════════════════════════════════════════
sessionInfo()
