# --------------------------------------------------------------------------------------
# 1ª Avaliação - Estatística Computacional
# Modelo: Compound Weibull
# Aluno: Yuliana Apaza Flores
# --------------------------------------------------------------------------------------

# Pacotes necessários
library(ggplot2)
library(patchwork)
library(knitr)
library(kableExtra)
library(dplyr)
library(magrittr)
library(tidyr)
library(tibble)
library(goftest) 
library(scales)  

# --------------------------------------------------------------------------------------
# 1. Especificação de parâmetros e tamanhos amostrais:
# --------------------------------------------------------------------------------------
# Definição dos 6 cenários para o modelo Compound Weibull, combinando diferentes
# valores de beta e rho, com um alpha fixo.
# Amostras de tamanhos n = 10, 100, 1.000, 10.000, 100.000.
# --------------------------------------------------------------------------------------

parametros <- expand.grid(n = c(10, 100, 1000, 10000, 100000),
                          alpha = c(2.0),
                          beta  = c(0.1, 0.5),
                          rho   = c(0.25, 0.5, 0.75))
# --------------------------------------------------------------------------------------
# 2. Geração dos valores pseudoaleatórios:
# --------------------------------------------------------------------------------------
# A função quantil do modelo Compound Weibull não possui uma forma analítica fechada.
# Portanto, a geração será feita pelo Método da Transformação Inversa,
# utilizando uma aproximação numérica para a função quantil via 'uniroot'.
# --------------------------------------------------------------------------------------

# Função de Densidade de Probabilidade (PDF)
dcompound_weibull <- function(x, rho, alpha, beta) {
  numerador <- alpha * beta^alpha * (1 - rho) * x^(alpha - 1) * exp(-(beta * x)^alpha)
  denominador <- (1 - rho * exp(-(beta * x)^alpha))^2
  return(numerador / denominador)
}

# Função de Distribuição Acumulada (CDF)
pcompound_weibull <- function(x, rho, alpha, beta) {
  ifelse(x <= 0, 0, {
    numerador <- 1 - exp(-(beta * x)^alpha)
    denominador <- 1 - rho * exp(-(beta * x)^alpha)
    return(numerador / denominador)
  })
}

# Função de Densidade de Probabilidade (PDF) de la propuesta (No se usa en el resto del script, pero se mantiene la estructura)
dweibull_generalizada <- function(x, rho, alpha, beta) {
  weibull_ge <- alpha * beta^alpha * (1 - rho) * x^(alpha - 1) * exp(-(beta * x)^alpha)
  return(weibull_ge)
}

# Função Quantil Numérica - F^-1(u)
# Encontra a raiz da equação F(x) - u = 0 para um dado u.
qcompound_weibullN <- function(u, rho, alpha, beta, lower = 0.0001, upper = 1000) {
  sapply(u, function(p) {
    if (p == 0) return(0)
    if (p >= 1) return(Inf)
    uniroot(function(x) pcompound_weibull(x, rho, alpha, beta) - p,
            interval = c(lower, upper))$root
  })
}

# Função Geradora de Valores Pseudoaleatórios via Inversão Numérica
rcompound_weibullN <- function(n, rho, alpha, beta) {
  u <- runif(n)
  qcompound_weibullN(u, rho, alpha, beta)
}

# --------------------------------------------------------------------------------------
# 3. Verificação dos geradores
# --------------------------------------------------------------------------------------
# - Testes de aderência (Kolmogorov-Smirnov e Cramér-von Mises).
# - Ilustração gráfica (Histograma com densidade e CDF empírica vs. teórica).
# --------------------------------------------------------------------------------------

# ----- SEÇÃO 3.1: Análise Gráfica -----
cat("### Gráficos: Análise Visual das Amostras Geradas\n\n")

params_unicos <- unique(parametros[, c("alpha", "beta", "rho")])
n_valores <- unique(parametros$n)

plot_rows <- list()

for (i in 1:nrow(params_unicos)) {
  a <- params_unicos$alpha[i]
  b <- params_unicos$beta[i]
  r <- params_unicos$rho[i]
  row_title_text <- bquote(paste(alpha, " = ", .(a), "  |  ", beta, " = ", .(b), "  |  ", rho, " = ", .(r)))
  g_row_title <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, 
             label = paste(deparse(row_title_text, width.cutoff = 500L), collapse = ""), 
             parse = TRUE, size = 3.5, fontface = "bold") +
    theme_void()
  row_plot_pairs <- list()
  
  for (n_val in n_valores) {
    set.seed(42)
    amostra <- rcompound_weibullN(n_val, r, a, b)
    # ---- CORREÇÃO 2 (de novo): Remover valores NA antes de plotar ----
    df_amostra <- data.frame(x = na.omit(amostra)) 
    
    g1 <- ggplot(df_amostra, aes(x = x)) +
      geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
      stat_function(fun = dcompound_weibull, args = list(rho = r, alpha = a, beta = b),
                    aes(color = "Teórica"), linewidth = 0.8) +
      scale_color_manual(name = "Densidade", values = c("Teórica" = "navy")) +
      labs(title = paste("n =", n_val), x = NULL, y = "Densidade") + 
      theme_minimal() +
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))
    
    g2 <- ggplot(df_amostra, aes(x = x)) +
      stat_ecdf(aes(color = "Empírica (Amostra)"), geom = "step", linewidth = 0.8) +
      stat_function(fun = pcompound_weibull, args = list(rho = r, alpha = a, beta = b),
                    aes(color = "Analítica (Teórica)"), linewidth = 0.8, linetype = "dashed") +
      scale_color_manual(name = "CDF", values = c("Empírica (Amostra)" = "darkgreen", "Analítica (Teórica)" = "red")) +
      labs(x = NULL, y = "Prob. Acumulada") +
      theme_minimal() +
      theme(legend.position = "none")
    
    plot_pair <- g1 + g2
    row_plot_pairs[[as.character(n_val)]] <- plot_pair
  }
  
  full_row_layout <- wrap_plots(c(list(g_row_title), row_plot_pairs), 
                                nrow = 1, 
                                widths = c(2.5, rep(5, length(n_valores))))
  plot_rows[[i]] <- full_row_layout
}

final_plot <- wrap_plots(plot_rows, ncol = 1)
final_plot_with_legend <- final_plot + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom', legend.box = 'vertical')

print(final_plot_with_legend)

# ----- Testes de Aderência -----

cat("### Tabela 1: Primeiros 50 Valores Pseudoaleatórios Gerados\n\n")



for (i in 1:nrow(params_unicos)) {
  a <- params_unicos$alpha[i]
  b <- params_unicos$beta[i]
  r <- params_unicos$rho[i]
  
  titulo_tabela <- paste0("Parâmetros: alpha = ", a, ", beta = ", b, ", rho = ", r)
  lista_amostras <- list()
  
  for (n_val in n_valores) {
    set.seed(42)
    amostra_head <- head(rcompound_weibullN(n_val, r, a, b), 50)
    lista_amostras[[paste0("n=", n_val)]] <- c(amostra_head, rep(NA, 50 - length(amostra_head)))
  }
  
  df_amostras <- as.data.frame(lista_amostras)
  
  print(
    kable(df_amostras, "html", caption = titulo_tabela, digits = 4) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)
  )
  cat("\n\n")
}

cat("### Tabela 2: Resultados dos Testes de Aderência\n\n")

resultados_testes <- data.frame()
for (i in 1:nrow(parametros)) {
  n_atual     <- parametros$n[i]
  alpha_atual <- parametros$alpha[i]
  beta_atual  <- parametros$beta[i]
  rho_atual   <- parametros$rho[i]
  set.seed(123)
  amostra <- rcompound_weibullN(n_atual, rho_atual, alpha_atual, beta_atual)
  amostra <- na.omit(amostra) 
  
  if(length(amostra) < 2) next

  ks_result <- ks.test(amostra, pcompound_weibull, rho = rho_atual, alpha = alpha_atual, beta = beta_atual)
  ks_formatado <- sprintf("D = %.4f, p-valor = %.4f", ks_result$statistic, ks_result$p.value)
  
  cvm_result <- cvm.test(amostra, pcompound_weibull, rho = rho_atual, alpha = alpha_atual, beta = beta_atual)
  cvm_formatado <- sprintf("ω² = %.4f, p-valor = %.4f", cvm_result$statistic, cvm_result$p.value)
  
  resultados_testes <- rbind(resultados_testes, data.frame(
    alpha = alpha_atual,
    beta = beta_atual,
    rho = rho_atual,
    n = n_atual,
    `Kolmogorov-Smirnov` = ks_formatado,
    `Cramér-von Mises` = cvm_formatado
  ))
}

print(
  kable(resultados_testes, "html", caption = "Resultados dos Testes de Aderência vs. Distribuição Compound Weibull Teórica", booktabs = TRUE) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE) %>%
    column_spec(1:4, bold = TRUE)
)
cat("\n\n")

# --------------------------------------------------------------------------------------
# Tarefa 4: Eficiência computacional
# --------------------------------------------------------------------------------------
# - Tabela com os tempos de execução da geração para cada cenário.
# - Gráfico ilustrando os resultados.
# --------------------------------------------------------------------------------------
cat("### 4. Análise de Eficiência Computacional\n\n")

resultados_tempo <- data.frame()

for (i in 1:nrow(parametros)) {
  n_atual     <- parametros$n[i]
  alpha_atual <- parametros$alpha[i]
  beta_atual  <- parametros$beta[i]
  rho_atual   <- parametros$rho[i]
  
  tempo_decorrido <- system.time({
    suppressWarnings({
      amostra <- rcompound_weibullN(n_atual, rho_atual, alpha_atual, beta_atual)
    })
  })["elapsed"]
  
  resultados_tempo <- rbind(resultados_tempo, data.frame(
    alpha = alpha_atual,
    beta = beta_atual,
    rho = rho_atual,
    n = n_atual,
    tempo_execucao_s = as.numeric(tempo_decorrido)
  ))
}

cat("#### Tabela 3: Tempo de Execução para Geração das Amostras\n\n")
print(
  kable(resultados_tempo, "html", 
        caption = "Tempo de Execução (em segundos) por Cenário", 
        digits = 4,
        col.names = c("alpha", "beta", "rho", "Tamanho da Amostra (n)", "Tempo (s)")) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE) %>%
    column_spec(4:5, bold = TRUE)
)
cat("\n\n")

cat("#### Gráfico: Relação entre Tamanho da Amostra e Tempo de Execução\n\n")

resultados_tempo$parametros_legenda <- factor(paste("beta =", resultados_tempo$beta, ", rho =", resultados_tempo$rho))

grafico_tempo <- ggplot(resultados_tempo, aes(x = n, y = tempo_execucao_s, color = parametros_legenda)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_log10(
    breaks = c(10, 100, 1000, 10000, 100000),
    labels = label_comma()
  ) +
  scale_y_continuous(labels = label_number(suffix = " s")) +
  labs(
    title = "Tempo de Execução para Geração de Amostras",
    subtitle = "Método da Transformação Inversa Numérica",
    x = "Tamanho da Amostra (n) - Escala Logarítmica",
    y = "Tempo de Execução (segundos)",
    color = "Combinação de Parâmetros"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(grafico_tempo)
