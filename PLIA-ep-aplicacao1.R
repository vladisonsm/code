library(statmod)
library(PLreg)
library(numDeriv)

dados = Planilha.sem.título...OphthalmologicalStudyData
dados$y = as.numeric(gsub(",", ".", dados$y))

dados$concentracao <- factor(dados$concentracao)
dados$concentracao <- relevel(dados$concentracao, ref = "-1") 

dados$log_t <- log(dados$tempo)
dados$GAS20 <- as.numeric(dados$concentracao == "0")
dados$GAS25 <- as.numeric(dados$concentracao == "1")

id <- as.numeric(factor(dados$ind))
y  <- dados$y

X <- cbind(1, dados$log_t, dados$GAS20, dados$GAS25,
           dados$GAS20 * dados$log_t, dados$GAS25 * dados$log_t)
colnames(X) <- c("Intercepto", "log_t", "GAS20", "GAS25", "GAS20:log_t", "GAS25:log_t")
p <- ncol(X)

# ==============================================================================
# 1. FUNÇÕES DE LOG-VEROSSIMILHANÇA (ADAPTADAS PARA EXPONENCIAL POTÊNCIA)
# ==============================================================================
# Modelo com 1 efeito aleatório (média)
log_vero_PLIA_vet <- function(theta, y, X, S, id, nos, pesos, k) {
  p <- ncol(X); q <- ncol(S)
  beta   <- theta[1:p]
  tau    <- theta[(p + 1):(p + q)]
  lambda <- exp(theta[p + q + 1])
  sigma2 <- exp(theta[p + q + 2])
  
  # Calcula o Zeta equivalente ao k escolhido[cite: 2]
  zeta_val <- 2 / (1 + k) 
  
  Q <- length(nos)
  b_0_vet <- sqrt(2 * sigma2) * nos
  pesos_prod <- pesos / sqrt(pi)
  
  ll_total <- 0
  for (i in unique(id)) {
    idx <- (id == i)
    y_i <- y[idx]; X_i <- X[idx, , drop = FALSE]; S_i <- S[idx, , drop = FALSE]
    n_i <- length(y_i)
    
    eta1_mat <- outer(as.vector(X_i %*% beta), b_0_vet, "+")
    eta2_mat <- matrix(rep(as.vector(S_i %*% tau), Q), nrow = n_i, ncol = Q)
    
    mu_mat    <- pmin(pmax(plogis(eta1_mat), 1e-6), 1 - 1e-6)
    sigma_mat <- pmin(pmax(exp(eta2_mat), 1e-4), 50)
    y_mat     <- matrix(rep(y_i, Q), nrow = n_i, ncol = Q)
    
    # Utilizando densidade PE e zeta_val[cite: 2]
    pdf_vec <- suppressWarnings(PLreg::dPL(as.vector(y_mat), as.vector(mu_mat),
                                           as.vector(sigma_mat), lambda, 
                                           zeta = zeta_val, family = "PE"))
    if (any(is.na(pdf_vec)) || any(is.infinite(pdf_vec))) return(1e10)
    
    pdf_mat <- matrix(pdf_vec, nrow = n_i, ncol = Q)
    log_pdf_mat <- log(pmax(pdf_mat, .Machine$double.xmin))
    log_densidade_conjunta <- colSums(log_pdf_mat)
    
    termos_soma <- log_densidade_conjunta + log(pesos_prod)
    max_termo <- max(termos_soma)
    log_soma_quadratura <- max_termo + log(sum(exp(termos_soma - max_termo)))
    ll_total <- ll_total + log_soma_quadratura
  }
  if (is.na(ll_total) || is.infinite(ll_total)) return(1e10)
  return(-ll_total)
}

# Modelo com 2 efeitos aleatórios (média e dispersão)
log_vero_PLIA_2RE <- function(theta, y, X, S, id, nos, pesos, k) {
  p <- ncol(X); q <- ncol(S)
  beta   <- theta[1:p]
  tau    <- theta[(p + 1):(p + q)]
  lambda <- exp(theta[p + q + 1])
  sigma2 <- exp(theta[p + q + 2])
  ups2   <- exp(theta[p + q + 3])
  
  # Calcula o Zeta equivalente ao k escolhido[cite: 2]
  zeta_val <- 2 / (1 + k) 
  
  b0_vet <- sqrt(2 * sigma2) * nos
  d0_vet <- sqrt(2 * ups2)  * nos
  
  grade  <- expand.grid(b0 = b0_vet, d0 = d0_vet)
  pesos2 <- expand.grid(wb = pesos, wd = pesos)
  pesos_prod <- (pesos2$wb * pesos2$wd) / pi
  K <- nrow(grade)
  
  ll_total <- 0
  for (i in unique(id)) {
    idx <- (id == i)
    y_i <- y[idx]; X_i <- X[idx, , drop = FALSE]; S_i <- S[idx, , drop = FALSE]
    n_i <- length(y_i)
    
    xb_i <- as.vector(X_i %*% beta)
    s_i  <- as.vector(S_i %*% tau)
    
    eta1_mat <- outer(xb_i, grade$b0, "+")
    eta2_mat <- outer(s_i,  grade$d0, "+")
    
    mu_mat    <- pmin(pmax(plogis(eta1_mat), 1e-6), 1 - 1e-6)
    sigma_mat <- pmin(pmax(exp(eta2_mat), 1e-4), 50)
    y_mat     <- matrix(rep(y_i, K), nrow = n_i, ncol = K)
    
    # Utilizando densidade PE e zeta_val[cite: 2]
    pdf_vec <- suppressWarnings(PLreg::dPL(as.vector(y_mat), as.vector(mu_mat),
                                           as.vector(sigma_mat), lambda, 
                                           zeta = zeta_val, family = "PE"))
    if (any(is.na(pdf_vec)) || any(is.infinite(pdf_vec))) return(1e10)
    
    pdf_mat <- matrix(pdf_vec, nrow = n_i, ncol = K)
    log_pdf_mat <- log(pmax(pdf_mat, .Machine$double.xmin))
    log_densidade_conjunta <- colSums(log_pdf_mat)
    
    termos_soma <- log_densidade_conjunta + log(pesos_prod)
    max_termo <- max(termos_soma)
    log_soma_quadratura <- max_termo + log(sum(exp(termos_soma - max_termo)))
    ll_total <- ll_total + log_soma_quadratura
  }
  if (is.na(ll_total) || is.infinite(ll_total)) return(1e10)
  return(-ll_total)
}

# Modelo nulo (sem efeitos aleatórios)
log_vero_PL_fixo <- function(theta, y, X, S, k) {
  p <- ncol(X); q <- ncol(S)
  beta   <- theta[1:p]
  tau    <- theta[(p + 1):(p + q)]
  lambda <- exp(theta[p + q + 1])
  
  zeta_val <- 2 / (1 + k) #[cite: 2]
  
  eta1 <- as.vector(X %*% beta)
  eta2 <- as.vector(S %*% tau)
  mu   <- pmin(pmax(plogis(eta1), 1e-6), 1 - 1e-6)
  sigma <- pmin(pmax(exp(eta2), 1e-4), 50)
  
  pdf_vec <- suppressWarnings(PLreg::dPL(y, mu, sigma, lambda, zeta = zeta_val, family = "PE"))
  if (any(is.na(pdf_vec)) || any(is.infinite(pdf_vec))) return(1e10)
  ll <- sum(log(pmax(pdf_vec, .Machine$double.xmin)))
  if (is.na(ll) || is.infinite(ll)) return(1e10)
  return(-ll)
}

# Modelo com 1 efeito aleatório (apenas dispersão)
log_vero_PL_1RE_disp <- function(theta, y, X, S, id, nos, pesos, k) {
  p <- ncol(X); q <- ncol(S)
  beta   <- theta[1:p]
  tau    <- theta[(p + 1):(p + q)]
  lambda <- exp(theta[p + q + 1])
  ups2   <- exp(theta[p + q + 2]) 
  
  zeta_val <- 2 / (1 + k) #[cite: 2]
  
  d0_vet <- sqrt(2 * ups2) * nos
  pesos_prod <- pesos / sqrt(pi)
  K <- length(nos)
  
  ll_total <- 0
  for (i in unique(id)) {
    idx <- (id == i)
    y_i <- y[idx]; X_i <- X[idx, , drop = FALSE]; S_i <- S[idx, , drop = FALSE]
    n_i <- length(y_i)
    
    eta1_mat <- matrix(rep(as.vector(X_i %*% beta), K), nrow = n_i, ncol = K)
    eta2_mat <- outer(as.vector(S_i %*% tau), d0_vet, "+")
    
    mu_mat    <- pmin(pmax(plogis(eta1_mat), 1e-6), 1 - 1e-6)
    sigma_mat <- pmin(pmax(exp(eta2_mat), 1e-4), 50)
    y_mat     <- matrix(rep(y_i, K), nrow = n_i, ncol = K)
    
    pdf_vec <- suppressWarnings(PLreg::dPL(as.vector(y_mat), as.vector(mu_mat),
                                           as.vector(sigma_mat), lambda, 
                                           zeta = zeta_val, family = "PE"))
    if (any(is.na(pdf_vec)) || any(is.infinite(pdf_vec))) return(1e10)
    
    pdf_mat <- matrix(pdf_vec, nrow = n_i, ncol = K)
    log_pdf_mat <- log(pmax(pdf_mat, .Machine$double.xmin))
    log_densidade_conjunta <- colSums(log_pdf_mat)
    
    termos_soma <- log_densidade_conjunta + log(pesos_prod)
    max_termo <- max(termos_soma)
    log_soma_quadratura <- max_termo + log(sum(exp(termos_soma - max_termo)))
    ll_total <- ll_total + log_soma_quadratura
  }
  if (is.na(ll_total) || is.infinite(ll_total)) return(1e10)
  return(-ll_total)
}

# ==============================================================================
# 2. CHUTES INICIAIS
# ==============================================================================
# (As funções de chutes iniciais permanecem as mesmas)
obter_chutes_iniciais <- function(y, X, S) {
  y_adj <- pmin(pmax(y, 1e-5), 1 - 1e-5)
  z_resp <- qlogis(y_adj)
  fit_lm <- lm(z_resp ~ X - 1)
  beta_init <- as.vector(coef(fit_lm))
  tau_init  <- c(log(sd(residuals(fit_lm))^2), rep(0, ncol(S) - 1))
  log_lambda_init <- log(1.0)
  log_sigma2_init <- log(0.1)
  c(beta_init, tau_init, log_lambda_init, log_sigma2_init)
}

obter_chutes_iniciais_2RE <- function(y, X, S) {
  base <- obter_chutes_iniciais(y, X, S)
  c(base, log(0.1))
}

obter_chutes_iniciais_fixo <- function(y, X, S) {
  base <- obter_chutes_iniciais(y, X, S)
  base[-length(base)]
}

obter_chutes_iniciais_1RE_disp <- function(y, X, S) {
  base <- obter_chutes_iniciais(y, X, S)
  c(base[1:(length(base)-1)], log(0.1))
}

# ==============================================================================
# 3. ROTINA DE VEROSSIMILHANÇA PERFILADA (BUSCA DO MELHOR k)
# ==============================================================================
# Malha de valores para k, onde -1 < k <= 1
grade_k <- c(-0.9, -0.7, -0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0)

ajustar_perfil_k <- function(obj_func, theta_init, lower_b, upper_b, 
                             y, X, S, id=NULL, nos=NULL, pesos=NULL, req_quad = TRUE) {
  melhor_obj <- Inf
  melhor_fit <- NULL
  melhor_k  <- NA
  
  for (k_val in grade_k) {
    if (req_quad) {
      fit_temp <- tryCatch({
        nlminb(start = theta_init, objective = obj_func,
               y = y, X = X, S = S, id = id, nos = nos, pesos = pesos, k = k_val,
               lower = lower_b, upper = upper_b,
               control = list(eval.max = 3000, iter.max = 3000, rel.tol = 1e-6))
      }, error = function(e) NULL)
    } else {
      # Para o modelo nulo que não usa quadratura
      fit_temp <- tryCatch({
        nlminb(start = theta_init, objective = obj_func,
               y = y, X = X, S = S, k = k_val,
               lower = lower_b, upper = upper_b,
               control = list(eval.max = 3000, iter.max = 3000, rel.tol = 1e-6))
      }, error = function(e) NULL)
    }
    
    if (!is.null(fit_temp) && fit_temp$convergence == 0) {
      if (fit_temp$objective < melhor_obj) {
        melhor_obj <- fit_temp$objective
        melhor_fit <- fit_temp
        melhor_k  <- k_val
      }
    }
  }
  return(list(fit = melhor_fit, k = melhor_k))
}

# ==============================================================================
# 4. TABELAS E TESTES (AUXILIARES)
# ==============================================================================
# (Permanecem inalteradas)
teste_lrt_fronteira <- function(logLik_completo, logLik_nulo, df_extra = 1) {
  LR <- 2 * (logLik_completo - logLik_nulo)
  LR <- max(LR, 0)
  p_valor <- 0.5 * (1 - pchisq(LR, df = df_extra))
  data.frame(LR = LR, p_valor = p_valor)
}

construir_tabela <- function(fit, hess, nomes, idx_lambda = NULL, idx_var = NULL) {
  vcov_par <- tryCatch(solve(hess), error = function(e) {
    warning("Hessiana singular.")
    matrix(NA, nrow(hess), ncol(hess))
  })
  se_par <- sqrt(diag(vcov_par))
  estim  <- fit$par
  
  Estimativa <- estim; EP <- se_par; Z <- estim / se_par; p_valor <- 2 * (1 - pnorm(abs(Z)))
  
  for (j in c(idx_lambda, idx_var)) {
    valor_original <- exp(estim[j])
    Estimativa[j]  <- valor_original
    EP[j]          <- valor_original * se_par[j]
    nomes[j]       <- sub("^log_", "", nomes[j])
  }
  for (j in idx_lambda) {
    Z[j]       <- (Estimativa[j] - 1) / EP[j]
    p_valor[j] <- 2 * (1 - pnorm(abs(Z[j])))
  }
  for (j in idx_var) { Z[j] <- NA; p_valor[j] <- NA }
  
  data.frame(Parametro = nomes, Estimativa = Estimativa, EP = EP, Z = Z, p_valor = p_valor)
}

# ==============================================================================
# 5. AJUSTE DOS MODELOS COM PERFILAMENTO DE K
# ==============================================================================
gh  <- statmod::gauss.quad(20, kind = "hermite")
gh2 <- statmod::gauss.quad(20, kind = "hermite") # Para integral dupla

# --- Modelo 1 ---
S_1 <- matrix(1, nrow = length(y), ncol = 1); colnames(S_1) <- "Intercepto"
q1 <- ncol(S_1)
res_1 <- ajustar_perfil_k(
  obj_func = log_vero_PLIA_vet, theta_init = obter_chutes_iniciais(y, X, S_1),
  lower_b = c(rep(-10, p), rep(-10, q1), -3, -5), upper_b = c(rep(10, p), rep(10, q1), 3, 5),
  y = y, X = X, S = S_1, id = id, nos = gh$nodes, pesos = gh$weights
)
fit_1 <- res_1$fit; k_1 <- res_1$k
hess_1 <- numDeriv::hessian(log_vero_PLIA_vet, x = fit_1$par, y = y, X = X, S = S_1, 
                            id = id, nos = gh$nodes, pesos = gh$weights, k = k_1)
tabela_1 <- construir_tabela(fit_1, hess_1, c(colnames(X), "tau0", "log_lambda", "log_sigma2_b"), 
                             idx_lambda = p + q1 + 1, idx_var = p + q1 + 2)

# --- Modelo 2 ---
S_2 <- cbind(1, dados$log_t); colnames(S_2) <- c("Intercepto", "log_t")
q2 <- ncol(S_2)
res_2 <- ajustar_perfil_k(
  obj_func = log_vero_PLIA_vet, theta_init = obter_chutes_iniciais(y, X, S_2),
  lower_b = c(rep(-10, p), rep(-10, q2), -3, -5), upper_b = c(rep(10, p), rep(10, q2), 3, 5),
  y = y, X = X, S = S_2, id = id, nos = gh$nodes, pesos = gh$weights
)
fit_2 <- res_2$fit; k_2 <- res_2$k
hess_2 <- numDeriv::hessian(log_vero_PLIA_vet, x = fit_2$par, y = y, X = X, S = S_2, 
                            id = id, nos = gh$nodes, pesos = gh$weights, k = k_2)
tabela_2 <- construir_tabela(fit_2, hess_2, c(colnames(X), "tau0", "tau1_logt", "log_lambda", "log_sigma2_b"), 
                             idx_lambda = p + q2 + 1, idx_var = p + q2 + 2)

# --- Modelo 3 ---
S_3 <- matrix(1, nrow = length(y), ncol = 1); colnames(S_3) <- "Intercepto"
q3 <- ncol(S_3)
res_3 <- ajustar_perfil_k(
  obj_func = log_vero_PLIA_2RE, theta_init = obter_chutes_iniciais_2RE(y, X, S_3),
  lower_b = c(rep(-10, p), rep(-10, q3), -3, -5, -5), upper_b = c(rep(10, p), rep(10, q3), 3, 5, 5),
  y = y, X = X, S = S_3, id = id, nos = gh2$nodes, pesos = gh2$weights
)
fit_3 <- res_3$fit; k_3 <- res_3$k
hess_3 <- numDeriv::hessian(log_vero_PLIA_2RE, x = fit_3$par, y = y, X = X, S = S_3, 
                            id = id, nos = gh2$nodes, pesos = gh2$weights, k = k_3)
tabela_3 <- construir_tabela(fit_3, hess_3, c(colnames(X), "tau0", "log_lambda", "log_sigma2_b", "log_upsilon2_d"), 
                             idx_lambda = p + q3 + 1, idx_var = c(p + q3 + 2, p + q3 + 3))

# ==============================================================================
# 6. MODELOS NULOS PARA LRT (FIXADOS NO K ESTIMADO PELOS RESPECTIVOS MODELOS)
# ==============================================================================
theta_init_1_nulo <- obter_chutes_iniciais_fixo(y, X, S_1)
fit_1_nulo <- nlminb(start = theta_init_1_nulo, objective = log_vero_PL_fixo,
                     y = y, X = X, S = S_1, k = k_1, lower = c(rep(-10, p), rep(-10, q1), -3),
                     upper = c(rep(10, p), rep(10, q1), 3), control = list(eval.max = 3000, iter.max = 3000))

theta_init_2_nulo <- obter_chutes_iniciais_fixo(y, X, S_2)
fit_2_nulo <- nlminb(start = theta_init_2_nulo, objective = log_vero_PL_fixo,
                     y = y, X = X, S = S_2, k = k_2, lower = c(rep(-10, p), rep(-10, q2), -3),
                     upper = c(rep(10, p), rep(10, q2), 3), control = list(eval.max = 3000, iter.max = 3000))

theta_init_3_nulo <- obter_chutes_iniciais_1RE_disp(y, X, S_3)
fit_3_nulo_sigma2b <- nlminb(start = theta_init_3_nulo, objective = log_vero_PL_1RE_disp,
                             y = y, X = X, S = S_3, id = id, nos = gh$nodes, pesos = gh$weights, k = k_3,
                             lower = c(rep(-10, p), rep(-10, q3), -3, -5), upper = c(rep(10, p), rep(10, q3), 3, 5),
                             control = list(eval.max = 3000, iter.max = 3000))

# Testes de Fronteira
lrt_sigma2b_1   <- teste_lrt_fronteira(-fit_1$objective, -fit_1_nulo$objective)
lrt_sigma2b_2   <- teste_lrt_fronteira(-fit_2$objective, -fit_2_nulo$objective)
lrt_upsilon2d_3 <- teste_lrt_fronteira(-fit_3$objective, -fit_1$objective) 
lrt_sigma2b_3   <- teste_lrt_fronteira(-fit_3$objective, -fit_3_nulo_sigma2b$objective)

tabela_lrt <- data.frame(
  Teste = c("H0: sigma2_b = 0  (Modelo 1)", "H0: sigma2_b = 0  (Modelo 2)",
            "H0: upsilon2_d = 0 (Modelo 3 | b0i)", "H0: sigma2_b = 0 (Modelo 3 | d0i)"),
  LR = c(lrt_sigma2b_1$LR, lrt_sigma2b_2$LR, lrt_upsilon2d_3$LR, lrt_sigma2b_3$LR),
  p_valor = c(lrt_sigma2b_1$p_valor, lrt_sigma2b_2$p_valor, lrt_upsilon2d_3$p_valor, lrt_sigma2b_3$p_valor)
)

tabela_1$p_valor[tabela_1$Parametro == "sigma2_b"] <- lrt_sigma2b_1$p_valor
tabela_2$p_valor[tabela_2$Parametro == "sigma2_b"] <- lrt_sigma2b_2$p_valor
tabela_3$p_valor[tabela_3$Parametro == "sigma2_b"] <- lrt_sigma2b_3$p_valor
tabela_3$p_valor[tabela_3$Parametro == "upsilon2_d"] <- lrt_upsilon2d_3$p_valor

# ==============================================================================
# 7. RESULTADOS FINAIS
# ==============================================================================
cat(sprintf("\n--- Modelo 1 (PE, k estimado: %.2f) ---\n", k_1)); print(tabela_1, row.names = FALSE)
cat(sprintf("\n--- Modelo 2 (PE, k estimado: %.2f) ---\n", k_2)); print(tabela_2, row.names = FALSE)
cat(sprintf("\n--- Modelo 3 (PE, k estimado: %.2f) ---\n", k_3)); print(tabela_3, row.names = FALSE)

cat("\n--- TRV corrigido para os efeitos aleatórios ---\n")
print(tabela_lrt, row.names = FALSE)

comparar_modelos <- function(fits, k_vals, nomes_modelos, n) {
  data.frame(
    Modelo       = nomes_modelos,
    k_est        = k_vals,
    npar         = sapply(fits, function(f) length(f$par) + 1), # +1 considerando a estimação do k
    logLik       = sapply(fits, function(f) -f$objective),
    AIC          = sapply(fits, function(f) 2 * (length(f$par) + 1) + 2 * f$objective),
    BIC          = sapply(fits, function(f) (length(f$par) + 1) * log(n) + 2 * f$objective),
    CAIC         = sapply(fits, function(f) (length(f$par) + 1) * log(n + 1) + 2 * f$objective),
    HQIC         = sapply(fits, function(f) (length(f$par) + 1) * log(log(n)) + 2 * f$objective),
    Convergencia = sapply(fits, function(f) f$convergence)
  )
}
tabela_comparativa <- comparar_modelos(list(fit_1, fit_2, fit_3), c(k_1, k_2, k_3),
                                       c("Modelo 1", "Modelo 2", "Modelo 3"), n = length(y))
cat("\n--- Comparação dos modelos (Exponencial Potência) ---\n"); print(tabela_comparativa, row.names = FALSE)

# ==============================================================================
# FUNÇÕES DE RESÍDUOS CORRIGIDAS (SEM ANOTAÇÕES DE CITAÇÃO)
# ==============================================================================
# Função para Modelos com 1 Efeito Aleatório (Modelo 1 e Modelo 2)
calcular_residuos_1RE <- function(fit, y, X, S, id, gh, k) {
  p <- ncol(X); q <- ncol(S)
  beta   <- fit$par[1:p]
  tau    <- fit$par[(p + 1):(p + q)]
  lambda <- exp(fit$par[p + q + 1])
  sigma2 <- exp(fit$par[p + q + 2])
  
  zeta_val <- 2 / (1 + k)
  
  nos <- gh$nodes; pesos <- gh$weights
  Q <- length(nos)
  b_0_vet <- sqrt(2 * sigma2) * nos
  pesos_prod <- pesos / sqrt(pi)
  
  ids_unicos <- unique(id)
  b0_ebe <- numeric(length(ids_unicos))
  names(b0_ebe) <- ids_unicos
  
  for (i in ids_unicos) {
    idx <- (id == i)
    y_i <- y[idx]; X_i <- X[idx, , drop = FALSE]; S_i <- S[idx, , drop = FALSE]
    n_i <- length(y_i)
    
    eta1_mat <- outer(as.vector(X_i %*% beta), b_0_vet, "+")
    eta2_mat <- matrix(rep(as.vector(S_i %*% tau), Q), nrow = n_i, ncol = Q)
    
    mu_mat    <- pmin(pmax(plogis(eta1_mat), 1e-6), 1 - 1e-6)
    sigma_mat <- pmin(pmax(exp(eta2_mat), 1e-4), 50)
    y_mat     <- matrix(rep(y_i, Q), nrow = n_i, ncol = Q)
    
    pdf_vec <- suppressWarnings(PLreg::dPL(as.vector(y_mat), as.vector(mu_mat),
                                           as.vector(sigma_mat), lambda, 
                                           zeta = zeta_val, family = "PE"))
    pdf_mat <- matrix(pdf_vec, nrow = n_i, ncol = Q)
    log_pdf_mat <- log(pmax(pdf_mat, .Machine$double.xmin))
    log_dens_conjunta <- colSums(log_pdf_mat)
    
    termos <- log_dens_conjunta + log(pesos_prod)
    w_post <- exp(termos - max(termos)) / sum(exp(termos - max(termos)))
    
    b0_ebe[as.character(i)] <- sum(b_0_vet * w_post)
  }
  
  b0_i_vet <- b0_ebe[as.character(id)]
  eta1_cond <- as.vector(X %*% beta) + b0_i_vet
  eta2_cond <- as.vector(S %*% tau)
  
  mu_cond    <- pmin(pmax(plogis(eta1_cond), 1e-6), 1 - 1e-6)
  sigma_cond <- pmin(pmax(exp(eta2_cond), 1e-4), 50)
  
  u_cond <- suppressWarnings(PLreg::pPL(y, mu_cond, sigma_cond, lambda, 
                                        zeta = zeta_val, family = "PE"))
  u_cond <- pmin(pmax(u_cond, 1e-7), 1 - 1e-7)
  
  r_cond <- qnorm(u_cond)
  return(list(r_cond = r_cond, b0_ebe = b0_ebe))
}

# Função para Modelo com 2 Efeitos Aleatórios (Modelo 3)
calcular_residuos_2RE <- function(fit, y, X, S, id, gh2, k) {
  p <- ncol(X); q <- ncol(S)
  beta   <- fit$par[1:p]
  tau    <- fit$par[(p + 1):(p + q)]
  lambda <- exp(fit$par[p + q + 1])
  sigma2 <- exp(fit$par[p + q + 2])
  ups2   <- exp(fit$par[p + q + 3])
  
  zeta_val <- 2 / (1 + k)
  
  nos <- gh2$nodes; pesos <- gh2$weights
  b0_vet <- sqrt(2 * sigma2) * nos
  d0_vet <- sqrt(2 * ups2) * nos
  
  grade      <- expand.grid(b0 = b0_vet, d0 = d0_vet)
  pesos2     <- expand.grid(wb = pesos, wd = pesos)
  pesos_prod <- (pesos2$wb * pesos2$wd) / pi
  K <- nrow(grade)
  
  ids_unicos <- unique(id)
  b0_ebe <- numeric(length(ids_unicos)); names(b0_ebe) <- ids_unicos
  d0_ebe <- numeric(length(ids_unicos)); names(d0_ebe) <- ids_unicos
  
  for (i in ids_unicos) {
    idx <- (id == i)
    y_i <- y[idx]; X_i <- X[idx, , drop = FALSE]; S_i <- S[idx, , drop = FALSE]
    n_i <- length(y_i)
    
    eta1_mat <- outer(as.vector(X_i %*% beta), grade$b0, "+")
    eta2_mat <- outer(as.vector(S_i %*% tau),  grade$d0, "+")
    
    mu_mat    <- pmin(pmax(plogis(eta1_mat), 1e-6), 1 - 1e-6)
    sigma_mat <- pmin(pmax(exp(eta2_mat), 1e-4), 50)
    y_mat     <- matrix(rep(y_i, K), nrow = n_i, ncol = K)
    
    pdf_vec <- suppressWarnings(PLreg::dPL(as.vector(y_mat), as.vector(mu_mat),
                                           as.vector(sigma_mat), lambda, 
                                           zeta = zeta_val, family = "PE"))
    pdf_mat <- matrix(pdf_vec, nrow = n_i, ncol = K)
    log_pdf_mat <- log(pmax(pdf_mat, .Machine$double.xmin))
    log_dens_conjunta <- colSums(log_pdf_mat)
    
    termos <- log_dens_conjunta + log(pesos_prod)
    w_post <- exp(termos - max(termos)) / sum(exp(termos - max(termos)))
    
    b0_ebe[as.character(i)] <- sum(grade$b0 * w_post)
    d0_ebe[as.character(i)] <- sum(grade$d0 * w_post)
  }
  
  b0_i_vet <- b0_ebe[as.character(id)]
  d0_i_vet <- d0_ebe[as.character(id)]
  
  eta1_cond <- as.vector(X %*% beta) + b0_i_vet
  eta2_cond <- as.vector(S %*% tau)  + d0_i_vet
  
  mu_cond    <- pmin(pmax(plogis(eta1_cond), 1e-6), 1 - 1e-6)
  sigma_cond <- pmin(pmax(exp(eta2_cond), 1e-4), 50)
  
  u_cond <- suppressWarnings(PLreg::pPL(y, mu_cond, sigma_cond, lambda, 
                                        zeta = zeta_val, family = "PE"))
  u_cond <- pmin(pmax(u_cond, 1e-7), 1 - 1e-7)
  
  r_cond <- qnorm(u_cond)
  return(list(r_cond = r_cond, b0_ebe = b0_ebe, d0_ebe = d0_ebe))
}

# ==============================================================================
# 2. CÁLCULO DIAGNÓSTICO E GERAÇÃO DOS GRÁFICOS
# ==============================================================================
par(mfrow = c(1, 3))
res_obj <- calcular_residuos_2RE(fit = fit_3, y = y, X = X, S = S_3, id = id, gh = gh, k = k_3)
r_cond  <- res_obj$r_cond
obs_id  <- 1:length(y)

# --- Gráfico 1: Resíduos por Paciente ---
plot(id, r_cond,
     pch = 3,                                      # Símbolo '+'
     main = "Modelo 3",
     cex.main = 1.7,
     xlab = "Índices dos pacientes",               # Rótulo eixo X
     ylab = "",     # Rótulo eixo Y[cite: 2]
     cex.lab = 1.4,
     las = 1,
     ylim=c(-3.5,3.5))

# Linhas limite inferior e superior de referência (-2.5 e 2.5)[cite: 2]
limite <- 2.5
abline(h = c(-limite, limite), lty = 3)

# Identificar observações fora dos limites [-2.5, 2.5][cite: 2]
idx_fora <- which(abs(r_cond) > limite)

if (length(idx_fora) > 0) {
  # Rótulo no formato 'ID_PACIENTE.NUM_OBSERVAÇÃO' (ex: 25.144)[cite: 2]
  rotulos <- paste0(id[idx_fora], ".", obs_id[idx_fora])
  
  text(x = id[idx_fora], 
       y = r_cond[idx_fora], 
       labels = rotulos, 
       pos = 2,                                    
       cex = 1.0)
  
  tabela_outliers <- data.frame(
    ID_Paciente      = id[idx_fora],
    Linha_Observacao = obs_id[idx_fora],
    Rotulo_Grafico   = rotulos,
    Residuo          = round(r_cond[idx_fora], 4)
  )
  cat("\n--- Observações fora dos limites (-2.5 ou +2.5) ---\n")
  print(tabela_outliers)
} else {
  cat("\nNenhuma observação ficou fora dos limites especificados.\n")
}

# --- Gráfico 2: QQ-Plot com Envelope Simulado ---
qqplot_envelope <- function(r_cond, B = 100, conf = 0.95, seed = 123) {
  set.seed(seed)
  
  # Remover valores ausentes caso existam
  r_clean <- r_cond[!is.na(r_cond)]
  n <- length(r_clean)
  
  # 1. Quantis teóricos da Normal Padrão
  q_teorico <- qnorm(ppoints(n))
  
  # 2. Ordenar os resíduos observados
  r_obs_ord <- sort(r_clean)
  
  # 3. Gerar B simulações sob a hipótese nula (N(0,1))
  sim_mat <- matrix(0, nrow = n, ncol = B)
  for (b in 1:B) {
    sim_mat[, b] <- sort(rnorm(n, mean = 0, sd = 1))
  }
  
  # 4. Calcular os limites do envelope (quantis 2.5%, média e 97.5%)
  alpha <- (1 - conf) / 2
  env_inf <- apply(sim_mat, 1, quantile, probs = alpha)
  env_med <- apply(sim_mat, 1, mean)
  env_sup <- apply(sim_mat, 1, quantile, probs = 1 - alpha)
  
  # 5. Calcular a porcentagem de pontos fora do envelope
  pontos_fora <- (r_obs_ord < env_inf) | (r_obs_ord > env_sup)
  perc_fora <- (sum(pontos_fora) / n) * 100
  
  # Limites do eixo Y
  ylim_val <- range(c(r_obs_ord, env_inf, env_sup), na.rm = TRUE)
  
  # 6. Criar o gráfico base
  plot(q_teorico, r_obs_ord,
       pch = 3,                                      # Símbolo '+'
       xlab = "Quantis teóricos",                    # Rótulo eixo X
       ylab = "",   # Rótulo eixo Y
       ylim = ylim_val,
       main = 'Modelo 3',
       cex.lab = 1.4,
       cex.main = 1.7,
       las = 1)
  
  # 7. Desenhar as linhas dos envelopes e da média esperada
  lines(q_teorico, env_inf, lty = 1, lwd = 1.2)      # Envelope inferior
  lines(q_teorico, env_med, lty = 2, lwd = 1.2)      # Linha de referência central
  lines(q_teorico, env_sup, lty = 1, lwd = 1.2)      # Envelope superior
  
  # 8. Adicionar a porcentagem no gráfico (canto superior esquerdo)
  texto_legenda <- sprintf("Fora do envelope: %.1f%%", perc_fora)
  legend("topleft", legend = texto_legenda, bty = "n", cex = 1.5)
  
  # Opcional: Retornar os índices dos resíduos que ficaram fora para análise extra
  return(invisible(list(perc_fora = perc_fora, 
                        is_out = pontos_fora, 
                        residuos_fora = r_obs_ord[pontos_fora])))
}

resultado <- qqplot_envelope(r_cond = r_cond, B = 100, conf = 0.95)

# --- Gráfico 3: Histograma com Curva Normal Padrão ---
#hist(r_cond, freq = FALSE, main = "", xlab = "Resíduos quantílicos condicionais", 
#     ylab = "Densidade", col = "gray90", border = "gray40", las = 1)
#x_grid <- seq(min(r_cond, na.rm = TRUE) - 0.5, max(r_cond, na.rm = TRUE) + 0.5, length.out = 200)
#lines(x_grid, dnorm(x_grid), col = "black", lwd = 1.5)

par(mfrow = c(1, 1))

# ==============================================================================
# 3. FUNÇÕES AUXILIARES PARA REAJUSTE E ANÁLISE DE SENSIBILIDADE
# ==============================================================================
calcular_var_perc <- function(tabela_est) {
  ref_completo <- tabela_est[, "Completo"]
  var_perc <- apply(tabela_est, 2, function(col) {
    ((col - ref_completo) / abs(ref_completo)) * 100
  })
  return(round(var_perc, 2))
}

calcular_var_perc_ep <- function(tabela_ep) {
  ref_completo <- tabela_ep[, "Completo"]
  var_perc <- apply(tabela_ep, 2, function(col) {
    ((col - ref_completo) / ref_completo) * 100
  })
  
  return(round(var_perc, 2))
}

reajustar_modelo_1RE <- function(linhas_remover = NULL, y, X, S, id, gh, init_par, k_val) {
  if (!is.null(linhas_remover) && length(linhas_remover) > 0) {
    y_sub <- y[-linhas_remover]
    X_sub <- X[-linhas_remover, , drop = FALSE]
    S_sub <- S[-linhas_remover, , drop = FALSE]
    id_sub <- id[-linhas_remover]
  } else {
    y_sub <- y; X_sub <- X; S_sub <- S; id_sub <- id
  }
  
  fit <- optim(par = init_par, 
               fn = log_vero_PLIA_vet, 
               y = y_sub, X = X_sub, S = S_sub, id = id_sub, 
               nos = gh$nodes, pesos = gh$weights,
               k = k_val, 
               method = "BFGS", 
               hessian = TRUE)
  
  cov_mat <- tryCatch(
    solve(fit$hessian), 
    error = function(e) MASS::ginv(fit$hessian)
  )
  
  se <- sqrt(abs(diag(cov_mat)))
  est <- fit$par
  
  idx_lambda <- length(est) - 1
  idx_sigma2 <- length(est)
  
  val_lambda <- exp(est[idx_lambda])
  val_sigma2 <- exp(est[idx_sigma2])
  
  est[idx_lambda] <- val_lambda
  est[idx_sigma2] <- val_sigma2
  
  se[idx_lambda] <- val_lambda * se[idx_lambda]
  se[idx_sigma2] <- val_sigma2 * se[idx_sigma2]
  
  z_stat <- est / se
  z_stat[idx_lambda] <- (est[idx_lambda] - 1) / se[idx_lambda]
  
  p_val <- 2 * (1 - pnorm(abs(z_stat)))
  
  wald_stat_sigma2 <- (est[idx_sigma2] / se[idx_sigma2])^2
  p_val[idx_sigma2] <- 0.5 * pchisq(wald_stat_sigma2, df = 1, lower.tail = FALSE)
  
  if (is.nan(p_val[idx_lambda]) || is.na(p_val[idx_lambda])) {
    if (est[idx_lambda] < 0.01 || est[idx_lambda] > 100) p_val[idx_lambda] <- 0.0000
  }
  
  resumo <- data.frame(
    Parametro = c(colnames(X_sub), paste0("tau_", colnames(S_sub)), "lambda", "sigma2"),
    Estimativa = est,
    EP = se,
    p_valor = p_val
  )
  return(list(fit = fit, resumo = resumo))
}

reajustar_modelo_2RE <- function(linhas_remover = NULL, y, X, S, id, gh2, init_par, k_val) {
  if (!is.null(linhas_remover) && length(linhas_remover) > 0) {
    y_sub <- y[-linhas_remover]
    X_sub <- X[-linhas_remover, , drop = FALSE]
    S_sub <- S[-linhas_remover, , drop = FALSE]
    id_sub <- id[-linhas_remover]
  } else {
    y_sub <- y; X_sub <- X; S_sub <- S; id_sub <- id
  }
  
  fit <- optim(par = init_par, 
               fn = log_vero_PLIA_2RE, 
               y = y_sub, X = X_sub, S = S_sub, id = id_sub, 
               nos = gh2$nodes, pesos = gh2$weights,
               k = k_val, 
               method = "BFGS", 
               hessian = TRUE)
  
  cov_mat <- tryCatch(
    solve(fit$hessian), 
    error = function(e) MASS::ginv(fit$hessian)
  )
  
  se <- sqrt(abs(diag(cov_mat)))
  est <- fit$par
  
  idx_lambda <- length(est) - 2
  idx_sigma2 <- length(est) - 1
  idx_ups2   <- length(est)
  
  val_lambda <- exp(est[idx_lambda])
  val_sigma2 <- exp(est[idx_sigma2])
  val_ups2   <- exp(est[idx_ups2])
  
  est[idx_lambda] <- val_lambda
  est[idx_sigma2] <- val_sigma2
  est[idx_ups2]   <- val_ups2
  
  se[idx_lambda] <- val_lambda * se[idx_lambda]
  se[idx_sigma2] <- val_sigma2 * se[idx_sigma2]
  se[idx_ups2]   <- val_ups2   * se[idx_ups2]
  
  z_stat <- est / se
  z_stat[idx_lambda] <- (est[idx_lambda] - 1) / se[idx_lambda]
  
  p_val <- 2 * (1 - pnorm(abs(z_stat)))
  
  wald_stat_sigma2 <- (est[idx_sigma2] / se[idx_sigma2])^2
  wald_stat_ups2   <- (est[idx_ups2]   / se[idx_ups2])^2
  p_val[idx_sigma2] <- 0.5 * pchisq(wald_stat_sigma2, df = 1, lower.tail = FALSE)
  p_val[idx_ups2]   <- 0.5 * pchisq(wald_stat_ups2,   df = 1, lower.tail = FALSE)
  
  if (is.nan(p_val[idx_lambda]) || is.na(p_val[idx_lambda])) {
    if (est[idx_lambda] < 0.01 || est[idx_lambda] > 100) p_val[idx_lambda] <- 0.0000
  }
  
  resumo <- data.frame(
    Parametro = c(colnames(X_sub), paste0("tau_", colnames(S_sub)), "lambda", "sigma2", "ups2"),
    Estimativa = est,
    EP = se,
    p_valor = p_val
  )
  return(list(fit = fit, resumo = resumo))
}

# ==============================================================================
# 4. REEXECUÇÃO DOS CENÁRIOS (MODELO 1)
# ==============================================================================
cenarios_m1 <- list(
  "Completo"               = NULL,
  "Sem 34"                 = 34,
  "Sem 92"                 = 92,
  "Sem 97"                 = 97,
  "Sem {34,92,97}"         = c(34, 92, 97)
)

resultados_m1 <- list()
for (nome in names(cenarios_m1)) {
  resultados_m1[[nome]] <- reajustar_modelo_1RE(
    linhas_remover = cenarios_m1[[nome]], 
    y = y, X = X, S = S_1, id = id, gh = gh, 
    init_par = fit_1$par, 
    k_val = k_1
  )
}

tabela_m1_est <- sapply(resultados_m1, function(res) round(res$resumo$Estimativa, 4))
tabela_m1_ep  <- sapply(resultados_m1, function(res) round(res$resumo$EP , 4))
tabela_m1_p   <- sapply(resultados_m1, function(res) round(res$resumo$p_valor, 4))
rownames(tabela_m1_est) <- resultados_m1$Completo$resumo$Parametro
rownames(tabela_m1_ep) <- resultados_m1$Completo$resumo$Parametro
rownames(tabela_m1_p)   <- resultados_m1$Completo$resumo$Parametro

tabela_m1_var_perc <- calcular_var_perc(tabela_m1_est)
tabela_m1_var_perc_ep <- calcular_var_perc_ep(tabela_m1_ep)

cat("\n=========================================================\n")
cat("MODELO 1: ESTIMATIVAS DOS PARÂMETROS\n")
cat("=========================================================\n")
print(tabela_m1_est)

cat("\n=========================================================\n")
cat("MODELO 1: ERRO PADRÃO DOS PARÂMETROS\n")
cat("=========================================================\n")
print(tabela_m1_ep)

cat("\n=========================================================\n")
cat("MODELO 1: VARIAÇÃO PERCENTUAL EM RELAÇÃO AO COMPLETO (%)\n")
cat("=========================================================\n")
print(tabela_m1_var_perc)

cat("\n=========================================================\n")
cat("MODELO 1: VARIAÇÃO DO EP EM RELAÇÃO AO COMPLETO (%)\n")
cat("=========================================================\n")
print(tabela_m1_var_perc_ep)

cat("\n=========================================================\n")
cat("MODELO 1: P-VALORES NOS DIFERENTES CENÁRIOS\n")
cat("=========================================================\n")
print(tabela_m1_p)

# ==============================================================================
# 4. REEXECUÇÃO DOS CENÁRIOS (MODELO 2)
# ==============================================================================
cenarios_m2 <- list(
  "Completo"               = NULL,
  "Sem 21"                 = 21,
  "Sem 54"                 = 54,
  "Sem 83"                 = 83,
  "Sem 97"                 = 97,
  "Sem 127"                = 127,
  "Sem {21,54,83,97,127}"  = c(21,54,83,97,127)
)

resultados_m2 <- list()
for (nome in names(cenarios_m2)) {
  resultados_m2[[nome]] <- reajustar_modelo_1RE(
    linhas_remover = cenarios_m2[[nome]], 
    y = y, X = X, S = S_2, id = id, gh = gh, 
    init_par = fit_2$par, 
    k_val = k_2
  )
}

tabela_m2_est <- sapply(resultados_m2, function(res) round(res$resumo$Estimativa, 4))
tabela_m2_ep  <- sapply(resultados_m2, function(res) round(res$resumo$EP , 4))
tabela_m2_p   <- sapply(resultados_m2, function(res) round(res$resumo$p_valor, 4))
rownames(tabela_m2_est) <- resultados_m2$Completo$resumo$Parametro
rownames(tabela_m2_ep)  <- resultados_m2$Completo$resumo$Parametro
rownames(tabela_m2_p)   <- resultados_m2$Completo$resumo$Parametro

# Calcular variação percentual do Modelo 2
tabela_m2_var_perc <- calcular_var_perc(tabela_m2_est)
tabela_m2_var_perc_ep <- calcular_var_perc_ep(tabela_m2_ep)

cat("\n=========================================================\n")
cat("MODELO 2: ESTIMATIVAS DOS PARÂMETROS\n")
cat("=========================================================\n")
print(tabela_m2_est)

cat("\n=========================================================\n")
cat("MODELO 2: ERRO PADRÃO DOS PARÂMETROS\n")
cat("=========================================================\n")
print(tabela_m2_ep)

cat("\n=========================================================\n")
cat("MODELO 2: VARIAÇÃO PERCENTUAL EM RELAÇÃO AO COMPLETO (%)\n")
cat("=========================================================\n")
print(tabela_m2_var_perc)

cat("\n=========================================================\n")
cat("MODELO 2: VARIAÇÃO DO EP EM RELAÇÃO AO COMPLETO (%)\n")
cat("=========================================================\n")
print(tabela_m2_var_perc_ep)

cat("\n=========================================================\n")
cat("MODELO 2: P-VALORES NOS DIFERENTES CENÁRIOS\n")
cat("=========================================================\n")
print(tabela_m2_p)

# ==============================================================================
# 5. REEXECUÇÃO DOS CENÁRIOS (MODELO 3)
# ==============================================================================
cenarios_m3 <- list(
  "Completo"               = NULL,
  "Sem 34"                 = 34,
  "Sem 92"                 = 92,
  "Sem 97"                 = 97,
  "Sem {34,92,97}"         = c(34, 92, 97)
)

resultados_m3 <- list()
for (nome in names(cenarios_m3)) {
  resultados_m3[[nome]] <- reajustar_modelo_2RE(
    linhas_remover = cenarios_m3[[nome]], 
    y = y, X = X, S = S_3, id = id, gh2 = gh2, 
    init_par = fit_3$par, 
    k_val = k_3
  )
}

tabela_m3_est <- sapply(resultados_m3, function(res) round(res$resumo$Estimativa, 4))
tabela_m3_ep  <- sapply(resultados_m3, function(res) round(res$resumo$EP , 4))
tabela_m3_p   <- sapply(resultados_m3, function(res) round(res$resumo$p_valor, 4))
rownames(tabela_m3_est) <- resultados_m3$Completo$resumo$Parametro
rownames(tabela_m3_ep)  <- resultados_m3$Completo$resumo$Parametro
rownames(tabela_m3_p)   <- resultados_m3$Completo$resumo$Parametro

# Calcular variação percentual do Modelo 3
tabela_m3_var_perc <- calcular_var_perc(tabela_m3_est)
tabela_m3_var_perc_ep <- calcular_var_perc_ep(tabela_m3_ep)

cat("\n=========================================================\n")
cat("MODELO 3: ESTIMATIVAS DOS PARÂMETROS\n")
cat("=========================================================\n")
print(tabela_m3_est)

cat("\n=========================================================\n")
cat("MODELO 3: ERRO PADRÃO DOS PARÂMETROS\n")
cat("=========================================================\n")
print(tabela_m3_ep)

cat("\n=========================================================\n")
cat("MODELO 3: VARIAÇÃO PERCENTUAL EM RELAÇÃO AO COMPLETO (%)\n")
cat("=========================================================\n")
print(tabela_m3_var_perc)

cat("\n=========================================================\n")
cat("MODELO 3: VARIAÇÃO DO EP EM RELAÇÃO AO COMPLETO (%)\n")
cat("=========================================================\n")
print(tabela_m3_var_perc_ep)

cat("\n=========================================================\n")
cat("MODELO 3: P-VALORES NOS DIFERENTES CENÁRIOS\n")
cat("=========================================================\n")
print(tabela_m3_p)

