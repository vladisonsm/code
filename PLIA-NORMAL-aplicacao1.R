library(statmod)
library(PLreg)
library(numDeriv)

dados = Planilha.sem.título...OphthalmologicalStudyData
dados$y = as.numeric(gsub(",", ".", dados$y)); dados$y

dados$concentracao <- factor(dados$concentracao)
dados$concentracao <- relevel(dados$concentracao, ref = "-1")  # <- ajuste o nível de referência

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
# 1. FUNÇÕES DE LOG-VEROSSIMILHANÇA
# ==============================================================================
log_vero_PLIA_vet <- function(theta, y, X, S, id, nos, pesos) {
  p <- ncol(X); q <- ncol(S)
  beta   <- theta[1:p]
  tau    <- theta[(p + 1):(p + q)]
  lambda <- exp(theta[p + q + 1])
  sigma2 <- exp(theta[p + q + 2])
  
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
    
    pdf_vec <- suppressWarnings(PLreg::dPL(as.vector(y_mat), as.vector(mu_mat),
                                           as.vector(sigma_mat), lambda, family = "NO"))
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

log_vero_PLIA_2RE <- function(theta, y, X, S, id, nos, pesos) {
  p <- ncol(X); q <- ncol(S)
  beta   <- theta[1:p]
  tau    <- theta[(p + 1):(p + q)]
  lambda <- exp(theta[p + q + 1])
  sigma2 <- exp(theta[p + q + 2])  # variância de b0i (média)
  ups2   <- exp(theta[p + q + 3])  # variância de d0i (dispersão)
  
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
    
    eta1_mat <- outer(xb_i, grade$b0, "+")   # n_i x K
    eta2_mat <- outer(s_i,  grade$d0, "+")   # n_i x K
    
    mu_mat    <- pmin(pmax(plogis(eta1_mat), 1e-6), 1 - 1e-6)
    sigma_mat <- pmin(pmax(exp(eta2_mat), 1e-4), 50)
    y_mat     <- matrix(rep(y_i, K), nrow = n_i, ncol = K)
    
    pdf_vec <- suppressWarnings(PLreg::dPL(as.vector(y_mat), as.vector(mu_mat),
                                           as.vector(sigma_mat), lambda, family = "NO"))
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
obter_chutes_iniciais <- function(y, X, S) {
  y_adj <- pmin(pmax(y, 1e-5), 1 - 1e-5)
  z_resp <- qlogis(y_adj)
  fit_lm <- lm(z_resp ~ X - 1)
  beta_init <- as.vector(coef(fit_lm))
  tau_init  <- c(log(sd(residuals(fit_lm))^2), rep(0, ncol(S) - 1))  # intercepto informativo + zeros
  log_lambda_init <- log(1.0)
  log_sigma2_init <- log(0.1)
  c(beta_init, tau_init, log_lambda_init, log_sigma2_init)
}

obter_chutes_iniciais_2RE <- function(y, X, S) {
  base <- obter_chutes_iniciais(y, X, S)
  c(base, log(0.1))  # log_upsilon2 inicial
}

# --- Modelo nulo (sem NENHUM efeito aleatório) — usado só para o TRV corrigido ---
# É o próprio modelo cheio no limite sigma2 -> 0: mesma família PL, mesmo X e S,
# mas sem quadratura (sem b0i), então a log-verossimilhança é só a soma direta.
log_vero_PL_fixo <- function(theta, y, X, S) {
  p <- ncol(X); q <- ncol(S)
  beta   <- theta[1:p]
  tau    <- theta[(p + 1):(p + q)]
  lambda <- exp(theta[p + q + 1])
  
  eta1 <- as.vector(X %*% beta)
  eta2 <- as.vector(S %*% tau)
  mu    <- pmin(pmax(plogis(eta1), 1e-6), 1 - 1e-6)
  sigma <- pmin(pmax(exp(eta2), 1e-4), 50)
  
  pdf_vec <- suppressWarnings(PLreg::dPL(y, mu, sigma, lambda, family = "NO"))
  if (any(is.na(pdf_vec)) || any(is.infinite(pdf_vec))) return(1e10)
  ll <- sum(log(pmax(pdf_vec, .Machine$double.xmin)))
  if (is.na(ll) || is.infinite(ll)) return(1e10)
  return(-ll)
}

obter_chutes_iniciais_fixo <- function(y, X, S) {
  base <- obter_chutes_iniciais(y, X, S)
  base[-length(base)]  # descarta o log_sigma2 (não existe no modelo nulo)
}

# Estatística e p-valor do TRV com a correção de fronteira de Self & Liang (1987):
# sob H0 (variância = 0, no limite do espaço paramétrico), a estatística do TRV
# segue assintoticamente uma mistura 0,5*chi2_0 + 0,5*chi2_1 em vez de chi2_1 pura.
teste_lrt_fronteira <- function(logLik_completo, logLik_nulo, df_extra = 1) {
  LR <- 2 * (logLik_completo - logLik_nulo)
  LR <- max(LR, 0)  # por construção o modelo cheio não pode ter logLik menor
  p_valor <- 0.5 * (1 - pchisq(LR, df = df_extra))
  data.frame(LR = LR, p_valor = p_valor)
}

# ==============================================================================
# 3. FUNÇÃO PARA MONTAR A TABELA (estimativa, EP, z, valor-p)
# ==============================================================================
construir_tabela <- function(fit, hess, nomes, idx_lambda = NULL, idx_var = NULL) {
  vcov_par <- tryCatch(solve(hess), error = function(e) {
    warning("Hessiana singular no ponto ótimo: tente rel.tol menor ou outro chute inicial")
    matrix(NA, nrow(hess), ncol(hess))
  })
  se_par <- sqrt(diag(vcov_par))
  estim  <- fit$par
  
  Estimativa <- estim
  EP         <- se_par
  Z          <- estim / se_par
  p_valor    <- 2 * (1 - pnorm(abs(Z)))
  
  # metodo delta: se g(theta) = exp(theta), Var(g) ~= exp(theta)^2 * Var(theta)
  for (j in c(idx_lambda, idx_var)) {
    valor_original <- exp(estim[j])
    Estimativa[j]  <- valor_original
    EP[j]          <- valor_original * se_par[j]
    nomes[j]       <- sub("^log_", "", nomes[j])
  }
  # para lambda, o teste substantivo é H0: lambda = 1 (não H0: lambda = 0)
  for (j in idx_lambda) {
    Z[j]       <- (Estimativa[j] - 1) / EP[j]
    p_valor[j] <- 2 * (1 - pnorm(abs(Z[j])))
  }
  # para variâncias, H0: variância = 0 é problema de fronteira -> Wald usual não vale
  for (j in idx_var) { Z[j] <- NA; p_valor[j] <- NA }
  
  data.frame(Parametro = nomes, Estimativa = Estimativa, EP = EP, Z = Z, p_valor = p_valor)
}

# ==============================================================================
# 4. AJUSTE DO MODELO 1  (log(delta) = tau0 ; sem efeito aleatório na dispersão)
# ==============================================================================
S_1 <- matrix(1, nrow = length(y), ncol = 1); colnames(S_1) <- "Intercepto"
q1 <- ncol(S_1)
gh <- statmod::gauss.quad(20, kind = "hermite")

theta_init_1 <- obter_chutes_iniciais(y, X, S_1)

fit_1 <- nlminb(
  start = theta_init_1, objective = log_vero_PLIA_vet,
  y = y, X = X, S = S_1, id = id, nos = gh$nodes, pesos = gh$weights,
  lower = c(rep(-10, p), rep(-10, q1), -3, -5),
  upper = c(rep(10, p),  rep(10, q1),  3,  5),
  control = list(eval.max = 3000, iter.max = 3000, rel.tol = 1e-8)
)
hess_1 <- numDeriv::hessian(log_vero_PLIA_vet, x = fit_1$par,
                             y = y, X = X, S = S_1, id = id, nos = gh$nodes, pesos = gh$weights)
nomes_1 <- c(colnames(X), "tau0", "log_lambda", "log_sigma2_b")
tabela_1 <- construir_tabela(fit_1, hess_1, nomes_1,
                              idx_lambda = p + q1 + 1, idx_var = p + q1 + 2)

# ==============================================================================
# 5. AJUSTE DO MODELO 2  (log(delta) = tau0 + tau1*log(t))
# ==============================================================================
S_2 <- cbind(1, dados$log_t); colnames(S_2) <- c("Intercepto", "log_t")
q2 <- ncol(S_2)

theta_init_2 <- obter_chutes_iniciais(y, X, S_2)

fit_2 <- nlminb(
  start = theta_init_2, objective = log_vero_PLIA_vet,
  y = y, X = X, S = S_2, id = id, nos = gh$nodes, pesos = gh$weights,
  lower = c(rep(-10, p), rep(-10, q2), -3, -5),
  upper = c(rep(10, p),  rep(10, q2),  3,  5),
  control = list(eval.max = 3000, iter.max = 3000, rel.tol = 1e-8)
)
hess_2 <- numDeriv::hessian(log_vero_PLIA_vet, x = fit_2$par,
                             y = y, X = X, S = S_2, id = id, nos = gh$nodes, pesos = gh$weights)
nomes_2 <- c(colnames(X), "tau0", "tau1_logt", "log_lambda", "log_sigma2_b")
tabela_2 <- construir_tabela(fit_2, hess_2, nomes_2,
                              idx_lambda = p + q2 + 1, idx_var = p + q2 + 2)

# ==============================================================================
# 6. AJUSTE DO MODELO 3  (log(delta) = d0i + tau0 ; efeito aleatório na dispersão)
# ==============================================================================
S_3 <- matrix(1, nrow = length(y), ncol = 1); colnames(S_3) <- "Intercepto"
q3 <- ncol(S_3)

# quadratura menor pois a integral dupla usa Q^2 pontos por indivíduo
gh2 <- statmod::gauss.quad(20, kind = "hermite")
theta_init_3 <- obter_chutes_iniciais_2RE(y, X, S_3)

fit_3 <- nlminb(
  start = theta_init_3, objective = log_vero_PLIA_2RE,
  y = y, X = X, S = S_3, id = id, nos = gh2$nodes, pesos = gh2$weights,
  lower = c(rep(-10, p), rep(-10, q3), -3, -5, -5),
  upper = c(rep(10, p),  rep(10, q3),  3,  5,  5),
  control = list(eval.max = 3000, iter.max = 3000, rel.tol = 1e-8)
)
hess_3 <- numDeriv::hessian(log_vero_PLIA_2RE, x = fit_3$par,
                             y = y, X = X, S = S_3, id = id, nos = gh2$nodes, pesos = gh2$weights)
nomes_3 <- c(colnames(X), "tau0", "log_lambda", "log_sigma2_b", "log_upsilon2_d")
tabela_3 <- construir_tabela(fit_3, hess_3, nomes_3,
                              idx_lambda = p + q3 + 1, idx_var = c(p + q3 + 2, p + q3 + 3))

# ==============================================================================
#  MODELO NULO PARA TESTAR sigma2_b NO MODELO 3 (apenas upsilon2_d presente)
# ==============================================================================
log_vero_PL_1RE_disp <- function(theta, y, X, S, id, nos, pesos) {
  p <- ncol(X); q <- ncol(S)
  beta   <- theta[1:p]
  tau    <- theta[(p + 1):(p + q)]
  lambda <- exp(theta[p + q + 1])
  ups2   <- exp(theta[p + q + 2])  # variância de d0i (dispersão) apenas
  
  d0_vet <- sqrt(2 * ups2) * nos
  pesos_prod <- pesos / sqrt(pi)
  K <- length(nos)
  
  ll_total <- 0
  for (i in unique(id)) {
    idx <- (id == i)
    y_i <- y[idx]; X_i <- X[idx, , drop = FALSE]; S_i <- S[idx, , drop = FALSE]
    n_i <- length(y_i)
    
    # eta1 é fixo (sem b0i), eta2 tem d0i
    eta1_mat <- matrix(rep(as.vector(X_i %*% beta), K), nrow = n_i, ncol = K)
    eta2_mat <- outer(as.vector(S_i %*% tau), d0_vet, "+")
    
    mu_mat    <- pmin(pmax(plogis(eta1_mat), 1e-6), 1 - 1e-6)
    sigma_mat <- pmin(pmax(exp(eta2_mat), 1e-4), 50)
    y_mat     <- matrix(rep(y_i, K), nrow = n_i, ncol = K)
    
    pdf_vec <- suppressWarnings(PLreg::dPL(as.vector(y_mat), as.vector(mu_mat),
                                           as.vector(sigma_mat), lambda, family = "NO"))
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

# Chute inicial (beta, tau, lambda, upsilon2_d)
obter_chutes_iniciais_1RE_disp <- function(y, X, S) {
  base <- obter_chutes_iniciais(y, X, S)
  c(base[1:(length(base)-1)], log(0.1)) # descarta log_sigma2_b e coloca log_upsilon2_d
}

theta_init_3_nulo <- obter_chutes_iniciais_1RE_disp(y, X, S_3)
fit_3_nulo_sigma2b <- nlminb(
  start = theta_init_3_nulo, objective = log_vero_PL_1RE_disp,
  y = y, X = X, S = S_3, id = id, nos = gh$nodes, pesos = gh$weights, # usando gh (1D)
  lower = c(rep(-10, p), rep(-10, q3), -3, -5),
  upper = c(rep(10, p),  rep(10, q3),  3,  5),
  control = list(eval.max = 3000, iter.max = 3000, rel.tol = 1e-8)
)

# ==============================================================================
# 7. TESTE DA RAZÃO DE VEROSSIMILHANÇAS CORRIGIDO (efeitos aleatórios)
# ==============================================================================
# H0: sigma2_b = 0 no Modelo 1 -> compara com o modelo nulo (mesmo X e S_1,
# sem NENHUM efeito aleatório na média)
theta_init_1_nulo <- obter_chutes_iniciais_fixo(y, X, S_1)
fit_1_nulo <- nlminb(
  start = theta_init_1_nulo, objective = log_vero_PL_fixo,
  y = y, X = X, S = S_1,
  lower = c(rep(-10, p), rep(-10, q1), -3),
  upper = c(rep(10, p),  rep(10, q1),  3),
  control = list(eval.max = 3000, iter.max = 3000, rel.tol = 1e-8)
)
lrt_sigma2b_1 <- teste_lrt_fronteira(-fit_1$objective, -fit_1_nulo$objective)

# H0: sigma2_b = 0 no Modelo 2 -> compara com o modelo nulo (mesmo X e S_2,
# sem NENHUM efeito aleatório na média)
theta_init_2_nulo <- obter_chutes_iniciais_fixo(y, X, S_2)
fit_2_nulo <- nlminb(
  start = theta_init_2_nulo, objective = log_vero_PL_fixo,
  y = y, X = X, S = S_2,
  lower = c(rep(-10, p), rep(-10, q2), -3),
  upper = c(rep(10, p),  rep(10, q2),  3),
  control = list(eval.max = 3000, iter.max = 3000, rel.tol = 1e-8)
)
lrt_sigma2b_2 <- teste_lrt_fronteira(-fit_2$objective, -fit_2_nulo$objective)

# H0: upsilon2_d = 0 no Modelo 3 -> ao remover d0i, o Modelo 3 se reduz
# EXATAMENTE ao Modelo 1 (mesmo X com b0i, mesma dispersão só com intercepto),
# então o modelo nulo aqui já está ajustado: é o próprio fit_1.
lrt_upsilon2d_3 <- teste_lrt_fronteira(-fit_3$objective, -fit_1$objective)
lrt_sigma2b_3 <- teste_lrt_fronteira(-fit_3$objective, -fit_3_nulo_sigma2b$objective)

tabela_lrt <- data.frame(
  Teste = c("H0: sigma2_b = 0  (Modelo 1)",
            "H0: sigma2_b = 0  (Modelo 2)",
            "H0: upsilon2_d = 0 (Modelo 3, dado b0i já presente)",
            "H0: sigma2_b = 0 (Modelo 3, dado d0i já presente)"),
  LR      = c(lrt_sigma2b_1$LR, lrt_sigma2b_2$LR, lrt_upsilon2d_3$LR, lrt_sigma2b_3$LR),
  p_valor = c(lrt_sigma2b_1$p_valor, lrt_sigma2b_2$p_valor, lrt_upsilon2d_3$p_valor, lrt_sigma2b_3$p_valor)
)
tabela_1$p_valor[tabela_1$Parametro == "sigma2_b"] <- lrt_sigma2b_1$p_valor
tabela_2$p_valor[tabela_2$Parametro == "sigma2_b"] <- lrt_sigma2b_2$p_valor
  tabela_3$p_valor[tabela_3$Parametro == "sigma2_b"] <- lrt_sigma2b_3$p_valor
tabela_3$p_valor[tabela_3$Parametro == "upsilon2_d"] <- lrt_upsilon2d_3$p_valor

# ==============================================================================
# 8. RESULTADOS
# ==============================================================================
cat("\n--- Modelo 1 ---\n"); print(tabela_1, row.names = FALSE)
cat("\n--- Modelo 2 ---\n"); print(tabela_2, row.names = FALSE)
cat("\n--- Modelo 3 ---\n"); print(tabela_3, row.names = FALSE)
cat("\n--- TRV corrigido para os efeitos aleatórios (mistura 0,5*chi2_0 + 0,5*chi2_1) ---\n")
print(tabela_lrt, row.names = FALSE)

comparar_modelos <- function(fits, nomes_modelos, n) {
  data.frame(
    Modelo       = nomes_modelos,
    npar         = sapply(fits, function(f) length(f$par)),
    logLik       = sapply(fits, function(f) -f$objective),
    AIC          = sapply(fits, function(f) 2 * (length(f$par)) + 2 * f$objective),
    BIC          = sapply(fits, function(f) length(f$par) * log(n) + 2 * f$objective),
    CAIC         = sapply(fits, function(f) length(f$par) * log(n + 1) + 2 * f$objective),
    HQIC         = sapply(fits, function(f) length(f$par) * log(log(n)) + 2 * f$objective),
    Convergencia = sapply(fits, function(f) f$convergence)
  )
}
tabela_comparativa <- comparar_modelos(list(fit_1, fit_2, fit_3),
                                       c("Modelo 1", "Modelo 2", "Modelo 3"),
                                       n = length(y))
cat("\n--- Comparação dos modelos ---\n"); print(tabela_comparativa, row.names = FALSE)

# ==============================================================================
# 1. FUNÇÕES PARA CÁLCULO DOS RESÍDUOS QUANTÍLICOS CONDICIONAIS
# ==============================================================================
# Função para Modelos com 1 Efeito Aleatório (Modelo 1 e Modelo 2)
calcular_residuos_1RE <- function(fit, y, X, S, id, gh) {
  p <- ncol(X); q <- ncol(S)
  beta   <- fit$par[1:p]
  tau    <- fit$par[(p + 1):(p + q)]
  lambda <- exp(fit$par[p + q + 1])
  sigma2 <- exp(fit$par[p + q + 2])
  
  nos <- gh$nodes; pesos <- gh$weights
  Q <- length(nos)
  b_0_vet <- sqrt(2 * sigma2) * nos
  pesos_prod <- pesos / sqrt(pi)
  
  ids_unicos <- unique(id)
  b0_ebe <- numeric(length(ids_unicos))
  names(b0_ebe) <- ids_unicos
  
  # Estimativa Bayesiana Empírica (EBE) de b0i para cada indivíduo
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
                                           as.vector(sigma_mat), lambda, family = "NO"))
    pdf_mat <- matrix(pdf_vec, nrow = n_i, ncol = Q)
    log_pdf_mat <- log(pmax(pdf_mat, .Machine$double.xmin))
    log_dens_conjunta <- colSums(log_pdf_mat)
    
    termos <- log_dens_conjunta + log(pesos_prod)
    w_post <- exp(termos - max(termos)) / sum(exp(termos - max(termos)))
    
    b0_ebe[as.character(i)] <- sum(b_0_vet * w_post)
  }
  
  # Predição condicional e resíduo quantílico
  b0_i_vet <- b0_ebe[as.character(id)]
  eta1_cond <- as.vector(X %*% beta) + b0_i_vet
  eta2_cond <- as.vector(S %*% tau)
  
  mu_cond    <- pmin(pmax(plogis(eta1_cond), 1e-6), 1 - 1e-6)
  sigma_cond <- pmin(pmax(exp(eta2_cond), 1e-4), 50)
  
  u_cond <- suppressWarnings(PLreg::pPL(y, mu_cond, sigma_cond, lambda, family = "NO"))
  u_cond <- pmin(pmax(u_cond, 1e-7), 1 - 1e-7)
  
  r_cond <- qnorm(u_cond)
  return(list(r_cond = r_cond, b0_ebe = b0_ebe))
}

# Função para Modelo com 2 Efeitos Aleatórios (Modelo 3)
calcular_residuos_2RE <- function(fit, y, X, S, id, gh) {
  p <- ncol(X); q <- ncol(S)
  beta   <- fit$par[1:p]
  tau    <- fit$par[(p + 1):(p + q)]
  lambda <- exp(fit$par[p + q + 1])
  sigma2 <- exp(fit$par[p + q + 2])
  ups2   <- exp(fit$par[p + q + 3])
  
  nos <- gh$nodes; pesos <- gh$weights
  b0_vet <- sqrt(2 * sigma2) * nos
  d0_vet <- sqrt(2 * ups2) * nos
  
  grade  <- expand.grid(b0 = b0_vet, d0 = d0_vet)
  pesos2 <- expand.grid(wb = pesos, wd = pesos)
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
                                           as.vector(sigma_mat), lambda, family = "NO"))
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
  
  u_cond <- suppressWarnings(PLreg::pPL(y, mu_cond, sigma_cond, lambda, family = "NO"))
  u_cond <- pmin(pmax(u_cond, 1e-7), 1 - 1e-7)
  
  r_cond <- qnorm(u_cond)
  return(list(r_cond = r_cond, b0_ebe = b0_ebe, d0_ebe = d0_ebe))
}

# ==============================================================================
# 2. CÁLCULO PARA O MODELO SELECIONADO E GERAÇÃO DO GRÁFICO
# ==============================================================================
par(mfrow=c(1,3))
# Para analisar o Modelo 1:
# res_obj <- calcular_residuos_1RE(fit = fit_1, y = y, X = X, S = S_1, id = id, gh = gh)

# Para analisar o Modelo 2:
# res_obj <- calcular_residuos_1RE(fit = fit_2, y = y, X = X, S = S_2, id = id, gh = gh)

# Para analisar o Modelo 3:
# 
res_obj <- calcular_residuos_2RE(fit = fit_3, y = y, X = X, S = S_3, id = id, gh = gh2)

r_cond <- res_obj$r_cond
obs_id <- 1:length(y) # Número da linha da observação no dataset

# Configurar o gráfico
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

# ==============================================================================
# FUNÇÃO PARA GERAR O QQ-PLOT COM ENVELOPES SIMULADOS E % DE PONTOS FORA
# ==============================================================================
par(mfrow=c(1,3))
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

# Para ver valores que ficaram de fora:
# print(resultado$residuos_fora)

# ==============================================================================
# FUNÇÃO AUXILIAR PARA REAJUSTAR MODELO E EXTRAIR INFERÊNCIA 
# (COM MISTURA 50:50 DE QUI-QUADRADOS PARA A VARIÂNCIA DO EFEITO ALEATÓRIO)
# ==============================================================================
reajustar_modelo_1RE <- function(linhas_remover = NULL, y, X, S, id, gh, init_par) {
  if (!is.null(linhas_remover) && length(linhas_remover) > 0) {
    y_sub <- y[-linhas_remover]
    X_sub <- X[-linhas_remover, , drop = FALSE]
    S_sub <- S[-linhas_remover, , drop = FALSE]
    id_sub <- id[-linhas_remover]
  } else {
    y_sub <- y; X_sub <- X; S_sub <- S; id_sub <- id
  }
  
  # Otimização
  fit <- optim(par = init_par, 
               fn = log_vero_PLIA_vet, 
               y = y_sub, X = X_sub, S = S_sub, id = id_sub, 
               nos = gh$nodes, pesos = gh$weights,
               method = "BFGS", 
               hessian = TRUE)
  
  # Extração da Matriz de Covariância usando Inversa Generalizada 
  # (mais robusta contra instabilidade numérica que o 'solve' padrão)
  cov_mat <- tryCatch(
    solve(fit$hessian), 
    error = function(e) MASS::ginv(fit$hessian)
  )
  
  # Força valores absolutos na diagonal para evitar variâncias negativas por arredondamento de máquina
  se <- sqrt(abs(diag(cov_mat)))
  est <- fit$par
  
  # Índices
  idx_lambda <- length(est) - 1
  idx_sigma2 <- length(est)
  
  # Transformação Exponencial
  val_lambda <- exp(est[idx_lambda])
  val_sigma2 <- exp(est[idx_sigma2])
  
  est[idx_lambda] <- val_lambda
  est[idx_sigma2] <- val_sigma2
  
  # Método Delta
  se[idx_lambda] <- val_lambda * se[idx_lambda]
  se[idx_sigma2] <- val_sigma2 * se[idx_sigma2]
  
  # Estatísticas Z
  z_stat <- est / se
  z_stat[idx_lambda] <- (est[idx_lambda] - 1) / se[idx_lambda] # H0: lambda = 1
  
  # P-valores normais
  p_val <- 2 * (1 - pnorm(abs(z_stat)))
  
  # P-valor para sigma2 (Mistura 50:50)
  wald_stat_sigma2 <- (est[idx_sigma2] / se[idx_sigma2])^2
  p_val[idx_sigma2] <- 0.5 * pchisq(wald_stat_sigma2, df = 1, lower.tail = FALSE)
  
  # TRATAMENTO DO NaN PARA O LAMBDA
  if (is.nan(p_val[idx_lambda]) || is.na(p_val[idx_lambda])) {
    # Se a variância explodiu, mas a estimativa está muito abaixo de 1 (tendendo a zero),
    # rejeitamos fortemente a hipótese nula de que lambda = 1.
    if (est[idx_lambda] < 0.01 || est[idx_lambda] > 100) {
      p_val[idx_lambda] <- 0.0000
    }
  }
  
  resumo <- data.frame(
    Parametro = c(colnames(X_sub), paste0("tau_", colnames(S_sub)), "lambda", "sigma2"),
    Estimativa = est,
    EP = se,
    p_valor = p_val
  )
  
  return(list(fit = fit, resumo = resumo))
}

# Função para calcular a variação percentual em relação ao modelo Completo
calcular_var_perc <- function(tabela_est) {
  ref_completo <- tabela_est[, "Completo"]
  var_perc <- apply(tabela_est, 2, function(col) {
    ((col - ref_completo) / abs(ref_completo)) * 100
  })
  return(round(var_perc, 2))
}

# Função para calcular a variação percentual do Erro Padrão (EP)
calcular_var_perc_ep <- function(tabela_ep) {
  ref_completo <- tabela_ep[, "Completo"]
  var_perc <- apply(tabela_ep, 2, function(col) {
    ((col - ref_completo) / ref_completo) * 100
  })
  
  return(round(var_perc, 2))
}

# ==============================================================================
# MODELO 1: EXECUÇÃO E TABELAS COMPARATIVAS
# ==============================================================================
cenarios_m1 <- list(
  "Completo"      = NULL,
  "Sem 21"        = 21,
  "Sem 94"        = 94,
  "Sem 97"        = 97,
  "Sem {21,94,97}"= c(21, 94, 97)
)

resultados_m1 <- list()
for (nome in names(cenarios_m1)) {
  resultados_m1[[nome]] <- reajustar_modelo_1RE(
    linhas_remover = cenarios_m1[[nome]], 
    y = y, X = X, S = S_1, id = id, gh = gh, init_par = fit_1$par
  )
}

tabela_m1_est <- sapply(resultados_m1, function(res) round(res$resumo$Estimativa, 4))
tabela_m1_ep  <- sapply(resultados_m1, function(res) round(res$resumo$EP , 4))
tabela_m1_p   <- sapply(resultados_m1, function(res) round(res$resumo$p_valor, 4))
rownames(tabela_m1_est) <- resultados_m1$Completo$resumo$Parametro
rownames(tabela_m1_ep) <- resultados_m1$Completo$resumo$Parametro
rownames(tabela_m1_p)   <- resultados_m1$Completo$resumo$Parametro

# Calcular variação percentual do Modelo 1
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
# MODELO 2: EXECUÇÃO E TABELAS COMPARATIVAS
# ==============================================================================
cenarios_m2 <- list(
  "Completo"              = NULL,
  "Sem 21"                = 21,
  "Sem 54"                = 54,
  "Sem 94"                = 94,
  "Sem 97"                = 97,
  "Sem 127"               = 127,
  "Sem {21,54,94,97,127}" = c(21, 54, 94, 97, 127)
)

resultados_m2 <- list()
for (nome in names(cenarios_m2)) {
  resultados_m2[[nome]] <- reajustar_modelo_1RE(
    linhas_remover = cenarios_m2[[nome]], 
    y = y, X = X, S = S_2, id = id, gh = gh, init_par = fit_2$par
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
# MODELO 3: EXECUÇÃO E TABELAS COMPARATIVAS
# ==============================================================================
cenarios_m3 <- list(
  "Completo"              = NULL,
  "Sem 21"                = 21,
  "Sem 54"                = 54,
  "Sem 94"                = 94,
  "Sem 97"                = 97,
  "Sem 127"               = 127,
  "Sem {21,54,94,97,127}" = c(21, 54, 94, 97, 127)
)

resultados_m3 <- list()
for (nome in names(cenarios_m3)) {
  resultados_m3[[nome]] <- reajustar_modelo_2RE(
    linhas_remover = cenarios_m3[[nome]], 
    y = y, X = X, S = S_3, id = id, gh = gh, init_par = fit_3$par
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
