library(statmod)
library(PLreg)

# ==============================================================================
# PLIA-t
# ==============================================================================

log_vero_PLIA_vet <- function(theta, y, X, S, id, nos, pesos, nu) {
  p <- ncol(X)
  q <- ncol(S)
  
  beta      <- theta[1:p]
  tau       <- theta[(p + 1):(p + q)]
  lambda    <- exp(theta[p + q + 1])  
  varsigma2 <- exp(theta[p + q + 2])  # Única variância de efeito aleatório
  
  Q <- length(nos) 
  
  b_0_vet <- sqrt(2 * varsigma2) * nos
  
  ll_total <- 0
  
  for (i in unique(id)) {
    idx <- (id == i)
    y_i <- y[idx]
    X_i <- X[idx, , drop = FALSE]
    S_i <- S[idx, , drop = FALSE]
    n_i <- length(y_i)
    
    # eta1 recebe o efeito aleatório, eta2 agora é fixo para o indivíduo
    eta1_mat <- outer(as.vector(X_i %*% beta), b_0_vet, "+")
    eta2_vec <- as.vector(S_i %*% tau)
    
    mu_mat    <- pmin(pmax(plogis(eta1_mat), 1e-6), 1 - 1e-6)
    sigma_vec <- pmin(pmax(exp(eta2_vec), 1e-4), 50)
    
    y_mat     <- matrix(rep(y_i, Q), nrow = n_i, ncol = Q)
    sigma_mat <- matrix(rep(sigma_vec, Q), nrow = n_i, ncol = Q)
    
    pdf_vec <- suppressWarnings(PLreg::dPL(as.vector(y_mat), as.vector(mu_mat),
                                           as.vector(sigma_mat), lambda, 
                                           zeta = nu, family = "TF"))
    
    if (any(is.na(pdf_vec)) || any(is.infinite(pdf_vec))) return(1e10)
    
    pdf_mat <- matrix(pdf_vec, nrow = n_i, ncol = Q)
    log_pdf_mat <- log(pmax(pdf_mat, .Machine$double.xmin))
    
    log_densidade_conjunta <- colSums(log_pdf_mat) 
    
    termos_soma <- log_densidade_conjunta + log(pesos)
    max_termo <- max(termos_soma)
    
    log_soma_quadratura <- max_termo + log(sum(exp(termos_soma - max_termo)))
    ll_total <- ll_total + log_soma_quadratura
  }
  
  if (is.na(ll_total) || is.infinite(ll_total)) return(1e10)
  
  return(-ll_total)
}

obter_chutes_iniciais <- function(y, X, S) {
  y_adj <- pmin(pmax(y, 1e-5), 1 - 1e-5)
  z_resp <- qlogis(y_adj)
  
  fit_lm <- lm(z_resp ~ X - 1)
  beta_init <- as.vector(coef(fit_lm))
  tau_init <- rep(0, ncol(S))
  log_lambda_init <- log(1.0) 
  log_varsigma2_init <- log(0.1) 
  
  return(c(beta_init, tau_init, log_lambda_init, log_varsigma2_init))
}

gerar_dados_PLIA <- function(N, n_i, beta, tau, lambda, varsigma2, nu) {
  n_total <- N * n_i
  id <- rep(1:N, each = n_i)
  t_ij <- rep(seq(0.1, 1, length.out = n_i), times = N)
  
  X <- cbind(1, t_ij)
  S <- cbind(1, t_ij)
  
  # Apenas o efeito aleatório na média
  b_0i <- rnorm(N, mean = 0, sd = sqrt(varsigma2))
  b_0ij <- rep(b_0i, each = n_i)
  
  eta1 <- as.vector(X %*% beta) + b_0ij
  eta2 <- as.vector(S %*% tau)
  
  mu    <- plogis(eta1)
  sigma <- exp(eta2)
  
  y <- PLreg::rPL(n_total, mu, sigma, lambda, zeta = nu, family = "TF")
  dados <- data.frame(id = factor(id), t_ij = t_ij, y = y)
  
  return(list(dados = dados, X = X, S = S))
}

# ==============================================================================
# 2. CONFIGURAÇÃO DA SIMULAÇÃO (Quadratura Unidimensional)
# ==============================================================================

rodar_simulacao <- function(M, varsigma2, n_medidas, N_amostra, n_nos, 
                            beta = c(-1.5, 1.5), tau = c(-0.8, 1.2), 
                            lambda = 1.5, nu = 4, seed = 123) {
  set.seed(seed)
  
  theta_verdadeiro <- c(beta, tau, lambda, varsigma2)
  
  # Quadratura Unidimensional
  gh <- statmod::gauss.quad(n_nos, kind = "hermite")
  nos <- gh$nodes
  pesos <- gh$weights / sqrt(pi) 
  
  resultados <- matrix(NA, nrow = M, ncol = 6)
  colnames(resultados) <- c("beta0", "beta1", "tau0", "tau1", "lambda", "varsigma2")
  
  sucessos <- 0
  tentativas <- 0
  
  cat("==============================================================\n")
  cat(sprintf("  INICIANDO SIMULAÇÃO (N = %d, var = %.2f)\n", N_amostra, varsigma2))
  cat("==============================================================\n\n")
  flush.console()
  
  while (sucessos < M) {
    tentativas <- tentativas + 1
    
    sim <- gerar_dados_PLIA(N = N_amostra, n_i = n_medidas, 
                            beta = beta, tau = tau, 
                            lambda = lambda, varsigma2 = varsigma2, nu = nu)
    
    theta_init <- obter_chutes_iniciais(sim$dados$y, sim$X, sim$S)
    
    fit_nosso <- tryCatch({
      nlminb(
        start     = theta_init, 
        objective = log_vero_PLIA_vet,
        y         = sim$dados$y, 
        X         = sim$X, 
        S         = sim$S, 
        id        = as.numeric(sim$dados$id), 
        nos       = nos,      
        pesos     = pesos,  
        nu        = nu,  
        # Limites ajustados para menos um parâmetro
        lower     = c(rep(-10, ncol(sim$X)), rep(-10, ncol(sim$S)), -3, -5),
        upper     = c(rep(10, ncol(sim$X)), rep(10, ncol(sim$S)),  3,  5),
        control   = list(eval.max = 2000, iter.max = 2000, rel.tol = 1e-5)
      )
    }, error = function(e) {
      cat(sprintf("\n[ERRO NA VEROSSIMILHANÇA]: %s\n", e$message))
      return(NULL)
    })
    
    if (!is.null(fit_nosso) && fit_nosso$convergence == 0) {
      
      vetor_nosso <- c(
        fit_nosso$par[1:4],         
        exp(fit_nosso$par[5]),      # lambda 
        exp(fit_nosso$par[6])       # varsigma2
      )
      
      sucessos <- sucessos + 1
      resultados[sucessos, ] <- vetor_nosso
      
      cat(sprintf("Réplica %d de %d concluída! (Tentativas totais: %d)\n", 
                  sucessos, M, tentativas))
      flush.console()
    }
  }
  
  taxa_sucesso <- (sucessos / tentativas) * 100
  limpo_nosso <- na.omit(resultados)
  
  cat("\n=========================================\n")
  cat(sprintf("  Taxa de Convergência (PL-t): %.2f%%\n", taxa_sucesso))
  cat("=========================================\n\n")
  flush.console()
  
  medias <- colMeans(limpo_nosso)
  vies_nosso <- medias - theta_verdadeiro
  rEQM_nosso <- sqrt(colMeans(sweep(limpo_nosso, 2, theta_verdadeiro, "-")^2))
  
  tabela_resumo <- data.frame(
    Verdadeiro = theta_verdadeiro,
    Media_Est  = round(medias, 4),
    Vies       = round(vies_nosso, 4),
    rEQM       = round(rEQM_nosso, 4),
    row.names  = c("beta0", "beta1", "tau0", "tau1", "lambda", "varsigma2")
  )
  
  return(list(
    limpo_nosso       = limpo_nosso,
    tabela_resumo     = tabela_resumo,
    taxa_convergencia = taxa_sucesso,
    tentativas_totais = tentativas
  ))
}

# ==============================================================================
# CENÁRIOS
# ==============================================================================
sim21 <- rodar_simulacao(M = 2000, varsigma2 = 0.05, n_medidas = 5, N_amostra = 10,  n_nos = 20)
sim22 <- rodar_simulacao(M = 2000, varsigma2 = 0.05, n_medidas = 5, N_amostra = 20,  n_nos = 20)
sim23 <- rodar_simulacao(M = 2000, varsigma2 = 0.05, n_medidas = 5, N_amostra = 40,  n_nos = 20)
sim24 <- rodar_simulacao(M = 2000, varsigma2 = 0.05, n_medidas = 5, N_amostra = 100, n_nos = 20)
sim25 <- rodar_simulacao(M = 2000, varsigma2 = 0.75, n_medidas = 5, N_amostra = 10,  n_nos = 20)
sim26 <- rodar_simulacao(M = 2000, varsigma2 = 0.75, n_medidas = 5, N_amostra = 20,  n_nos = 20)
sim27 <- rodar_simulacao(M = 2000, varsigma2 = 0.75, n_medidas = 5, N_amostra = 40,  n_nos = 20)
sim28 <- rodar_simulacao(M = 2000, varsigma2 = 0.75, n_medidas = 5, N_amostra = 100, n_nos = 20)

print(sim21$tabela_resumo)
print(sim22$tabela_resumo)
print(sim23$tabela_resumo)
print(sim24$tabela_resumo)
print(sim25$tabela_resumo)
print(sim26$tabela_resumo)
print(sim27$tabela_resumo)
print(sim28$tabela_resumo)

print(sim21$taxa_convergencia)
print(sim22$taxa_convergencia)
print(sim23$taxa_convergencia)
print(sim24$taxa_convergencia)
print(sim25$taxa_convergencia)
print(sim26$taxa_convergencia)
print(sim27$taxa_convergencia)
print(sim28$taxa_convergencia)

saveRDS(sim21$limpo_nosso, file = "plt7_1.rds"); plt7_1 <- readRDS("plt7_1.rds")
saveRDS(sim22$limpo_nosso, file = "plt7_2.rds"); plt7_2 <- readRDS("plt7_2.rds")
saveRDS(sim23$limpo_nosso, file = "plt7_3.rds"); plt7_3 <- readRDS("plt7_3.rds")
saveRDS(sim24$limpo_nosso, file = "plt7_4.rds"); plt7_4 <- readRDS("plt7_4.rds")
saveRDS(sim25$limpo_nosso, file = "plt7_5.rds"); plt7_5 <- readRDS("plt7_5.rds")
saveRDS(sim26$limpo_nosso, file = "plt7_6.rds"); plt7_6 <- readRDS("plt7_6.rds")
saveRDS(sim27$limpo_nosso, file = "plt7_7.rds"); plt7_7 <- readRDS("plt7_7.rds")
saveRDS(sim28$limpo_nosso, file = "plt7_8.rds"); plt7_8 <- readRDS("plt7_8.rds")
