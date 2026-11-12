library(statmod)
library(PLreg)

# ==============================================================================
# PLIA-no
# ==============================================================================

log_vero_PLIA_vet <- function(theta, y, X, S, id, nos, pesos) {
  p <- ncol(X)
  q <- ncol(S)
  
  beta       <- theta[1:p]
  tau        <- theta[(p + 1):(p + q)]
  lambda     <- exp(theta[p + q + 1])  # Traz log(lambda) de volta para lambda
  sigma2     <- exp(theta[p + q + 2])  # Traz log(sigma2_b) de volta para variância
  
  Q <- length(nos)
  b_0_vet <- sqrt(2 * sigma2) * nos
  pesos_prod <- pesos / sqrt(pi)
  
  ll_total <- 0
  
  for (i in unique(id)) {
    idx <- (id == i)
    y_i <- y[idx]
    X_i <- X[idx, , drop = FALSE]
    S_i <- S[idx, , drop = FALSE]
    n_i <- length(y_i)
    
    eta1_mat <- outer(as.vector(X_i %*% beta), b_0_vet, "+")
    eta2_mat <- matrix(rep(as.vector(S_i %*% tau), Q), nrow = n_i, ncol = Q)
    
    # TRAVAS DE SEGURANÇA EXTREMA PARA MU E SIGMA
    mu_mat    <- pmin(pmax(plogis(eta1_mat), 1e-6), 1 - 1e-6)
    sigma_mat <- pmin(pmax(exp(eta2_mat), 1e-4), 50)
    
    y_mat     <- matrix(rep(y_i, Q), nrow = n_i, ncol = Q)
    
    pdf_vec <- suppressWarnings(PLreg::dPL(as.vector(y_mat), as.vector(mu_mat),
                                           as.vector(sigma_mat), lambda, family = "NO"))
    
    # Se o otimizador testar um ponto absurdo, penaliza com valor alto
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

obter_chutes_iniciais <- function(y, X, S) {
  y_adj <- pmin(pmax(y, 1e-5), 1 - 1e-5)
  z_resp <- qlogis(y_adj)
  
  # Chutes conservadores baseados em modelo linear
  fit_lm <- lm(z_resp ~ X - 1)
  beta_init <- as.vector(coef(fit_lm))
  tau_init <- rep(0, ncol(S))
  log_lambda_init <- log(1.0) # Assume lambda = 1 (simétrico) no começo
  log_sigma2_init <- log(0.1) # Variância pequena no começo
  
  return(c(beta_init, tau_init, log_lambda_init, log_sigma2_init))
}

gerar_dados_PLIA <- function(N, n_i, beta, tau, lambda, sigma2_b) {
  n_total <- N * n_i
  id <- rep(1:N, each = n_i)
  t_ij <- rep(seq(0.1, 1, length.out = n_i), times = N)
  
  X <- cbind(1, t_ij)
  S <- cbind(1, t_ij)
  
  b_0i <- rnorm(N, mean = 0, sd = sqrt(sigma2_b))
  b_0ij <- rep(b_0i, each = n_i)
  
  eta1 <- as.vector(X %*% beta) + b_0ij
  eta2 <- as.vector(S %*% tau)
  
  mu    <- plogis(eta1)
  sigma <- exp(eta2)
  
  y <- PLreg::rPL(n_total, mu, sigma, lambda, family = "NO")
  dados <- data.frame(id = factor(id), t_ij = t_ij, y = y)
  
  return(list(dados = dados, X = X, S = S))
}

# ==============================================================================
# 2. CONFIGURAÇÃO DA SIMULAÇÃO (Sem gradiente analítico)
# ==============================================================================

rodar_simulacao <- function(M, sigma2, n_medidas, N_amostra, n_nos, 
                            beta = c(-1.5, 1.5), tau = c(-0.8, 1.2), 
                            lambda = 1.5, seed = 123) {
  set.seed(seed)
  
  theta_verdadeiro <- c(beta, tau, lambda, sigma2)
  gh <- statmod::gauss.quad(n_nos, kind = "hermite")
  
  resultados <- matrix(NA, nrow = M, ncol = 6)
  colnames(resultados) <- c("beta0", "beta1", "tau0", "tau1", "lambda", "sigma2_b")
  
  sucessos <- 0
  tentativas <- 0
  
  cat("=========================================\n")
  cat(sprintf("  INICIANDO SIMULAÇÃO (N = %d, sigma2 = %.2f)\n", N_amostra, sigma2))
  cat("=========================================\n\n")
  flush.console()
  
  while (sucessos < M) {
    tentativas <- tentativas + 1
    
    sim <- gerar_dados_PLIA(N = N_amostra, n_i = n_medidas, 
                            beta = beta, tau = tau, 
                            lambda = lambda, sigma2_b = sigma2)
    
    theta_init <- obter_chutes_iniciais(sim$dados$y, sim$X, sim$S)
    
    fit_nosso <- tryCatch({
      nlminb(
        start     = theta_init, 
        objective = log_vero_PLIA_vet,
        # AQUI É O PULO DO GATO: GRADIENTE ANALÍTICO REMOVIDO!
        y         = sim$dados$y, 
        X         = sim$X, 
        S         = sim$S, 
        id        = as.numeric(sim$dados$id), 
        nos       = gh$nodes, 
        pesos     = gh$weights,
        lower     = c(rep(-10, ncol(sim$X)), rep(-10, ncol(sim$S)), -3, -5),
        upper     = c(rep(10, ncol(sim$X)), rep(10, ncol(sim$S)),  3,  5),
        control   = list(eval.max = 2000, iter.max = 2000, rel.tol = 1e-5)
      )
    }, error = function(e) return(NULL))
    
    if (!is.null(fit_nosso) && fit_nosso$convergence == 0) {
      
      # Reverte log(lambda) e log(sigma2) para a escala original!
      vetor_nosso <- c(
        fit_nosso$par[1:4],         
        exp(fit_nosso$par[5]),      # lambda = exp(log_lambda)
        exp(fit_nosso$par[6])       # sigma2 = exp(log_sigma2)
      )
      
      sucessos <- sucessos + 1
      resultados[sucessos, ] <- vetor_nosso
      
      cat(sprintf("Réplica %d de %d concluída! (Tentativas: %d)\n", 
                  sucessos, M, tentativas))
      flush.console()
    }
  }
  
  taxa_sucesso <- (sucessos / tentativas) * 100
  limpo_nosso <- na.omit(resultados)
  
  cat("\n=========================================\n")
  cat(sprintf("  Taxa de Convergência: %.2f%%\n", taxa_sucesso))
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
    row.names  = c("beta0", "beta1", "tau0", "tau1", "lambda", "sigma2_b")
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
sim1 <- rodar_simulacao(M = 2000, sigma2 = 0.05, n_medidas = 5, N_amostra = 10,  n_nos = 100)
sim2 <- rodar_simulacao(M = 2000, sigma2 = 0.05, n_medidas = 5, N_amostra = 20,  n_nos = 100)
sim3 <- rodar_simulacao(M = 2000, sigma2 = 0.05, n_medidas = 5, N_amostra = 40,  n_nos = 100)
sim4 <- rodar_simulacao(M = 2000, sigma2 = 0.05, n_medidas = 5, N_amostra = 100, n_nos = 100)
sim5 <- rodar_simulacao(M = 2000, sigma2 = 0.75, n_medidas = 5, N_amostra = 10,  n_nos = 100)
sim6 <- rodar_simulacao(M = 2000, sigma2 = 0.75, n_medidas = 5, N_amostra = 20,  n_nos = 100)
sim7 <- rodar_simulacao(M = 2000, sigma2 = 0.75, n_medidas = 5, N_amostra = 40,  n_nos = 100)
sim8 <- rodar_simulacao(M = 2000, sigma2 = 0.75, n_medidas = 5, N_amostra = 100, n_nos = 100)

print(sim1$tabela_resumo)
print(sim2$tabela_resumo)
print(sim3$tabela_resumo)
print(sim4$tabela_resumo)
print(sim5$tabela_resumo)
print(sim6$tabela_resumo)
print(sim7$tabela_resumo)
print(sim8$tabela_resumo)

print(sim1$taxa_convergencia)
print(sim2$taxa_convergencia)
print(sim3$taxa_convergencia)
print(sim4$taxa_convergencia)
print(sim5$taxa_convergencia)
print(sim6$taxa_convergencia)
print(sim7$taxa_convergencia)
print(sim8$taxa_convergencia)

saveRDS(sim1$limpo_nosso, file = "plno6_1.rds"); plno6_1 <- readRDS("plno6_1.rds")
saveRDS(sim2$limpo_nosso, file = "plno6_2.rds"); plno6_2 <- readRDS("plno6_2.rds")
saveRDS(sim3$limpo_nosso, file = "plno6_3.rds"); plno6_3 <- readRDS("plno6_3.rds")
saveRDS(sim4$limpo_nosso, file = "plno6_4.rds"); plno6_4 <- readRDS("plno6_4.rds")
saveRDS(sim5$limpo_nosso, file = "plno6_5.rds"); plno6_5 <- readRDS("plno6_5.rds")
saveRDS(sim6$limpo_nosso, file = "plno6_6.rds"); plno6_6 <- readRDS("plno6_6.rds")
saveRDS(sim7$limpo_nosso, file = "plno6_7.rds"); plno6_7 <- readRDS("plno6_7.rds")
saveRDS(sim8$limpo_nosso, file = "plno6_8.rds"); plno6_8 <- readRDS("plno6_8.rds")
