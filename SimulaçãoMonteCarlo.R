# Simulação de Monte Carlo

## Importação dos Pacotes Necessários ############################################################################
library(tibble)
library(dplyr)
library(magrittr)

## Definição do Método Numérico para Estimação ###################################################################

### Definir uma função em parâmetro alpha, dada uma amostra, para aplicar o processo iterativo
f_alpha <- function(alpha, amostra) {
    log(alpha) - digamma(alpha) - log(mean(amostra)) + mean(log(amostra))
}

### Determinar a estimativa do parâmetro alpha pelo método de Brent aplicado a função f_alpha
estimar_alpha <- function(amostra, intervalo) {
    estimativa <- uniroot(
        f = f_alpha,
        interval = intervalo,
        check.conv = TRUE,
        amostra = amostra
    )
    return(estimativa$root)
}

#### A estimativa do parâmetro beta é expressa como a média amostral dividida pela estimativa do parâmetro alpha

## Geração de Diversas Amostras a partir da Distribuição Gama ####################################################

### Fixar a semente para reprodutibilidade
### set.seed(2026)

### Selecionar os verdadeiros valores dos parâmetros da distribuição Gama
alpha <- 4
beta <- 3

### Determinar todos os tamanhos amostrais utilizados
tamanhos_amostrais <- c(20, 30, 50, 100, 200)

### Definir uma tabela das amostras com seus tamanhos amostrais e estimativas dos parâmetros
amostras <- tibble(
    n = numeric(),
    est_alpha = numeric(),
    est_beta = numeric()
)

### Gerar 10000 amostras para cada tamanho amostral
for (t in tamanhos_amostrais) {
    for (r in 1:10000) {
        #### Criar r-ésima amostra de tamanho amostral t
        amostra <- rgamma(n = t, shape = alpha, scale = beta)

        #### Salvar o tamanho amostral utilizado
        tamanho_amostral <- t

        #### Estimar o parâmetro alpha pelo método numérico
        intervalo_inicial <- c(1, 10)
        while (TRUE) {
            estimativa_alpha <- try(
                estimar_alpha(amostra, intervalo = intervalo_inicial),
                silent = TRUE
            )
            ifelse(
                is.numeric(estimativa_alpha),
                break,
                intervalo_inicial <- intervalo_inicial * c(0.5, 2)
            )
        }

        #### Estimar o parâmetro beta
        estimativa_beta <- mean(amostra) / estimativa_alpha

        #### Adicionar a r-ésima amostra na tabela criada
        amostras <- amostras %>%
            add_row(
                n = tamanho_amostral,
                est_alpha = estimativa_alpha,
                est_beta = estimativa_beta
            )
    }
}

amostras %>%
    summarise(med_est_alpha = mean(est_alpha), med_est_beta = mean(est_beta))

## Avaliação dos Estimadores #####################################################################################
#
# Em progresso...
#
