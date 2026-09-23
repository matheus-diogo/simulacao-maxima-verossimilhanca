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

### Selecionar os verdadeiros valores dos parâmetros da distribuição Gama
alpha <- 4
beta <- 3

### Determinar todos os tamanhos amostrais utilizados e quantidade R de amostras para cada tamanho amostral
tamanhos_amostrais <- c(20, 30, 50, 100, 200)
R <- 10000

### Definir uma tabela das amostras com seus tamanhos amostrais e estimativas dos parâmetros
amostras_estimativas <- tibble(
    n = numeric(),
    est_alpha = numeric(),
    est_beta = numeric()
)

### Fixar a semente para garantir reprodutibilidade
set.seed(2026)

### Gerar 10000 amostras para cada tamanho amostral
for (t in tamanhos_amostrais) {
    for (r in 1:R) {
        #### Criar r-ésima amostra de tamanho amostral t
        amostra <- rgamma(n = t, shape = alpha, scale = beta)

        #### Salvar o tamanho amostral utilizado
        tamanho_amostral <- t

        #### Estimar o parâmetro alpha pelo método numérico
        intervalo_inicial <- c(1, 1)
        while (TRUE) {
            estimativa_alpha <- try(
                estimar_alpha(amostra, intervalo = intervalo_inicial),
                silent = TRUE
            )
            ifelse(
                is.numeric(estimativa_alpha),
                break,
                intervalo_inicial <- intervalo_inicial * c(1 / 2, 2)
            )
        }

        #### Estimar o parâmetro beta
        estimativa_beta <- mean(amostra) / estimativa_alpha

        #### Adicionar a r-ésima amostra na tabela criada
        amostras_estimativas <- amostras_estimativas %>%
            add_row(
                n = tamanho_amostral,
                est_alpha = estimativa_alpha,
                est_beta = estimativa_beta
            )
    }
}

## Avaliação dos Estimadores #####################################################################################

### Gerar a tabela com a avaliação dos estimadores para cada tamanho amostral
for (t in tamanhos_amostrais) {
    #### Calcular a média das estimativas
    avaliacao_parcial_estimadores <- amostras_estimativas %>%
        filter(n == t) %>%
        summarise(
            tamanho_amostral = t,
            med_est_alpha = mean(est_alpha),
            med_est_beta = mean(est_beta)
        )

    #### Calcular o viés
    avaliacao_parcial_estimadores <- avaliacao_parcial_estimadores %>%
        mutate(
            vies_est_alpha = med_est_alpha - alpha,
            vies_est_beta = med_est_beta - beta
        )

    #### Calcular o viés relativo percentual
    avaliacao_parcial_estimadores <- avaliacao_parcial_estimadores %>%
        mutate(
            vies_rel_est_alpha = (100 / alpha) * vies_est_alpha,
            vies_rel_est_beta = (100 / beta) * vies_est_beta
        )

    #### Calcular a variância Monte Carlo
    var_monte_carlo <- amostras_estimativas %>%
        filter(n == t) %>%
        mutate(
            var_mc_est_alpha = var(est_alpha - mean(est_alpha)),
            var_mc_est_beta = var(est_beta - mean(est_beta))
        ) %>%
        select(var_mc_est_alpha, var_mc_est_beta) %>%
        slice(1)

    avaliacao_parcial_estimadores <- bind_cols(
        avaliacao_parcial_estimadores,
        var_monte_carlo
    )

    #### Calcular o erro quadrático médio
    erro_quad_medio <- amostras_estimativas %>%
        filter(n == t) %>%
        mutate(
            eqm_est_alpha = var(est_alpha - alpha) * ((t) / (t - 1)),
            eqm_est_beta = var(est_beta - beta) * ((t) / (t - 1))
        ) %>%
        select(eqm_est_alpha, eqm_est_beta) %>%
        slice(1)

    avaliacao_parcial_estimadores <- bind_cols(
        avaliacao_parcial_estimadores,
        erro_quad_medio
    )

    if (t == 20) {
        avaliacao_estimadores <- avaliacao_parcial_estimadores
    } else {
        avaliacao_estimadores <- bind_rows(
            avaliacao_estimadores,
            avaliacao_parcial_estimadores
        )
    }
}

### Exibir a tabela com a avaliação dos estimadores
avaliacao_estimadores
