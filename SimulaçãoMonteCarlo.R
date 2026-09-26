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
        intervalo_inicial <- c(0.1, 10)
        while (TRUE) {
            estimativa_alpha <- try(
                estimar_alpha(amostra, intervalo = intervalo_inicial),
                silent = TRUE
            )
            ifelse(
                is.numeric(estimativa_alpha),
                break,
                intervalo_inicial <- intervalo_inicial * c(1/2, 2)
            )
        }

        #### Estimar o parâmetro beta
        estimativa_beta <- mean(amostra) / estimativa_alpha

        #### Adicionar a r-ésima amostra da tabela amostras_estimativas
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

### Salvar a tabela em um arquivo CSV
avaliacao_estimadores %>%
    write.csv(file = 'DadosAvaliaçãoEstimadores.csv', row.names = FALSE)

### Gerar os histogramas das estimativas do parâmetro alpha para cada tamanho amostral

#### Organizar os histogramas em 2 por 3
layout(matrix(c(
  1, 1, 2, 2, 3, 3,
  0, 4, 4, 5, 5, 0
), 2, 6, byrow = TRUE))

for (t in tamanhos_amostrais) {
    #### Filtrar as amostras e estimativas por tamanho amostral t
    n_amostras_estimativas <- amostras_estimativas %>% filter(n == t)

    #### Gerar os histogramas por tamanho amostral t
    hist(
        x = n_amostras_estimativas$est_alpha,
        main = paste0('De amostras com ',as.character(t), ' observações'),
        col = '#d8caa3',
        xlab = '',
        ylab = '',
        freq = FALSE
    )

    #### Calcular variância teórica da distribuição do estimador de alpha
    var_teo_est_alpha <- 1 / (t * trigamma(alpha) - t / alpha)

    #### Adicionar a curva da f.d.p. da distribuição assintótica
    curve(
        dnorm(x, mean = alpha, sd = sqrt(var_teo_est_alpha)),
        col = 'red',
        lty = 2,
        add = TRUE
    )
}

### Gerar os histogramas das estimativas do parâmetro beta para cada tamanho amostral

#### Organizar os histogramas em 2 por 3
layout(matrix(c(
  1, 1, 2, 2, 3, 3,
  0, 4, 4, 5, 5, 0
), 2, 6, byrow = TRUE))

for (t in tamanhos_amostrais) {
    #### Filtrar as amostras e estimativas por tamanho amostral t
    n_amostras_estimativas <- amostras_estimativas %>% filter(n == t)

    #### Gerar os histogramas por tamanho amostral t
    hist(
        x = n_amostras_estimativas$est_beta,
        main = paste0('De amostras com ',as.character(t), ' observações'),
        col = '#d8caa3',
        xlab = '',
        ylab = '',
        freq = FALSE
    )

    #### Calcular variância teórica da distribuição do estimador de beta
    var_teo_est_beta <- (trigamma(alpha) * beta**2) /
        (t * alpha * trigamma(alpha) - t)

    #### Adicionar a curva da f.d.p. da distribuição assintótica
    curve(
        dnorm(x, mean = beta, sd = sqrt(var_teo_est_beta)),
        col = 'red',
        lty = 2,
        add = TRUE
    )
}
