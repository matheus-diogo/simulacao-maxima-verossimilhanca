# Simulação de Monte Carlo

## Importação dos Pacotes Necessários ############################################################################
library(tibble)
library(dplyr)
library(magrittr)

## Definição do Método Numérico para Estimação ###################################################################
#
# Em progresso...
#

## Geração de Diversas Amostras a partir da Distribuição Gama ####################################################

### Fixar a semente para reprodutibilidade
set.seed(2026)

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
        #
        # Em progresso...
        #
        estimativa_alpha <- NA

        #### Estimar o parâmetro beta pelo método numérico
        #
        # Em progresso...
        #
        estimativa_beta <- NA

        #### Adicionar a r-ésima amostra na tabela criada
        amostras <- amostras %>%
            add_row(
                n = tamanho_amostral,
                est_alpha = estimativa_alpha,
                est_beta = estimativa_beta
            )
    }
}

## Avaliação dos Estimadores #####################################################################################
#
# Em progresso...
#