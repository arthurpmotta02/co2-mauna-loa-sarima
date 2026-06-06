# CO₂ Atmosférico em Mauna Loa — Análise SARIMA

**Trabalho 1 — Análise de Séries Temporais**  
Arthur Motta & Catarine Martins

📄 **[Ver relatório completo](https://arthurpmotta02.github.io/co2-mauna-loa-sarima/)**

---

## Descrição

Análise completa da série temporal `co2` (disponível nativamente no R), contendo as concentrações médias mensais de CO₂ atmosférico medidas pelo **Observatório de Mauna Loa**, Havaí, de janeiro de 1959 a dezembro de 1997 (468 observações).

A série é conhecida como **Curva de Keeling** e é um caso de referência clássico na literatura de séries temporais por combinar uma tendência crescente de longo prazo com sazonalidade anual regular — estrutura ideal para a metodologia SARIMA.

---

## Estrutura do Repositório

```
.
├── .github/
│   └── workflows/
│       └── publish.yml                 # GitHub Actions: renderiza e publica automaticamente
├── Trabalho1_SeriesTemporais_CO2.qmd   # Documento principal (Quarto)
├── references.bib                      # Referências bibliográficas
└── README.md
```

> O relatório HTML é gerado automaticamente pelo GitHub Actions a cada `push` na branch `main` e publicado em [arthurpmotta02.github.io/co2-mauna-loa-sarima](https://arthurpmotta02.github.io/co2-mauna-loa-sarima/). Não há arquivo HTML versionado no repositório — ele é sempre produzido a partir do `.qmd`.

---

## Metodologia

A análise segue quatro etapas:

1. **Análise descritiva** — visualização da série, decomposição STL, padrão sazonal e taxa de crescimento anual
2. **Identificação e ajuste** — testes de estacionariedade (ADF/KPSS), leitura de ACF/PACF, comparação de 7 modelos SARIMA candidatos por AIC, AICc e BIC
3. **Diagnóstico** — análise dos resíduos (autocorrelação via Ljung-Box, normalidade via Shapiro-Wilk e Jarque-Bera, homoscedasticidade visual)
4. **Previsão** — projeções para 1998–1999 com intervalos de confiança de 80% e 95%

Como complemento, ajusta-se também um **modelo de regressão com erros ARMA** (tendência linear + harmônicos de Fourier), comparando com o SARIMA por critérios de informação.

---

## Modelo Selecionado

**SARIMA(1,1,1)(0,1,1)₁₂**

Selecionado por parcimônia diante de empate técnico (ΔAICc < 0,3) entre os três melhores candidatos. Com apenas 3 parâmetros livres, oferece capacidade preditiva equivalente a modelos mais complexos.

| Parâmetro | Estimativa | p-valor |
|-----------|-----------|---------|
| φ₁ (AR não-sazonal) | 0.2399 | 0.094 |
| θ₁ (MA não-sazonal) | −0.5710 | < 0.001 |
| θ̃₁ (SMA sazonal) | −0.8516 | < 0.001 |

**Diagnósticos:** resíduos sem autocorrelação (Ljung-Box p > 0,10 para h ≤ 48), normalidade não rejeitada (Shapiro-Wilk p = 0,532), modelo causal e invertível.

---

## Como Reproduzir Localmente

### Pré-requisitos

R (≥ 4.1) e [Quarto](https://quarto.org/) instalados, com os seguintes pacotes:

```r
install.packages(c(
  "forecast", "tseries", "TSA",
  "ggplot2", "patchwork", "scales", "ggfortify",
  "gt", "broom", "tidyverse"
))
```

### Renderizar o relatório

```bash
quarto render Trabalho1_SeriesTemporais_CO2.qmd
```

O arquivo `Trabalho1_SeriesTemporais_CO2.html` será gerado na mesma pasta — autocontido (`embed-resources: true`), sem dependências externas.

---

## Principais Resultados

- Taxa de crescimento médio: **~1,23 ppm/ano** ao longo de 39 anos
- Sazonalidade anual de amplitude **~6 ppm** (pico em maio, vale em setembro)
- Previsão para dezembro/1998: **365,60 ppm** | dezembro/1999: **367,14 ppm** (IC 95%: 365,34–368,94 ppm)
- Incerteza de previsão: de ±0,5 ppm (h = 1) a ±2,0 ppm (h = 24) — erro relativo < 0,6%

---

## Referências

- Brockwell, P.J. & Davis, R.A. (2010). *Introduction to Time Series and Forecasting*. Springer.
- Morettin, P.A. & Toloi, C.M.C. (2004). *Análise de Séries Temporais*. Edgard Blücher.
- Box, G.E.P. et al. (2015). *Time Series Analysis: Forecasting and Control*. Wiley.
- Keeling, C.D. et al. (1976). Atmospheric carbon dioxide variations at Mauna Loa Observatory. *Tellus*, 28(6), 538–551.
- Cleveland, R.B. et al. (1990). STL: A seasonal-trend decomposition procedure based on LOESS. *Journal of Official Statistics*, 6(1), 3–73.
- Burnham, K.P. & Anderson, D.R. (2002). *Model Selection and Multimodel Inference*. Springer.