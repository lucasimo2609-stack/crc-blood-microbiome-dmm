# =========================================================
# ANALISI MICROBIOMA CON DIRICHLET MULTINOMIAL MIXTURE
# =========================================================

# =========================================================
# INSTALLAZIONE PACCHETTI
# =========================================================

if (!require("haven")) install.packages("haven")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("DirichletMultinomial"))
  BiocManager::install("DirichletMultinomial")

if (!require("ggplot2"))
  install.packages("ggplot2")

if (!require("lattice"))
  install.packages("lattice")

# =========================================================
# LIBRERIE
# =========================================================

library(haven)
library(DirichletMultinomial)
library(ggplot2)
library(lattice)

# =========================================================
# IMPORT FILE SAS
# =========================================================

df_microbioma <- read_sas(file.choose())

# Controllo struttura
str(df_microbioma)

# Prime righe
head(df_microbioma)

# =========================================================
# PREPARAZIONE MATRICE COUNT
# =========================================================

# Prima colonna = ID pazienti
row_names_pazienti <- df_microbioma[[1]]

# Rimuove colonna ID
count_matrix <- df_microbioma[, -1]

# Converte tutto in numerico
count_matrix <- apply(
  count_matrix,
  2,
  as.numeric
)

# Converte in matrice
count_matrix <- as.matrix(count_matrix)

# Assegna rownames
rownames(count_matrix) <- row_names_pazienti

# Sostituisce eventuali NA
count_matrix[is.na(count_matrix)] <- 0

# =========================================================
# FILTRO TAXA RARI
# =========================================================

# evito conti inutili che tanto non mi danno prevalenza 
keep_taxa <- colSums(count_matrix) > 10

count_matrix <- count_matrix[, keep_taxa]

# =========================================================
# CONVERSIONE A INTERI
# =========================================================

storage.mode(count_matrix) <- "integer"

# =========================================================
# CONTROLLO DATI
# =========================================================

cat("\nDimensioni matrice:\n")
print(dim(count_matrix))

cat("\nPrime colonne:\n")
print(head(count_matrix[, 1:5]))

# =========================================================
# DISTRIBUZIONE TAXA
# =========================================================

# cnts <- log10(colSums(count_matrix) + 1)

# densityplot(
#  cnts,
#  xlab = "Taxon representation (log10 count)",
#  main = "Distribuzione abbondanza taxa"
# )

# =========================================================
# FIT DMM SOLO K = 1 E K = 2
# =========================================================

set.seed(123)

fit <- lapply(
  1:2,
  dmn,
  count = count_matrix,
  verbose = TRUE
)

# =========================================================
# VALUTAZIONE MODELLI
# =========================================================

laplace_values <- sapply(
  fit,
  laplace
)

aic_values <- sapply(
  fit,
  AIC
)

bic_values <- sapply(
  fit,
  BIC
)

model_eval <- data.frame(
  K = 1:2,
  Laplace = laplace_values,
  AIC = aic_values,
  BIC = bic_values
)

cat("\nValutazione modelli:\n")

print(model_eval)

# =========================================================
# GRAFICO LAPPLACE
# =========================================================

plot(
  1:2,
  laplace_values,
  type = "b",
  pch = 19,
  xlab = "Numero cluster",
  ylab = "Laplace",
  main = "Confronto modelli DMM"
)

# =========================================================
# USA SEMPRE K = 2
# =========================================================

fit_k2 <- fit[[2]]

cat("\nModello k=2:\n")
print(fit_k2)

# =========================================================
# PESI CLUSTER
# =========================================================

mixture_weights <- mixturewt(fit_k2)

cat("\nPesi cluster:\n")
print(mixture_weights)

# =========================================================
# PROBABILITA POSTERIORI
# =========================================================

posterior_prob <- mixture(fit_k2)

cat("\nProbabilità posteriori:\n")
print(head(posterior_prob))

# =========================================================
# ASSEGNAZIONE CLUSTER
# =========================================================

clusters <- apply(
  posterior_prob,
  1,
  which.max
)

cat("\nNumero pazienti per cluster:\n")
print(table(clusters))

# =========================================================
# DATAFRAME RISULTATI
# =========================================================

risultati_cluster <- data.frame(
  Paziente = rownames(count_matrix),
  Cluster = clusters
)

cat("\nPrime assegnazioni:\n")
print(head(risultati_cluster))

# =========================================================
# COMPOSIZIONE TAXA
# =========================================================

taxa_fit <- fitted(
  fit_k2,
  scale = FALSE
)

cat("\nPrime frequenze taxa:\n")
print(head(taxa_fit))


# =========================================================
# TAXA DIFFERENZIALI
# =========================================================

cluster1 <- taxa_fit[,1]

cluster2 <- taxa_fit[,2]

diff_taxa <- abs(
  cluster1 - cluster2
)

ord <- order(
  diff_taxa,
  decreasing = TRUE
)

top_taxa <- data.frame(
  Taxa = rownames(taxa_fit)[ord],
  Cluster1 = cluster1[ord],
  Cluster2 = cluster2[ord],
  Difference = diff_taxa[ord]
)

cat("\nTop 20 taxa discriminanti:\n")
print(head(top_taxa, 20))

# =========================================================
# HEATMAP
# =========================================================

# cat("\nGenerazione heatmap...\n")

# heatmapdmn(
#   count = count_matrix,
#   fit1 = fit[[1]],
#   fitN = fit_k2,
#   ntaxa = 30
# )

# =========================================================
# MESSAGGIO FINALE
# =========================================================

cat("\n====================================\n")
cat("ANALISI COMPLETATA\n")
cat("====================================\n")

cat("\nFile esportati:\n")
cat("- cluster_pazienti_k2.csv\n")
cat("- top_taxa_k2.csv\n")
cat("- valutazione_modelli_DMM.csv\n")

# =========================================================
# CARICAMENTO DEI DATI CLINICI E INCROCIO CON I CLUSTER DMM
# =========================================================

cat("\nCaricamento del dataset clinico MTCR USA...\n")

# 1. Importa il secondo dataset (cambia l'estensione se è .sas7bdat, .csv o .xlsx)
df_clinico <- read_sas(file.choose())

# Controlliamo al volo i nomi delle colonne del file clinico per essere sicuri di quale sia l'ID
cat("\nColonne del dataset clinico:\n")
print(colnames(df_clinico))

# 2. Prepariamo il dataframe con i cluster ottenuti dal DMM
# Abbiniamo l'ID 'v1' al Cluster calcolato (1 o 2)
df_cluster_risultati <- data.frame(
  v1 = df_microbioma$v1,
  Cluster_DMM = as.factor(clusters)
)

# 3. Uniamo i due dataset usando come chiave l'ID del paziente
# NOTA: se nel file clinico l'ID non si chiama 'v1' ma ad esempio 'ID', 
# sostituisci con: by.x = "v1", by.y = "ID"
df_completo <- merge(df_cluster_risultati, df_clinico, by = "v1")

# Recodifichiamo la variabile V2 per renderla leggibile nei grafici e nelle tabelle
df_completo$v2 <- factor(df_completo$v2, 
                                   levels = c(1, 2, 3), 
                                   labels = c("Sano", "Adenoma", "Cancro"))

# =========================================================
# TABELLA DI CONFRONTO E STATISTICA SIGNIFICATIVA
# =========================================================

# 4. Creiamo la tabella di contingenza incrociando i 2 Cluster con le 3 classi cliniche (V2)
tabella_confronto_reale <- table(df_completo$Cluster_DMM, df_completo$v2)

cat("\nTabella di confronto REALE (Cluster vs Stato Patologico):\n")
print(tabella_confronto_reale)

# 5. Test del Chi-quadrato per verificare la significatività statistica dell'associazione
test_chi2 <- chisq.test(tabella_confronto_reale)

cat("\nRisultato del Test del Chi-quadrato:\n")
print(test_chi2)

# =========================================================
# GRAFICO DI PRESENTAZIONE
# =========================================================

# Generiamo un grafico a barre impilate per mostrare visivamente la scomposizione clinica nei due cluster
barplot(tabella_confronto_reale, 
        main="Distribuzione degli Stati Patologici nei Cluster DMM",
        xlab="Stato di Salute (V2)", 
        col=c("#4287f5", "#f54242"), # Blu per Cluster 1, Rosso per Cluster 2
        legend=paste("Cluster", rownames(tabella_confronto_reale)), 
        beside=FALSE)

# 6. Esportazione finale del file unito per eventuali analisi successive (anche in SAS)
write.csv(
  df_completo,
  "cluster_dmm_con_clinica_completa.csv",
  row.names = FALSE
)

cat("\nAnalisi clinica completata con successo!\n")

# =========================================================
# ANALISI 2: VARIABILE Cancer (Dicotomica 0 / 1)
# =========================================================
cat("\n--- Analisi rispetto a 'Cancer' ---\n")

# Convertiamo in fattore la variabile Cancer
df_completo$Cancer_Fattore <- as.factor(df_completo$cancer)

# Tabella di contingenza per Cancer
tabella_confronto_cancer <- table(df_completo$Cluster_DMM, df_completo$Cancer_Fattore)
cat("\nTabella di confronto REALE (Cluster vs Cancer):\n")
print(tabella_confronto_cancer)

# Test del Chi-quadrato per Cancer
test_chi2_cancer <- chisq.test(tabella_confronto_cancer)
cat("\nRisultato del Test del Chi-quadrato (Cancer):\n")
print(test_chi2_cancer)


# =========================================================
# GRAFICI DI PRESENTAZIONE (Affiancati)
# =========================================================
# Configura R per mostrare due grafici nella stessa riga
par(mfrow = c(1, 2)) 

# Grafico v2
barplot(tabella_confronto_reale, 
        main="Progressione Clinica (v2)",
        xlab="Stato di Salute", 
        col=c("#4287f5", "#f54242"), 
        legend=paste("Cluster", rownames(tabella_confronto_reale)), 
        beside=FALSE)

# Grafico Cancer
barplot(tabella_confronto_cancer, 
        main="Presenza Tumore (Cancer)",
        xlab="Presenza di Cancro", 
        col=c("#4287f5", "#f54242"), 
        beside=FALSE)

# Ripristina il layout del grafico singolo
par(mfrow = c(1, 1)) 


# =========================================================
# ESPORTAZIONE FINALE DEL FILE COMPLETO
# =========================================================
write.csv(
  df_completo,
  "cluster_dmm_con_clinica_completa.csv",
  row.names = FALSE
)

cat("\n====================================\n")
cat("TUTTE LE ANALISI CLINICHE COMPLETATE CON SUCCESSO!\n")
cat("====================================\n")


# =========================================================
# REGRESSIONE LOGISTICA MULTIVARIATA PER LA VARIABILE 'cancer'
# =========================================================

cat("\n====================================\n")
cat("MODELLO DI REGRESSIONE LOGISTICA\n")
cat("====================================\n")

# 1. Pulizia e preparazione rapida delle covariate
# Ci assicuriamo che la variabile target sia numerica (0/1) per la regressione logistica
df_completo$cancer_num <- as.numeric(as.character(df_completo$cancer))

df_completo$age_catt1=as.factor(df_completo$age_catt)
levels(df_completo$age_catt1)

#-------------------------------------------------------------------------------
# 2. Fit del modello di regressione logistica
# Modellizziamo la probabilità di avere il cancro (cancer_num = 1)
modello_logistico <- glm(
  cancer_num ~ Cluster_DMM + centro + v4 + edu_cat_an + age_catt, 
  data = df_completo, 
  family = binomial(link = "logit")
)
#-------------------------------------------------------------------------------

# 2. Fit del modello di regressione logistica
# Modellizziamo la probabilità di avere il cancro (cancer_num = 1)
modello_logistico <- glm(
  cancer_num ~ Cluster_DMM + centro + v4 + edu_cat_an + age_catt1, 
  data = df_completo, 
  family = binomial(link = "logit")
)

# 3. Visualizzazione dei risultati standard (Coeff, Errori Standard, p-value)
cat("\n--- Sintesi del Modello (Summary) ---\n")
print(summary(modello_logistico))

# 4. Calcolo degli Odds Ratio (OR) e degli intervalli di confidenza al 95%
cat("\n--- Odds Ratio (OR) e Intervalli di Confidenza (95%) ---\n")

or_valori <- exp(coef(modello_logistico))
ci_valori <- exp(confint(modello_logistico))

# Uniamo tutto in una tabella leggibile
tabella_or <- cbind(Odds_Ratio = or_valori, ci_valori)
print(tabella_or)


cat("\n====================================\n")
cat("REGRESSIONE LOGISTICA COMPLETATA!\n")
cat("====================================\n")












# ==============================================================================
# ADATTAMENTO DEL MODELLO DI HOLMES (DMM) PER TUTTE LE VARIABILI CLINICHE
# ==============================================================================

cat("\n==================================================\n")
cat("PREPARAZIONE VARIABILI CLINICHE PER IL MODELLO DMM\n")
cat("==================================================\n")

# 1. Selezione automatica delle variabili cliniche utili (solo numeriche)
colonne_da_escludere <- c("v1", "Cluster_DMM", "ID", "v2", "cancer", "Cancer_Fattore", "cancer_num", "lesione")
tutte_le_altre <- setdiff(colnames(df_completo), colonne_da_escludere)

# Estraiamo solo le colonne che contengono numeri
colonne_numeriche <- sapply(df_completo[, tutte_le_altre], is.numeric)
df_clinico_numerico <- df_completo[, tutte_le_altre][, colonne_numeriche]

# Filtro di qualità: rimuoviamo colonne con troppi NA (es. > 40%) o varianza zero (tutti valori identici)
colonne_da_tenere <- sapply(df_clinico_numerico, function(x) {
  pct_na <- sum(is.na(x)) / length(x)
  dev_standard <- sd(x, na.rm = TRUE)
  return(pct_na < 0.4 && !is.na(dev_standard) && dev_standard > 0)
})

df_clinico_pulito <- df_clinico_numerico[, colonne_da_tenere]
cat("Numero di variabili cliniche valide selezionate per il modello:", ncol(df_clinico_pulito), "\n")


# 2. Gestione dei dati mancanti (Imputazione tramite Mediana)
# DMM non tollera i valori NA. Imputiamo la mediana per mantenere tutti i 300 pazienti a bordo.
for(colonna in colnames(df_clinico_pulito)) {
  if(any(is.na(df_clinico_pulito[[colonna]]))) {
    mediana_val <- median(df_clinico_pulito[[colonna]], na.rm = TRUE)
    df_clinico_pulito[[colonna]][is.na(df_clinico_pulito[[colonna]])] <- mediana_val
  }
}


# 3. Funzione di Normalizzazione e Discretizzazione Equilibrata
# Trasforma le variabili in un range comune (1-50) di numeri interi per dare a tutte lo stesso peso
normalizza_per_dmm <- function(x, nuovo_min = 1, nuovo_max = 50) {
  x_scaled <- (x - min(x)) / (max(x) - min(x)) # Min-Max standard
  x_intero <- round(x_scaled * (nuovo_max - nuovo_min) + nuovo_min) # Proiezione su range intero
  return(as.integer(x_intero))
}

# Creazione della matrice finale di "conteggi clinici"
matrice_clinica_holmes <- apply(df_clinico_pulito, 2, normalizza_per_dmm)
rownames(matrice_clinica_holmes) <- df_completo$v1

# Controllo di sicurezza prima del fit
cat("Verifica matrice: Ci sono NA?", any(is.na(matrice_clinica_holmes)), 
    "| Valori negativi?", any(matrice_clinica_holmes < 0), "\n")


# 4. Esecuzione del Modello di Holmes (DMM) sulle Variabili Cliniche
library(DirichletMultinomial)

cat("\nEsecuzione del modello DMM clinico (K = 1 e K = 2)...\n")
fit_clinico_k1 <- dmn(matrice_clinica_holmes, k = 1)
fit_clinico_k2 <- dmn(matrice_clinica_holmes, k = 2)

# Visualizzazione delle metriche di fit (AIC, BIC, Laplace)
cat("\n--- Fit Modello Clinico K=1 ---\n")
print(fit_clinico_k1)
cat("\n--- Fit Modello Clinico K=2 ---\n")
print(fit_clinico_k2)

cluster_clinici_assegnati <- mixture(fit_clinico_k2, assign = TRUE)

# Riga corretta senza l'argomento fantasma:
df_completo$Cluster_Clinico_Holmes <- as.factor(cluster_clinici_assegnati)


cat("\n==================================================\n")
cat("CONFRONTO TRA CLUSTER CLINICI E PRESENZA DI CANCRO\n")
cat("==================================================\n")

tabella_clinica_cancer <- table(df_completo$Cluster_Clinico_Holmes, as.factor(df_completo$cancer))
cat("\nTabella di contingenza (Cluster Clinici vs Cancer):\n")
print(tabella_clinica_cancer)

test_chi2_clinico <- chisq.test(tabella_clinica_cancer)
cat("\nRisultato del Test del Chi-quadrato:\n")
print(test_chi2_clinico)

# Sovrascrittura finale del file CSV con tutte le elaborazioni incluse
write.csv(df_completo, "cluster_dmm_con_clinica_completa.csv", row.names = FALSE)
cat("\nAnalisi conclusa. I nuovi gruppi clinici sono pronti e salvati!\n")

aggregate(df_completo[, c("age", "bmi")], by = list(Cluster = df_completo$Cluster_Clinico_Holmes), FUN = mean, na.rm = TRUE)


################################################################################

# ==============================================================================
# SCRIPT COMPLETO: MODELLO DMM IBRIDO (MATRICE BATTERI + VARIABILI CLINICHE)
# K = 2 CLUSTER & VERIFICA ASSOCIAZIONE CON CANCER
# ==============================================================================

cat("\n==================================================\n")
cat("1. PREPARAZIONE E SELEZIONE DELLE VARIABILI CLINICHE\n")
cat("==================================================\n")

# A. Identificazione dinamica del nome della variabile età (evita errori tra Age/age)
nome_eta <- if("Age" %in% colnames(df_clinico)) "Age" else "age"
target_vars <- c("bmi", nome_eta, "v12")

cat("Variabili cliniche rilevate nel dataset:", paste(target_vars, collapse=", "), "\n")

# B. Estrazione delle colonne target dal dataframe principale
df_estrazione <- df_clinico[, target_vars]

# C. Controllo e conversione forzata di v12 a numerica (necessario per la matrice del DMM)
if(!is.numeric(df_estrazione$v12)) {
  cat("Nota: 'v12' è una variabile categoriale/testo. Viene convertita in numerica.\n")
  df_estrazione$v12 <- as.numeric(as.factor(df_estrazione$v12))
}

cat("\n==================================================\n")
cat("2. IMPUTAZIONE DEI DATI MANCANTI (NA)\n")
cat("==================================================\n")

# Sostituzione dei valori NA con la mediana della rispettiva colonna
for(col in colnames(df_estrazione)) {
  if(any(is.na(df_estrazione[[col]]))) {
    mediana_v <- median(df_estrazione[[col]], na.rm = TRUE)
    df_estrazione[[col]][is.na(df_estrazione[[col]])] <- mediana_v
    cat("Imputati valori mancanti per la variabile:", col, "usando la mediana (", mediana_v, ")\n")
  }
}

cat("\n==================================================\n")
cat("3. RISCALAMENTO PROPORZIONALE (PESO CLINICO)\n")
cat("==================================================\n")

# Funzione per proiettare le variabili cliniche su un range 1-500,
# evitando che vengano matematicamente "annegate" dai grandi conteggi dei batteri
normalizza_ibrido <- function(x, nuovo_min = 1, nuovo_max = 500) {
  if(max(x) == min(x)) return(rep(as.integer(nuovo_min), length(x)))
  x_scaled <- (x - min(x)) / (max(x) - min(x))
  x_intero <- round(x_scaled * (nuovo_max - nuovo_min) + nuovo_min)
  return(as.integer(x_intero))
}

matrice_clinica_scaled <- apply(df_estrazione, 2, normalizza_ibrido)

cat("\n==================================================\n")
cat("4. UNIONE SICURA (MERGE) TRAMITE ID PAZIENTI\n")
cat("==================================================\n")

# Trasformiamo i dati clinici scalati in dataframe e assegniamo gli ID presenti in df_clinico$v1
df_clinico_fissato <- data.frame(matrice_clinica_scaled)
rownames(df_clinico_fissato) <- df_clinico$v1

# Convertiamo temporaneamente la matrice dei batteri in dataframe per effettuare il merge
df_batteri <- data.frame(count_matrix)

# Unione formale per riga (by = 'row.names'): R accoppia i dati solo se gli ID coincidono
df_ibrido_join <- merge(df_batteri, df_clinico_fissato, by = "row.names")

# Ripristiniamo la struttura a matrice pura richiesta dal pacchetto DirichletMultinomial
rownames(df_ibrido_join) <- df_ibrido_join$Row.names
matrice_ibrida_holmes <- as.matrix(df_ibrido_join[, -1]) # Escludiamo la colonna di testo degli ID

cat("Merge completato! Nuove dimensioni matrice mista:", dim(matrice_ibrida_holmes)[1], "x", dim(matrice_ibrida_holmes)[2], "\n")

cat("\n==================================================\n")
cat("5. ESECUZIONE ALGORITMO DI HOLMES (DMM) PER K = 2\n")
cat("==================================================\n")

library(DirichletMultinomial)
set.seed(123) # Garantisce la riproducibilità esatta dei cluster ad ogni lancio

fit_ibrido_k2 <- dmn(matrice_ibrida_holmes, k = 2)
print(fit_ibrido_k2)

# Estrazione dell'assegnazione dei cluster per ogni paziente
cluster_ibridi_assegnati <- mixture(fit_ibrido_k2, assign = TRUE)

cat("\n==================================================\n")
cat("6. RIALLINEAMENTO CLUSTER E DIAGNOSTICA FINALE\n")
cat("==================================================\n")

# Mappatura dei cluster dentro il dataframe ORIGINALE (df_clinico)
# Il comando match garantisce che ad ogni riga di df_clinico corrisponda il cluster corretto
df_clinico$Cluster_Ibrido_K2 <- cluster_ibridi_assegnati[match(df_clinico$v1, names(cluster_ibridi_assegnati))]
df_clinico$Cluster_Ibrido_K2 <- as.factor(df_clinico$Cluster_Ibrido_K2)

# Creazione della tabella di contingenza (entrambe le colonne risiedono in df_clinico -> STESSA LUNGHEZZA)
tabella_ibrida_cancer <- table(df_clinico$Cluster_Ibrido_K2, as.factor(df_clinico$cancer))

cat("--- Tabella di Contingenza (Cluster Ibridi vs Stato Oncologico) ---\n")
print(tabella_ibrida_cancer)

# Calcolo del Test del Chi-quadrato
test_chi2_ibrido <- chisq.test(tabella_ibrida_cancer)

cat("\n--- Risultato del Test del Chi-quadrato ---\n")
print(test_chi2_ibrido)
cat("==================================================\n")



