library(here)
library(readxl)
library(dplyr)
library(survival)
library(survminer)
library(parameters)
library(GGally)
library(ggsurvfit)
#library(cmprsk)
library(tidycmprsk)

# Charger les données depuis le fichier XLSX
data <- read_excel(here("Copy_of_Oncology_dataset_for_R.xlsx"))
str(data)
data_structure <- data.frame(
  Column = names(data),
  Type = sapply(data, class),
  Example_Values = sapply(data, function(x) paste(head(x, 3), collapse = ", "))
)
data_structure

# Définir les types des colonnes
data$Cell <- as.factor(data$Cell)          # Type de cellule (catégoriel)
data$Therapy <- as.factor(data$Therapy)    # Thérapie (catégoriel)
data$Prior <- as.factor(data$Prior)        # Thérapie préalable (binaire : 0/10)
data$Age <- as.numeric(data$Age)           # Âge (numérique)
data$DiagTime <- as.numeric(data$DiagTime) # Temps depuis le diagnostic (numérique)
data$Kps <- as.numeric(data$Kps)           # Statut de performance (numérique)
data$Censor <- as.factor(data$Censor)     # Indicateur de censure (binaire)
data$Treatment <- as.factor((data$Treatment))  # binaire (Standard ou Test)

str(data)

# 1. What was the maximum survival time for the cell type adeno?
#---------------------------------------------------------------

# Filtrer les données pour le type de cellule "adeno"
adeno_data <- subset(data, Cell == "adeno")
#adeno_data <- data %>% filter(Cell == "adeno")

# Trouver le temps de survie maximum
max_survival_time <- max(adeno_data$SurvTime, na.rm = TRUE)
#max_survival_time <- max(data$SurvTime[data$Cell == "adeno"], na.rm = TRUE)

# Afficher le résultat
print(max_survival_time)


# 2. What is the average age of subjects in this study?
#------------------------------------------------------

# Calculer l'âge moyen
average_age <- mean(data$Age, na.rm = TRUE)

# Afficher le résultat
print(average_age)


# 3. Which cell type appeared the most during this study?
#--------------------------------------------------------

# Trouver le type de cellule le plus fréquent
most_common_cell <- names(which.max(table(data$Cell)))

# Afficher le résultat
print(most_common_cell)


# 4. Calculate descriptive statistics for all numeric variables within this dataset?
#-----------------------------------------------------------------------------------

# Résumé descriptif pour toutes les colonnes numériques
#numeric_summary <- summary(data)
summary(data)
# Afficher le résultat
#print(numeric_summary)


# 5. Perform a survival analysis to assess the survival time (variable SurvTime)?
# based on the cancerous cells (var Cell)? Consider applying survival
# functions/kaplan meier quartiles/cumulative incidence function/cox regression etc.
#-------------------------------------------------------------------------------

# Vérifier et préparer les données
data$Censor <- as.numeric(data$Censor)  # Si nécessaire, convertir en numérique
data$Cell <- as.factor(data$Cell)       # Convertir Cell en facteur si ce n'est pas déjà fait

# Créer un objet de survie
surv_obj <- Surv(time = data$SurvTime, event = data$Censor)

# 1. Courbes de survie de Kaplan-Meier pour chaque type de cellule :
# comparer visuellement la survie entre les différents types de cellules cancéreuses.
km_fit <- survfit(surv_obj ~ Cell, data = data)

# Résumé des résultats
summary(km_fit)

# Tracer les courbes de survie
plot(km_fit, col = 1:length(levels(data$Cell)),
     xlab = "Temps en jours", ylab = "Probabilité de survie",
     main = "Courbes de survie par type de cellule")
legend("bottomleft", legend = levels(data$Cell), col = 1:length(levels(data$Cell)), lty = 1)

# Tracer les courbes de Kaplan-Meier avec ggsurplot
ggsurvplot(
  km_fit,
  data = data,
  pval = TRUE,                # Afficher la p-valeur du test Log-Rank
  conf.int = TRUE,            # Intervalle de confiance
  risk.table = TRUE,          # Ajouter une table des risques
  legend.title = "Cell Types",
  xlab = "Time (days)",
  ylab = "Survival Probability",
  ggtheme = theme_minimal()
)


# Test du log-rank
surv_diff <- survdiff(surv_obj ~ Cell, data = data)
print(surv_diff)

# Calculer la p-value
p_value <- 1 - pchisq(surv_diff$chisq, df = length(levels(data$Cell)) - 1)
cat("P-value du test du log-rank :", p_value, "\n", "Une p-value inférieure à 0,05 indique une différence significative entre les groupes.")

# 2. Régression de Cox pour évaluer l'effet du type de cellule
cox_fit <- coxph(surv_obj ~ Cell, data = data)

# Résumé du modèle
summary(cox_fit)
parameters(cox_fit, exponentiate = TRUE)
cat("Le résumé fournira les coefficients, les erreurs standards, les valeurs z, les p-values, et les hazard ratios (HR).")

# Représentation graphique du modèle de Cox
ggcoef_model(cox_fit, exponentiate = TRUE)

# Ajout des covariables
# Cela te permettra de voir si l'association entre 'Cell' et la survie persiste après ajustement pour les autres facteurs.
cox_fit_cov <- coxph(surv_obj ~ Cell + Age + Prior + Kps + Treatment + DiagTime, data = data)
summary(cox_fit_cov)
parameters(cox_fit_cov, exponentiate = TRUE)
ggcoef_model(cox_fit_cov, exponentiate = TRUE)

# Evaluation des conditions d’application
## Risques Proportionnels
# Test des risques proportionnels
cox_zph <- cox.zph(cox_fit_cov)
print(cox_zph)
# ==> Aucune p-value n’est inférieure à 0,05, ainsi rien ne permet de conclure que l’hypothèse des risques proportionnels n’est pas respectée pour aucune des covariables incluses dans le modèle.

# Visualiser les résidus de Schoenfeld
par(mfrow=c(2,3)) # division de la fenêtre graphique en 2 lignes et 3 colonnes
plot(cox_zph)




# Ajuster une courbe Kaplan-Meier
fit <- survfit(Surv(SurvTime, Censor) ~ 1, data = data)

# Calcul des quartiles
quartiles <- quantile(fit, probs = c(0.25, 0.50, 0.75))$quantile

# Annoter manuellement les quartiles sur le graphique
ggsurv$plot <- ggsurv$plot +
  geom_hline(yintercept = c(0.75, 0.50, 0.25), linetype = "dashed", color = "red") +
  geom_vline(xintercept = quartiles, linetype = "dotted", color = "blue") +
  annotate("text", x = quartiles, y = c(0.75, 0.50, 0.25), 
           label = paste0("Q", 1:3), color = "blue", hjust = -0.1)

# Afficher le graphique mis à jour
print(ggsurv)


# 3. Fonction de densité cumulative
cum_inc_fit <- survfit(Surv(SurvTime, Censor) ~ Cell, data = data)
summary(cum_inc_fit)

# Tracer la fonction cumulative
ggsurvplot(
  cum_inc_fit,
  data = data,
  fun = "cumhaz",         # Risque cumulé
  conf.int = TRUE,        # Ajouter l'intervalle de confiance
  risk.table = TRUE,      # Afficher un tableau des risques
  pval = TRUE,            # Ajouter la p-valeur du test Log-Rank
  legend.title = "Cell Types",
  xlab = "Time (days)",
  ylab = "Cumulative Hazard",
  ggtheme = theme_minimal()
)


ggsurvfit(surv_fit) +
  add_confidence_interval() + # Ajouter l'intervalle de confiance
  add_risktable() +           # Ajouter un tableau des risques
  scale_y_continuous(limits = c(0, 1)) + # Fixer les limites de l'axe Y
  labs(
    title = "Courbes de survie selon les types de cellules",
    x = "Temps de survie (jours)",
    y = "Probabilité de survie"
  ) +
  theme_minimal()

# 6. Perform an appropriate multivariable analysis to analyze the effect of
# independent variables age on the hazard ratio between the different cancerous
# cells (var Cell)?

# Ajuster le modèle de régression de Cox avec l'âge comme covariable
cox_model <- coxph(surv_obj ~ Cell + Age, data = data)

# Résumé des résultats
summary(cox_model)
parameters(cox_model, exponentiate = TRUE)

ggsurvplot(
  adjusted_surv_fit,
  data = data,
  conf.int = TRUE,
  ggtheme = theme_minimal(),
  xlab = "Time (days)",
  ylab = "Survival Probability",
  legend.title = "Cancer Cell Type"
)


# Test proportional hazards assumption
ph_test <- cox.zph(cox_model)
print(ph_test)
# Visualize Schoenfeld residuals for each covariate
cat(" greater than 0.05, the proportional hazards assumption is not violated for that variable.")

# Visualiser les résultats avec forest plot
ggforest(cox_model, data = data, main = "Hazard Ratios for Cancerous Cells with Age Adjustment")






