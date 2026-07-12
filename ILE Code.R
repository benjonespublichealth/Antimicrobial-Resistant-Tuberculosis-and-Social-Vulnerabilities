#####
# Packages

library(janitor)
library(tidyverse)
library(dplyr)
library(visdat)
library(car)
library(nnet)
library(pROC)
library(ggplot2)
library(broom)
library(gtsummary)
library(gt)
library(naniar)
library(forcats)

#####
# Read in data

df <- read.csv("Moldova Dataset Final.csv", na.strings = (""))
df <- clean_names(df)

#####
# Inclusion/Exclusion Criteria + Create Mortality

df_clean <- df %>%
  mutate(tb_mortality = case_when(treatment_result == 51 ~ 1, is.na(treatment_result) ~ NA_real_, TRUE ~ 0)) %>%
  mutate(tb_mortality = factor(tb_mortality, levels = c(0, 1))) %>% 
  filter(!is.na(treatment_result)) %>% 
  filter(!is.na(dst_date))

#####
# Create Drug Resistance Type Variable

df_clean$resistance_num <- rowSums(df_clean[, c("isoniazid", 
                                                "rifampicin", 
                                                "ethambutol",
                                                "streptomycin", 
                                                "pyrazinamide", 
                                                "ethionamide", 
                                                "kanamycin", 
                                                "ofloxacin", 
                                                "cycloserine", 
                                                "paminosalicylic")] == T)
df_clean$mdr <- df_clean$rifampicin & df_clean$isoniazid
df_clean$drug_resistance_type <- with(df_clean, ifelse(resistance_num == 0, "Susceptible", ifelse(mdr == TRUE, "MDR", ifelse(resistance_num == 1, "Mono", "Poly"))))

#####
# Look at Missing Data

vis_miss(df_clean, warn_large_data = FALSE)
gg_miss_var(df_clean, show_pct = TRUE)
# Largest Missing Variable Rate is <2.5% (Homeless)

table(is.na(df_clean$homeless), df_clean$tb_mortality)
# Data with N/As: 5.97%
# Data without N/As: 4.94%
# Similar rates, missingness not dependent on mortality

table(is.na(df_clean$homeless), df_clean$drug_resistance)
# Data with N/As: 22.9%
# Data without N/As: 25.9%
# Similar rates, missingness not dependent on drug resistance

df_complete <- na.omit(df_clean)
missing_data <- ((nrow(df_clean) - nrow(df_complete))/6727)*100
print(missing_data)
# Participants with missing data is 4.7% (<5%) so complete case analysis should introduce minimal bias


#####
# Factor Variables

df_complete$homeless <- factor(df_complete$homeless)
df_complete$homeless <- fct_relevel(df_complete$homeless, "N")

df_complete$gender <- factor(df_complete$gender)
df_complete$gender <- fct_relevel(df_complete$gender, "M")

df_complete$citizenship <- factor(df_complete$citizenship)
df_complete$citizenship <- fct_recode(df_complete$citizenship,
                                     "Moldova" = "Moldova",
                                     "Other" = "other")
df_complete$citizenship <- fct_relevel(df_complete$citizenship, "Moldova")


df_complete$salaried <- factor(df_complete$salaried)
df_complete$salaried <- fct_relevel(df_complete$salaried, "Y")

df_complete$ever_in_prison <- factor(df_complete$ever_in_prison)
df_complete$ever_in_prison <- fct_relevel(df_complete$ever_in_prison, "N")


df_complete$drug_resistance_type <- factor(df_complete$drug_resistance_type)
df_complete$drug_resistance_type <- fct_relevel(as.factor(df_complete$drug_resistance_type), "None")

df_complete$urban_or_rural <- factor(df_complete$urban_or_rural)
levels(df_complete$urban_or_rural) = c("Rural", "Urban")
levels(df_complete$urban_or_rural)

df_complete$occupation <- factor(df_complete$occupation)
df_complete$occupation <- fct_recode(df_complete$occupation,
                                     "Disabled" = "disabled",
                                     "Pensioner" = "pensioner",
                                     "Student" = "student",
                                     "Unemployed" = "unemployed",
                                     "Worker" = "worker")
df_complete$occupation <- fct_relevel(df_complete$occupation, "Worker")

df_complete$education_level <- factor(df_complete$education_level)
df_complete$education_level <- fct_recode(df_complete$education_level,
                                          "Primary School" = "1",
                                          "Secondary School" = "2",
                                          "Specialized Secondary Education" = "3",
                                          "Higher Education" = "4",
                                          "No Education" = "5")
df_complete$education_level <- fct_relevel(df_complete$education_level, "Secondary School")

df_complete$tb_mortality <- fct_recode(df_complete$tb_mortality,
                                       "All Other Treatment Outcomes" = "0",
                                       "Death From TB Progression" = "1")

levels(df_complete$tb_mortality)

#####
# Association Between Drug Resistance & TB Mortality

resistance_mortality_model <- glm(tb_mortality ~ drug_resistance_type +
                                    urban_or_rural +
                                    homeless +
                                    gender +
                                    age_at_diagnosis +
                                    citizenship +
                                    salaried +
                                    ever_in_prison +
                                    education_level +
                                    occupation,
                                  data = df_complete,
                                  family = binomial)
exp(cbind(OR = coef(resistance_mortality_model), confint(resistance_mortality_model)))

# ***** #
vif(resistance_mortality_model)
plot(fitted(resistance_mortality_model))
resist_roc <- roc(df_complete$tb_mortality, fitted(resistance_mortality_model))
plot(resist_roc, legacy.axes = T)
auc(resist_roc)
# ***** #

#####
# Association Between Social Factors & Drug Resistance

social_factors_model <- multinom(drug_resistance_type ~ urban_or_rural +
                                   homeless +
                                   gender +
                                   age_at_diagnosis +
                                   citizenship +
                                   salaried +
                                   ever_in_prison +
                                   education_level +
                                   occupation,
                                 data = df_complete)
exp(coef(social_factors_model))
exp(confint(social_factors_model))

#####
# Effect Modification
model_interaction <- glm(tb_mortality ~ drug_resistance_type * 
                           (urban_or_rural +
                           homeless +
                           gender +
                           age_at_diagnosis +
                           citizenship +
                           salaried +
                           ever_in_prison +
                           education_level +
                           occupation),
                         family = binomial,
                         data = df_complete)
step(model_interaction, direction = "both")
summary(model_interaction)

anova(resistance_mortality_model, model_interaction, test = "Chisq")

#####
# Visuals

tidy_model <- tidy(resistance_mortality_model, conf.int = TRUE, exponentiate = TRUE)

forest_data <- tidy_model %>%
  filter(term %in% c("drug_resistance_typeMono", 
                     "drug_resistance_typePoly", 
                     "drug_resistance_typeMDR")) %>%
  mutate(term = dplyr::recode(term, 
                              "drug_resistance_typeMono" = "Mono-resistance",
                              "drug_resistance_typePoly" = "Poly-resistance",
                              "drug_resistance_typeMDR" = "MDR-TB"))

ggplot(forest_data, aes(x = term, y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  labs(title = "Association Between Drug Resistance and TB Mortality",
       subtitle = "Reference Group: Drug-Sensitive (None)",
       x = "",
       y = "Adjusted Odds Ratio (95% CI)") +
  theme_minimal()

#

table1 <- df_complete %>%
  select(drug_resistance_type,
         age_at_diagnosis,
         gender,
         urban_or_rural,
         salaried,
         citizenship,
         homeless,
         ever_in_prison,
         education_level,
         occupation,
         tb_mortality) %>%
  tbl_summary(
    by = drug_resistance_type,
    type = list(
      age_at_diagnosis ~ "continuous"),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})", 
      all_categorical() ~ "{n} ({p}%)"),
    label = list(
      age_at_diagnosis ~ "Age (Years)",
      gender ~ "Gender",
      urban_or_rural ~ "Type of Residence",
      salaried ~ "Salaried",
      citizenship ~"Citizenship",
      homeless ~ "Ever Experienced Homelessness",
      ever_in_prison ~ "Ever Incarcerated",
      education_level ~ "Education Level",
      occupation ~ "Occupation",
      tb_mortality ~ "TB Mortality"),
    missing = "no") %>%
  add_overall() %>%
  bold_labels()

table1

as_gt(table1) %>%
  gtsave("Table1_Descriptive.docx")


table(df_complete$drug_resistance_type, df_complete$occupation, df_complete$tb_mortality)

car::vif(model_interaction)

model_interaction_revised <- update(model_interaction, . ~ . - drug_resistance_type:occupation)
summary(model_interaction_revised)
anova(model_interaction_revised, model_interaction, test = "Chisq")

