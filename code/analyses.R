# ---
# title: "Analyses"
# output: html_document
# date: "2024-09-11"
# ---
# 
# # Set-Up
# 
# We load the necessary R packages, install additional packages, and define custom functions that will be used throughout the analysis.
# 
## Load Packages
# 

dir.create("../results/all_chat_exports/", showWarnings = FALSE)

#setwd("~/Dropbox (Princeton)/rechat pilots/papers/introducing rechat/replicable_analyses/v11_codeocean/code")

#library(tidyverse)
library(dplyr)
library(lubridate)
library(tidyr)
library(stringr)
library(purrr)
library(readr)

library(haven)
library(stargazer)
library(sandwich)
library(lmtest)
library(webshot)
#library(stm)
library(sentimentr)
library(moments)
library(pdftools)
library(psych)

# 
## Install rechat R package
# 
# Install the rechat package from GitHub repository
# library(devtools)
# install_github("willschulz/rechat")
library(rechat)
# 
# 
## Create Additional Functions

# Function credit: pewmethods/pewmethods (couldn't install pewmethods package on codeocean)
fct_case_when <- function(...)
{
  local({library(rlang)})
  default_env <- caller_env()
  arguments <- list2(...)
  arguments <- Filter(function(elt) !is.null(elt), arguments)
  arg_len <- length(arguments)
  output_levels <- purrr::map(arguments, function(a) {
    out <- f_rhs(a) %>% eval_tidy(env = default_env)
    return(levels(out) %||% out)
  }) %>%
    squash_chr()
  for (i in 1:arg_len) {
    f_rhs(arguments[[i]]) <- as.character(f_rhs(arguments[[i]]) %>%
                                            eval_tidy(env = default_env))
  }
  cw <- do.call(case_when, arguments)
  cw <- factor(cw, levels = unique(output_levels))
  return(cw)
}

# Custom function to resize stargazer tables for LaTeX output
resizebox.stargazer = function(..., tab.width = "!", tab.height = "!"
                               ){
  input_list <- as.list(substitute(list(...)))

  # Load stringr package for string manipulation
  require(stringr) 

  #Extract the code returned from stargazer()
  res = capture.output(
    stargazer::stargazer(...)
    )

  #Render the arguments:
  tab.width = tab.width
  tab.height = tab.height

  #Attach "}" between \end{tabular} and \end{table}
  res = 
    prepend(res, "}", before = length(res))

  #Input \resizebox before \begin{tabular}
  res = 
    c(res[1:str_which(res, "^\\\\begin\\{tabular\\}.*")-1],
      paste0("\\resizebox{",tab.width,"}{",tab.height,"}{%"),
      #paste0("\\scalebox{",tab.width,"}{",tab.height,"}{%"),
      res[str_which(res, "^\\\\begin\\{tabular\\}.*"):length(res)]
      )

  #Produce the whole strings
  cat(res, sep = "\n", file = input_list$out)
}

# 
# 
# Custom function to print significance levels with p-values
print_sig_level <- function(x, max_places = 50){
  levels <- c(.1, .05, 1/(10^seq(2:max_places)))
  levels_greater <- levels[which(levels>x)]
  level_to_print <- levels_greater[which.min(levels_greater)]
  paste0("p<",level_to_print)
}
# 
# 
# Custom function to make robust standard errors
make_robust_se <- function(mod){
  cov1         <- vcovCL(mod, cluster = ~room_id)
  robust_se    <- sqrt(diag(cov1))
  return(robust_se)
}
# 
# 
# # Prompt Selection Analyses
# 
# This section re-analyzes data from Druckman et al (2022). The analysis explores polarization in trait ratings between in-group and out-group members, to identify traits characterized by high polarization.
# 
# Load the dataset and filter by specific waves
druckman_traits <- read_dta(file = "../data/druckman/Data Misestimating AP September 5 2020.dta") %>%
  filter(wave1pid %in% c(1,2,3,5,6,7)) %>%
  dplyr::select(ends_with("_ingroup") | ends_with("_outgroup"), wave1pid) %>%
  dplyr::select(-starts_with("feelingtherm_"), -starts_with("marry_"), -starts_with("category_"), -starts_with("trust_"), -starts_with("friends_"), -starts_with("neighbors_")) %>%
  drop_na()

# Extract the party identification variable
pid <- druckman_traits %>% pull(wave1pid)

# Calculate the difference between in-group and out-group ratings
in_out <- (druckman_traits %>% dplyr::select(ends_with("_ingroup")) %>% as.matrix) - (druckman_traits %>% dplyr::select(ends_with("_outgroup")) %>% as.matrix)

# Clean up column names by removing the suffix "_ingroup"
colnames(in_out) <- colnames(in_out) %>% str_remove_all("_ingroup")

# Recode so that positive differences reflect more polarization (favoring in-group)
infavor_num <- in_out
infavor_num[,which(colMeans(in_out)<0)] <- -infavor_num[,which(colMeans(in_out)<0)]

# Create a binary indicator for in-group favoritism
infavor_bin <- (infavor_num>0)
# 
# 
# Plot traits' polarization by party identification
means_all <- colMeans(infavor_bin)
means_rep <- colMeans(infavor_bin[which(pid>4),])
means_dem <- colMeans(infavor_bin[which(pid<4),])

pdf("../results/trait_favoritism.pdf", width = 8.5, height = 5)

plot(1:length(means_all), means_all, ylim = c(.4,.8), pch = 16, ylab = "Mean In-Group Favoritism", main = "Mean In-Group Favoritism in Trait Ratings", xlab = "", xaxt = "n")
segments(x0 = 1:length(means_all), y0 = means_rep, y1 = means_dem)
points(1:length(means_all), means_rep, col = "red", pch = 16)
points(1:length(means_all), means_dem, col = "blue", pch = 16)
axis(1, at = 1:length(means_all), labels = names(means_all), cex.axis = 1)
abline(h = .5, col = "black", lty=3)

legend("topleft", pch = 16, col = c("black", "red", "blue"), legend = c("All Partisans", "Republicans", "Democrats"), bty = "n")

dev.off()
# 
# 
# # Study 1
# 
# This section conducts all analyses focused on Study 1.
# 
## Prep Data
# 
### Load Data
# 
# Load raw survey data for Study 1
s1_raw_survey_data <- readRDS("../data/study_1/raw_survey_data.rds")
# Recode the treatment variable (t = 1-t) so that neutrality is the treatment, per reviewer request
s1_raw_survey_data <- s1_raw_survey_data %>% mutate(t = 1-t)
# 
# 
# Load and parse chat data for Study 1
chat_data <- parseChat("../data/study_1/47cb1084-4cf5-4697-be97-f6447c9f186a.csv")
chat_data_republicans <- parseChat("../data/study_1/3583fe12-0f20-44b6-a08d-878ce3842b58.csv") #(a handful of republicans were recruited despite our filtering criteria; they are excluded from analysis but their data is included for completeness)
# 
# 
# Drop test chats from chat_data
chat_data <- chat_data[-(1)]
chat_data <- c(chat_data, chat_data_republicans)
# 
# 
### Clean Data
# 
# Drop empty chats (may result from technical issues)
empty <- c()
for(i in 1:length(chat_data)){
  empty[i] <- nrow(chat_data[[i]]$messages)==0
}

chat_data <- chat_data[-which(empty)]

# Drop "chats" where partner did not send any messages (may result from non-compliance/technical issues)
single <- c()
for(i in 1:length(chat_data)){
  single[i] <- chat_data[[i]]$messages %>% pull(participantCode) %>% {length(unique(.))<2}
}

chat_data <- chat_data[-which(single)]
# 
# 
# Define a function to recode news consumption variables
news_recode <- function(x){
  return(case_when(x == "Never" ~ 0,
          x == "Less often" ~ 1,
          x == "1-2 days a week" ~ 2,
          x == "3-6 days a week" ~ 3,
          x == "About once a day" ~ 4,
          x == "Several times a day" ~ 5))
}
# 
# 
# Clean and rename/recode variables in the survey data

s1_survey_data <- s1_raw_survey_data %>%
  rename(ego_code = number_entry, #participant confirmation code
         mot_get_along = motivations_1, #motivation: to get along
         mot_say_exactly = motivations_2, #motivation: to say exactly what I think
         mot_finish_quick = motivations_3, # motivation: to finish quickly
         mot_good_rating = motivations_4) %>% #motivation: to get a good rating
  mutate(
    chat_consent = case_when((chat_consent=="No, I do not want to participate") ~ FALSE, #recode chat consent as binary
                             (chat_consent=="Yes, I want to participate") ~ TRUE),
    affpol_pre = case_when(partisanship == "Democrat" ~ thermo_pre_4 - thermo_pre_7, #pre-chat affective polarization
                           partisanship == "Republican" ~ thermo_pre_7 - thermo_pre_4),
    affpol_post = case_when(partisanship == "Democrat" ~ thermo_post_4 - thermo_post_7, #post-chat affective polarization
                            partisanship == "Republican" ~ thermo_post_7 - thermo_post_4),
    affpol_change = affpol_post - affpol_pre, #change in affective polarization
    demthermo_change = thermo_post_4 - thermo_pre_4, #change in thermometer ratings for Democrats
    repthermo_change = thermo_post_7 - thermo_pre_7, #change in thermometer ratings for Republicans
    PID_5 = case_when( #recode party identification to 5-point scale
      PID_dem_strength == "Strong" ~ -2,
      PID_dem_strength == "Not very strong" ~ -1,
      PID_leaners == "Closer to Democratic Party" ~ 0,
      PID_leaners == "Closer to Republican Party" ~ 0,
      PID_rep_strength == "Not very strong" ~ 1,
      PID_rep_strength == "Strong" ~ 2
    ),
    PID_6 = case_when( #recode party identification to 6-point scale
      PID_dem_strength == "Strong" ~ 1,
      PID_dem_strength == "Not very strong" ~ 2,
      PID_leaners == "Closer to Democratic Party" ~ 3,
      PID_leaners == "Closer to Republican Party" ~ 4,
      PID_rep_strength == "Not very strong" ~ 5,
      PID_rep_strength == "Strong" ~ 6
    ),
    
    id_important = case_when( #recode identity importance
      huddy1 =="Not important at all" ~ 0,
      huddy1 =="Not very important" ~ 1,
      huddy1 =="Very important" ~ 2,
      huddy1 =="Extremely important" ~ 3
    ),
    id_describe = case_when( #recode identity descriptiveness
      huddy2 =="Not at all" ~ 0,
      huddy2 =="Not very well" ~ 1,
      huddy2 =="Very well" ~ 2,
      huddy2 =="Extremely well" ~ 3
    ),
    id_wethey = case_when( #recode identity we-they
      huddy3 =="Never" ~ 0,
      huddy3 =="Rarely" ~ 1,
      huddy3 =="Some of the time" ~ 2,
      huddy3 =="Most of the time" ~ 3,
      huddy3 =="All of the time" ~ 4
    ),
    
    # recode media consumption variables
    tv = news_recode(news_media_1),
    newspapers = news_recode(news_media_2),
    radio = news_recode(news_media_3),
    internet_sm = news_recode(news_media_4),
    discussions = news_recode(news_media_5),
    podcasts = news_recode(news_media_6),
    
    # recode social media usage variables
    twitter = !is.na(social_media_1),
    facebook = !is.na(social_media_2),
    instagram = !is.na(social_media_3),
    reddit = !is.na(social_media_9),
    tiktok = !is.na(social_media_10),
    youtube = !is.na(social_media_11),
    
    #primary_sm = primary_sm, #keep primary_sm as coded
    
    age = (2021 - birthyr), #calculate age
    
    male = (gender=="Male"), #binary indicator for male gender
    
    #education = education, #keep education as coded
    college = as.numeric(education)>=4, #binary indicator for college education
    
    expressor = !(primary_sm %in% c("I don't use any of these services to express my opinions", NA)), #binary indicator for whether one expresses opinions about politics or current events on social media
    
    sm_postpol_num = as.numeric(sm_postpol), #numeric version of self-reported frequency of political posting on social media
    
    ideo_7 = as.numeric(ideology), #numeric version of 7-point ideology scale
    PID_strength = (3-PID_6)
  ) %>% filter(ideo_7<=4) #filter out conservative participants (screening at recruitment stage let in some conservatives)

# recode guesses of chat partner's ideology
s1_survey_data$guess_partner_ideo <- as.numeric(s1_survey_data$guess_partner_ideo)
s1_survey_data$guess_partner_ideo[which(s1_survey_data$guess_partner_ideo==8)] <- NA

# recode enjoyment of chat
s1_survey_data$enjoyment <- as.numeric(s1_survey_data$enjoyment)
s1_survey_data$enjoyment[which(s1_survey_data$enjoyment==5)] <- NA
s1_survey_data$enjoyment <- max(s1_survey_data$enjoyment, na.rm = T)-s1_survey_data$enjoyment

# clean up posting variable
s1_survey_data$sm_postpol_num[which(s1_survey_data$sm_postpol_num == 6)] <- NA
s1_survey_data$sm_postpol_num_nazero <- s1_survey_data$sm_postpol_num
s1_survey_data$sm_postpol_num_nazero[which(is.na(s1_survey_data$sm_postpol_num_nazero))] <- 0

#notes:
#thermo_pre_4 is Democrats
#thermo_pre_7 is Republicans

s1_survey_data <- s1_survey_data %>% filter(partisanship!="Republican") #filter out Republicans (screening at recruitment stage let in some Republicans)

# 
# 
# Make a variable for "ideological extremity" as requested by reviewer
s1_survey_data <- s1_survey_data %>% mutate(ideo_ext = 4-ideo_7)

# 
# 
### PCA Scales
# 
#Recode agree-disagree scales
ag_disag_recode <- function(x){
  return(case_when(x == "Strongly disagree" ~ 0,
          x == "Disagree" ~ 1,
          x == "Neither agree nor disagree" ~ 2,
          x == "Agree" ~ 3,
          x == "Strongly agree" ~ 4))
}

extro_cols <- s1_survey_data %>% dplyr::select(ends_with("_extroversion"))
extro_cols <- apply(extro_cols, MARGIN = 2, ag_disag_recode)
# 
# 
# PCA to create extroversion scale
extro_scale <- psych::principal(extro_cols, nfactors = 1, rotate = "varimax", missing=TRUE, impute = "mean")$scores[,1]

# Flip scale automatically if necessary
if (cor(extro_scale, extro_cols[,5])<0) {
  extro_scale <- (-extro_scale)
}

s1_survey_data$extro_scale <- extro_scale # Add extroversion scale to dataframe
# 
# 
# PCA to create identity scale
identity_scale <- psych::principal(s1_survey_data %>% dplyr::select(starts_with("id_")), 
                              nfactors = 1, rotate = "varimax", missing=TRUE, impute = "mean")$scores[,1]
if (cor(identity_scale, s1_survey_data$id_important)<0) {
  identity_scale <- (-identity_scale)
}
s1_survey_data$identity_scale <- identity_scale
# 
# 
# PCA to create social media scale
sm_scale <- psych::principal(s1_survey_data %>% dplyr::select(twitter, facebook, instagram, reddit, tiktok, youtube),
                              nfactors = 1, rotate = "varimax", missing=TRUE, impute = "mean")$scores[,1]
s1_survey_data$sm_scale <- sm_scale
# 
# 
# PCA to create media scale
media_scale <- psych::principal(s1_survey_data %>% dplyr::select(tv, newspapers, radio, internet_sm, discussions, podcasts),
                              nfactors = 1, rotate = "varimax", missing=TRUE, impute = "mean")$scores[,1]
s1_survey_data$media_scale <- media_scale
# 
# 
# PCA to create self-monitoring scale
selfmon_scale <- psych::principal(s1_survey_data %>% transmute(as.numeric(selfMon1), as.numeric(selfMon2), as.numeric(selfMon3)),
                              nfactors = 1, rotate = "varimax", missing=TRUE, impute = "mean")$scores[,1]

if(cor(selfmon_scale, as.numeric(s1_survey_data$selfMon2))>0){selfmon_scale <- (-selfmon_scale)}

s1_survey_data$selfmon_scale <- selfmon_scale
# 
# 
# Combine extroversion and self-monitoring into a single scale
extroselfmon_scale <- psych::principal(cbind(extro_cols, s1_survey_data %>% transmute(as.numeric(selfMon1), as.numeric(selfMon2), as.numeric(selfMon3))),
                              nfactors = 1, rotate = "varimax", missing=TRUE, impute = "mean")$scores[,1]

if (cor(extroselfmon_scale, extro_cols[,5])<0) {
  extroselfmon_scale <- (-extroselfmon_scale)
}
s1_survey_data$extroselfmon_scale <- extroselfmon_scale
# 
# 
## Analyses
# 
### Summarize Sample Size (Recruitment/Attrition)
# 
# Check N at each stage of study
s1_survey_data %>% filter(consent_response=="Yes, I want to participate") %>% nrow #483
s1_survey_data %>% filter(consent_response=="Yes, I want to participate") %>% filter(chat_consent) %>% nrow #349

s1_survey_data %>% filter(consent_response=="Yes, I want to participate") %>% pull(chat_consent) %>% mean

# Export percentage consenting to txt
s1_consent_percentage <- s1_survey_data %>% filter(consent_response=="Yes, I want to participate") %>% pull(chat_consent) %>% mean %>% {.*100} %>% round

s1_consent_percentage %>% paste0(., "\\%") %>% writeLines(., "../results/s1_consent_percentage.txt")

# 
# 
### DV: Self-Selection into Chat Participation
# 
# Estimate logistic regression models for self-selection into chat participation

log_ex_base <- glm(chat_consent ~ extro_scale, data = s1_survey_data, family = "binomial")

log_sm_base <- glm(chat_consent ~ selfmon_scale, data = s1_survey_data, family = "binomial")

log_ex <- glm(chat_consent ~ extro_scale + selfmon_scale + as.numeric(polInt) + PID_strength + affpol_pre + ideo_7 + identity_scale + media_scale, data = s1_survey_data, family = "binomial")

log_ex_ks <- glm(chat_consent ~ extro_scale + selfmon_scale + as.numeric(polInt) + PID_strength + affpol_pre + ideo_7 + identity_scale + media_scale + male + age + college, data = s1_survey_data, family = "binomial")

log_ex_sm <- glm(chat_consent ~ extro_scale + as.numeric(polInt) + PID_strength + affpol_pre + ideo_7 + identity_scale + media_scale + selfmon_scale + male + age + college + expressor
                             , data = s1_survey_data, family = "binomial")

# 
# 
# 
# Export regression tables
stargazer(log_ex_base, log_sm_base, log_ex, log_ex_ks, log_ex_sm, type = "latex", label = "tab:s1_selection_into_chat", title = "Self-Selection into Chat Participation (Study 1)", out = "../results/s1_selection_into_chat.tex", table.placement = "H",
          single.row = T)

resizebox.stargazer(log_ex_base, log_sm_base, log_ex, log_ex_ks, log_ex_sm, 
                    type = "latex", label = "tab:s1_selection_into_chat", title = "Self-Selection into Chat Participation (Study 1)",
                    out = "../results/s1_selection_into_chat_resized.tex",
                    table.placement = "H",
                    covariate.labels = c("Extroversion", "Self-Monitoring", "Political Interest", "Partisanship Strength", "Affective Polarization", "7-Point Ideology", "Ideological Identity Strength", "Media Consumption Scale", "Male", "Age", "College", "Social Media Expressor"),
                    dep.var.labels = "Chat Consent",
                    tab.width = "\\textwidth", omit.stat = c("f"), single.row = TRUE)

# 
# 
#Export in-text p-values
summary(log_ex_base)$coefficients[2,4] %>% round(digits = 3) %>% str_replace(., "0.", "($p=.") %>% paste0(.,"$)") %>% writeLines(., "../results/s1_extro_p.txt")

summary(log_sm_base)$coefficients[2,4] %>% round(digits = 3) %>% str_replace(., "0.", "($p=.") %>% paste0(.,"$)") %>% writeLines(., "../results/s1_sm_p.txt")

summary(log_ex)$coefficients[4,4] %>% round(digits = 3) %>% str_replace(., "0.", "($p=.") %>% paste0(.,"$)") %>% writeLines(., "../results/s1_polint_p.txt")

summary(log_ex_ks)$coefficients[11,4] %>% round(digits = 3) %>% str_replace(., "0.", "($p=.") %>% paste0(.,"$)") %>% writeLines(., "../results/s1_college_p.txt")

summary(log_ex_sm)$coefficients[12,4] %>% round(digits = 3) %>% str_replace(., "0.", "($p=.") %>% paste0(.,"$)") %>% writeLines(., "../results/s1_expressor_p.txt")

# 
# 
#### Was extroversion correlated with something else?
# 
# 
# Investigate whether extroversion was correlated with other predictors

extro_corr_otherpredictors <- lm(data = s1_survey_data, extro_scale ~ as.numeric(polInt) + college + expressor)

stargazer(extro_corr_otherpredictors, type = "latex", label = "tab:s1_extro_corr_otherpredictors", omit.stat = c("f"), title = "Predictors of Extroversion", out = "../results/s1_extro_corr_otherpredictors.tex", table.placement = "H", single.row = T, dep.var.labels = "Extroversion",
          covariate.labels = c("Political Interest", "College", "Social Media Expressor"))

# 
# 
### Loquaciousness Analyses
# 
# Function to count words in a message
word_counter <- function(message){
  return(str_count(message, "\\w+"))
}

chat_data <- featurizeChat(chat_data, featurization_function = nchar) #featurize chat data: character count
chat_data <- featurizeChat(chat_data, featurization_function = word_counter) #featurize chat data: word count

# Summarize features and add to survey data
s1_survey_data <- summarizeChat(s1_survey_data, chat_data, chat_feature_name = "nchar", summary_function = sum, na.rm = T, summary_feature_name = "char_count")
s1_survey_data <- summarizeChat(s1_survey_data, chat_data, chat_feature_name = "word_counter", summary_function = sum, na.rm = T, summary_feature_name = "word_count")
# 
# 
# Function to calculate other loquaciousness metrics
more_thoroughness <- function(chat){
  #link from participant codd to color code
  out <- chat$messages %>%
    group_by(participantCode) %>%
    summarize(message_count = length(message),
              unique_word_count = length(unique(unlist(strsplit(message, " ")))),
              median_word_length = median(nchar(unlist(strsplit(message, " ")))),
              mean_word_length = mean(nchar(unlist(strsplit(message, " "))))
              )
  out <- chat$participants %>% dplyr::select(receiptCode, colorCode) %>% left_join(., out, by = c("colorCode" = "participantCode"))
  out <- out %>% mutate(room_id = chat$room_id, started_at = chat$started_at + 12*60*60) #manually adjust for 24-hour format
  return(out)
}

# Function to apply the above function to all chats
more_thoroughnesses <- function(chats, simplify = TRUE){
  out <- lapply(chats, more_thoroughness)
  if (simplify) {out <- data.table::rbindlist(l = out)}
  return(out)
}

# Implement function to calculate other loquaciousness metrics and add to survey data
more_thoroughnesses <- more_thoroughnesses(chat_data)
s1_survey_data_nona <- left_join(s1_survey_data %>% filter(!is.na(ego_code)), more_thoroughnesses %>% filter(!is.na(receiptCode)) %>% dplyr::select(-c("colorCode", "started_at")), by = c("ego_code" = "receiptCode"))
# 
# 
# Estimate linear regression models for loquaciousness

s1_lm_mcount <- lm(message_count ~ t + extro_scale + selfmon_scale + as.numeric(polInt) + PID_strength + affpol_pre + ideo_ext + identity_scale + media_scale + male + age + college
                       + expressor
                       ,
                       data = s1_survey_data_nona %>% filter(!is.na(room_id)))

s1_lm_ccount <- lm(char_count ~ t + extro_scale + selfmon_scale + as.numeric(polInt) + PID_strength + affpol_pre + ideo_ext + identity_scale + media_scale + male + age + college
                    + expressor
                    ,
                    data = s1_survey_data_nona %>% filter(!is.na(room_id)))

s1_lm_wcount <- lm(word_count ~ t + extro_scale + selfmon_scale + as.numeric(polInt) + PID_strength + affpol_pre + ideo_ext + identity_scale + media_scale + male + age + college
                    + expressor
                    ,
                    data = s1_survey_data_nona %>% filter(!is.na(room_id)))

s1_lm_uwcount <- lm(unique_word_count ~ t + extro_scale + selfmon_scale + as.numeric(polInt) + PID_strength + affpol_pre + ideo_ext + identity_scale + media_scale + male + age + college
                           + expressor
                           , data = s1_survey_data_nona %>% filter(!is.na(room_id)))

# Export regression tables
resizebox.stargazer(s1_lm_ccount, s1_lm_mcount, s1_lm_wcount, s1_lm_uwcount, type = "latex",
                     se = list(make_robust_se(s1_lm_ccount),make_robust_se(s1_lm_mcount),make_robust_se(s1_lm_wcount),make_robust_se(s1_lm_uwcount)), title = "Loquaciousness (Study 1)", label = "tab:s1_loquaciousness", out = "../results/s1_loquaciousness_resized.tex", table.placement = "H", omit.stat = c("f"), single.row = TRUE,
                    covariate.labels = c("Treatment\\\\(Neutral Accountability)", "Extroversion", "Self-Monitoring", "Political Interest", "Partisanship Strength", "Affective Polarization", "Ideological Extremity", "Ideological Identity Strength", "Media Consumption Scale", "Male", "Age", "College", "Social Media Expressor"),
                    dep.var.labels = c("Character Count", "Message Count", "Word Count", "Unique Word Count"),
                    tab.width = "\\textwidth"
                    )


# 
# 
# ### Analysis of Free Responses
# # This has been commented out because STM proved difficult to install on Code Ocean, and the analyses that depend on STM are purely supplementary.
# # 
# # Here, we analyze free responses offered by participants who refused to participate in a chat.  We use a Structural Topic Model to explore the content of these free responses.  These results are not reported in the paper, but are provided here as a supplement.
# # 
# reasons <- s1_survey_data$why_not[which(!is.na(s1_survey_data$why_not))]#free responses from participants who refused to chat

# reasons_processed <- textProcessor(reasons, verbose = F)#process text
# stm_0 <- stm(documents = reasons_processed$documents, vocab = reasons_processed$vocab, K = 0, init.type = "Spectral", verbose = F)#fit stm, using spectral initialization to choose topic number

# stm_0
# labelTopics(stm_0)

# stm_2 <- stm(documents = reasons_processed$documents, vocab = reasons_processed$vocab, K = 2, init.type = "Spectral", verbose = F) #fit stm with 2 topics
# stm_3 <- stm(documents = reasons_processed$documents, vocab = reasons_processed$vocab, K = 3, init.type = "Spectral", verbose = F) #fit stm with 3 topics
# stm_4 <- stm(documents = reasons_processed$documents, vocab = reasons_processed$vocab, K = 4, init.type = "Spectral", verbose = F) #fit stm with 4 topics
# stm_5 <- stm(documents = reasons_processed$documents, vocab = reasons_processed$vocab, K = 5, init.type = "Spectral", verbose = F) #fit stm with 5 topics
# stm_6 <- stm(documents = reasons_processed$documents, vocab = reasons_processed$vocab, K = 6, init.type = "Spectral", verbose = F) #fit stm with 6 topics

# #View topic content

# labelTopics(stm_2)
# labelTopics(stm_3) #this seems like the most interpretable one
# labelTopics(stm_4)
# labelTopics(stm_5)
# labelTopics(stm_6)

# findThoughts(stm_2, texts = reasons, n = 2, topics = 1)$docs[[1]]
# findThoughts(stm_3, texts = reasons, n = 2, topics = 3)$docs[[1]]
# findThoughts(stm_4, texts = reasons, n = 2, topics = 4)$docs[[1]]

# 
# 
## Export Chats
# 
# We export the Study 1 chats, for the online appendix.
# 

pb <- progress::progress_bar$new(
  format = "[:bar] :percent ETA: :eta",
  total = length(chat_data),
  clear = FALSE)

for (i in 1:length(chat_data)){
  pb$tick()
  printChat(chat = chat_data[[i]], file = paste0("../results/all_chat_exports/s1_chat_exports_", chat_data[[i]]$room_id))
}

# 
# 
## Rate of Recruitment
# 
# Note: Study 1 was run starting 3/2/2022 1:41pm EST, with a target N of 500.
# 
# Plot recruitment rate
s1_ylim <- c(0,30)
s1_xlim <- s1_survey_data %>% filter(consent_response=="Yes, I want to participate") %>% pull(StartDate) %>% range %>% round_date(., "10 minutes")
border_color = "white"

pdf(file = "../results/s1_recruitment_timehist.pdf", width = 6, height = 4)

s1_survey_data %>% filter(Finished) %>% filter(consent_response=="Yes, I want to participate") %>% pull(StartDate) %>% hist(breaks = "min", freq = T, ylim = s1_ylim, xlim = s1_xlim, main = "", xlab = "Survey Start Time (EST)", ylab = "Recruits Per Minute", border = border_color, xaxt = "n")
s1_survey_data %>% filter(Finished) %>%
  filter(consent_response=="Yes, I want to participate") %>%
  filter(chat_consent) %>%
  pull(StartDate) %>% hist(breaks = "min", freq = T, col = "forestgreen", ylim = s1_ylim, xlim = s1_xlim, add = T, xaxt = "n", yaxt = "n", border = border_color)
s1_survey_data %>% filter(Finished) %>%
  filter(consent_response=="Yes, I want to participate") %>%
  filter(chat_consent) %>%
  filter(is.na(ego_code) | str_starts(ego_code, "2CC", negate = T)) %>%
  pull(StartDate) %>% hist(breaks = "min", freq = T, col = "red3", ylim = s1_ylim, xlim = s1_xlim, add = T, xaxt = "n", yaxt = "n", border = border_color)

axis(1, at = seq(s1_xlim[1], s1_xlim[2], by = "10 min"), labels = seq(s1_xlim[1], s1_xlim[2], by = "10 min") %>% strftime(., format="%H:%M", tz = "EST"))

legend("topright", legend = c("Took Survey", "Consented to Chat", "Failed to Chat"), col = c("gray", "forestgreen", "red3"), pch = 15, bty = "n")

dev.off()
# 
# 
# count consented
s1_n_consented <- s1_survey_data %>% filter(Finished) %>%
  filter(consent_response=="Yes, I want to participate") %>%
  filter(chat_consent) %>% nrow

# count chatted
s1_n_chatted <- s1_survey_data %>% filter(Finished) %>%
  filter(consent_response=="Yes, I want to participate") %>%
  filter(chat_consent) %>%
  filter(!(is.na(ego_code) | str_starts(ego_code, "2CC", negate = T))) %>% nrow

(s1_n_chatted/s1_n_consented)*100

# Export completion percentages
(s1_n_chatted/s1_n_consented) %>% {.*100} %>% round %>% paste0(., "\\%") %>% writeLines(., "../results/s1_consent_complete_perc.txt")
(1-s1_n_chatted/s1_n_consented) %>% {.*100} %>% round %>% paste0(., "\\%") %>% writeLines(., "../results/s1_consent_noncomplete_perc.txt")

# 
# 
# # Study 2
# 
# This section conducts all analyses focused on Study 2.
# 
## Prep Data
# 
### Load Data
# 
# Load the raw survey data for Study 2
s2_raw_survey_data <- readRDS(file = "../data/study_2/raw_survey_data.rds")
# Recode the treatment variable (t = 1-t) so that neutrality is the treatment, per reviewer request
s2_raw_survey_data <- s2_raw_survey_data %>% mutate(t = 1-t)
# 
# Load chat data
chat_data <- parseChat("../data/study_2/dem_chats.csv")
chat_data_republicans <- parseChat("../data/study_2/rep_chats.csv") #(a handful of republicans were recruited despite our filtering criteria; they are excluded from analysis but their data is included for completeness)
# 
chat_data <- chat_data[-(1:4)]# chats 1:4 are tests
chat_data <- c(chat_data, chat_data_republicans)
# 
# 
### Clean Data
# 
# Drop empty chats (may result from technical issues)
empty <- c()
for(i in 1:length(chat_data)){
  empty[i] <- nrow(chat_data[[i]]$messages)==0
}
if (any(empty)){chat_data <- chat_data[-which(empty)]}

# Drop "chats" where partner did not send any messages (may result from non-compliance/technical issues)
single <- c()
for(i in 1:length(chat_data)){
  single[i] <- chat_data[[i]]$messages %>% pull(participantCode) %>% {length(unique(.))<2}
}
if (any(single)){chat_data <- chat_data[-which(single)]}

# Drop chats with fewer than 2 messages from each participant (inclusion threshold specified in pre-analysis plan)
min_mess <- c()
for(i in 1:length(chat_data)){
  min_mess[i] <- chat_data[[i]]$messages %>% group_by(participantCode) %>% summarise(mess_counts = length(message)) %>% pull(mess_counts) %>% {min(.)}
}
if (any(min_mess<2)){chat_data <- chat_data[-which(min_mess<2)]}
# 
# 
# Function to recode news consumption variables
news_recode <- function(x){
  return(case_when(x == "Never" ~ 0,
          x == "Less often" ~ 1,
          x == "1-2 days a week" ~ 2,
          x == "3-6 days a week" ~ 3,
          x == "About once a day" ~ 4,
          x == "Several times a day" ~ 5))
}
# 
# 
s2_survey_data <- s2_raw_survey_data %>%
  rename(ego_code = number_entry #participant confirmation code
  ) %>%
  mutate(
    chat_consent = case_when((chat_consent=="No, I do not want to participate") ~ FALSE, #recode chat consent as binary
                             (chat_consent=="Yes, I want to participate") ~ TRUE),
    affpol_pre = case_when(partisanship == "Democrat" ~ thermo_pre_4 - thermo_pre_7,
                           partisanship == "Republican" ~ thermo_pre_7 - thermo_pre_4), #pre-chat affective polarization
    affpol_post = case_when(partisanship == "Democrat" ~ thermo_post_4 - thermo_post_7,
                            partisanship == "Republican" ~ thermo_post_7 - thermo_post_4), #post-chat affective polarization
    affpol_change = affpol_post - affpol_pre, #change in affective polarization
    demthermo_change = thermo_post_4 - thermo_pre_4, #change in thermometer ratings for Democrats
    repthermo_change = thermo_post_7 - thermo_pre_7, #change in thermometer ratings for Republicans
    PID_5 = case_when( #recode party identification to 5-point scale
      PID_dem_strength == "Strong" ~ -2,
      PID_dem_strength == "Not very strong" ~ -1,
      PID_leaners == "Closer to Democratic Party" ~ 0,
      PID_leaners == "Closer to Republican Party" ~ 0,
      PID_rep_strength == "Not very strong" ~ 1,
      PID_rep_strength == "Strong" ~ 2
    ),
    PID_6 = case_when( #recode party identification to 6-point scale
      PID_dem_strength == "Strong" ~ 1,
      PID_dem_strength == "Not very strong" ~ 2,
      PID_leaners == "Closer to Democratic Party" ~ 3,
      PID_leaners == "Closer to Republican Party" ~ 4,
      PID_rep_strength == "Not very strong" ~ 5,
      PID_rep_strength == "Strong" ~ 6
    ),
    
    id_important = case_when( #recode identity importance
      huddy1 =="Not important at all" ~ 0,
      huddy1 =="Not very important" ~ 1,
      huddy1 =="Very important" ~ 2,
      huddy1 =="Extremely important" ~ 3
    ),
    id_describe = case_when( #recode identity descriptiveness
      huddy2 =="Not at all" ~ 0,
      huddy2 =="Not very well" ~ 1,
      huddy2 =="Very well" ~ 2,
      huddy2 =="Extremely well" ~ 3
    ),
    id_wethey = case_when( #recode identity we-they
      huddy3 =="Never" ~ 0,
      huddy3 =="Rarely" ~ 1,
      huddy3 =="Some of the time" ~ 2,
      huddy3 =="Most of the time" ~ 3,
      huddy3 =="All of the time" ~ 4
    ),
    
    # recode media consumption variables
    tv = news_recode(news_media_1),
    newspapers = news_recode(news_media_2),
    radio = news_recode(news_media_3),
    internet_sm = news_recode(news_media_4),
    discussions = news_recode(news_media_5),
    podcasts = news_recode(news_media_6),
    
    # recode social media usage variables
    twitter = !is.na(social_media_1),
    facebook = !is.na(social_media_2),
    instagram = !is.na(social_media_3),
    reddit = !is.na(social_media_9),
    tiktok = !is.na(social_media_10),
    youtube = !is.na(social_media_11),
    
    age = (2021 - birthyr), #calculate age
    
    gender = fct_case_when( #recode gender
      gender == "Female" ~ "Female",
      gender == "Male" ~ "Male",
      gender == "Other (please specify):" ~ "Other"
    ),
    
    male = (gender=="Male"), #binary indicator for male gender
    
    college = as.numeric(education)>=4,#binary indicator for college education
    
    expressor = !(primary_sm %in% c("I don't use any of these services to express my opinions", NA)), #binary indicator for whether one expresses opinions about politics or current events on social media
    
    sm_postpol_num = as.numeric(sm_postpol), #numeric version of self-reported frequency of political posting on social media
    
    ideo_7 = as.numeric(ideology) #numeric version of 7-point ideology scale
  ) %>% filter(ideo_7<=4) #filter out conservative participants (screening at recruitment stage let in some conservatives)

# recode guesses of chat partner's ideology
s2_survey_data$guess_partner_ideo <- as.numeric(s2_survey_data$guess_partner_ideo)
s2_survey_data$guess_partner_ideo[which(s2_survey_data$guess_partner_ideo==8)] <- NA

# recode enjoyment of chat
s2_survey_data$enjoyment <- as.numeric(s2_survey_data$enjoyment)
s2_survey_data$enjoyment[which(s2_survey_data$enjoyment==5)] <- NA
s2_survey_data$enjoyment <- max(s2_survey_data$enjoyment, na.rm = T)-s2_survey_data$enjoyment

# clean up posting variable
s2_survey_data$sm_postpol_num[which(s2_survey_data$sm_postpol_num == 6)] <- NA
s2_survey_data$sm_postpol_num_nazero <- s2_survey_data$sm_postpol_num
s2_survey_data$sm_postpol_num_nazero[which(is.na(s2_survey_data$sm_postpol_num_nazero))] <- 0

#notes:
#thermo_pre_4 is Democrats
#thermo_pre_7 is Republicans

s2_survey_data <- s2_survey_data %>% filter(partisanship!="Republican") #filter out Republicans (screening at recruitment stage let in some Republicans)
# 
# Make a variable for "ideological extremity" as requested by reviewer
s2_survey_data <- s2_survey_data %>% mutate(ideo_ext = 4-ideo_7)

# 
# 
extro_cols <- s2_survey_data %>% dplyr::select(ends_with("_extroversion"))
extro_cols <- apply(extro_cols, MARGIN = 2, ag_disag_recode)

# PCA to create extroversion scale
extro_scale <- psych::principal(extro_cols, nfactors = 1, rotate = "varimax", missing=TRUE, impute = "mean")$scores[,1]

if (cor(extro_scale, extro_cols[,5], use = "complete.obs")<0) {
  extro_scale <- (-extro_scale)
}
s2_survey_data$extro_scale <- extro_scale
# 
# 
# PCA to create identity scale
identity_scale <- psych::principal(s2_survey_data %>% dplyr::select(starts_with("id_")), 
                              nfactors = 1, rotate = "varimax", missing=TRUE, impute = "mean")$scores[,1]

if (cor(identity_scale, s2_survey_data$id_important)<0) {
  identity_scale <- (-identity_scale)
}
s2_survey_data$identity_scale <- identity_scale
# 
# 
# PCA to create social media scale
sm_scale <- psych::principal(s2_survey_data %>% dplyr::select(twitter, facebook, instagram, reddit, tiktok, youtube),
                              nfactors = 1, rotate = "varimax", missing=TRUE, impute = "mean")$scores[,1]
s2_survey_data$sm_scale <- sm_scale
# 
# 
# PCA to create media scale
media_scale <- psych::principal(s2_survey_data %>% dplyr::select(tv, newspapers, radio, internet_sm, discussions, podcasts),
                              nfactors = 1, rotate = "varimax", missing=TRUE, impute = "mean")$scores[,1]
s2_survey_data$media_scale <- media_scale
# 
# 
# PCA to create self-monitoring scale
selfmon_scale <- psych::principal(s2_survey_data %>% transmute(as.numeric(selfMon1), as.numeric(selfMon2), as.numeric(selfMon3)),
                              nfactors = 1, rotate = "varimax", missing=TRUE, impute = "mean")$scores[,1]

if(cor(selfmon_scale, as.numeric(s2_survey_data$selfMon2))>0){selfmon_scale <- (-selfmon_scale)}
s2_survey_data$selfmon_scale <- selfmon_scale
# 
# 
# PCA to create combined extroversion-self-monitoring scale
extroselfmon_scale <- psych::principal(cbind(extro_cols, s2_survey_data %>% transmute(as.numeric(selfMon1), as.numeric(selfMon2), as.numeric(selfMon3))),
                              nfactors = 1, rotate = "varimax", missing=TRUE, impute = "mean")$scores[,1]

if (cor(extroselfmon_scale, extro_cols[,5], use = "complete.obs")<0) {
  extroselfmon_scale <- (-extroselfmon_scale)
}
s2_survey_data$extroselfmon_scale <- extroselfmon_scale
# 
# 
## Analyses
# 
### Summarize Sample Size
# 
# s2_survey_data %>% filter(consent_response=="Yes, I want to participate") %>% nrow #629
# s2_survey_data %>% filter(consent_response=="Yes, I want to participate") %>% filter(chat_consent) %>% nrow #483
# s2_survey_data %>% filter(consent_response=="Yes, I want to participate") %>% filter(chat_consent) %>% filter(!is.na(ego_code)) %>% nrow #438
# 
# 
### Loquaciousness Analyses
# 
# Function to count words in a message
word_counter <- function(message){
  return(str_count(message, "\\w+"))
}

chat_data <- featurizeChat(chat_data, featurization_function = nchar) #featurize chat data: character count
chat_data <- featurizeChat(chat_data, featurization_function = word_counter) #featurize chat data: word count

# Summarize features and add to survey data
s2_survey_data <- summarizeChat(s2_survey_data, chat_data, chat_feature_name = "nchar", summary_function = sum, na.rm = T, summary_feature_name = "char_count")
s2_survey_data <- summarizeChat(s2_survey_data, chat_data, chat_feature_name = "word_counter", summary_function = sum, na.rm = T, summary_feature_name = "word_count")
# 
# 
# Functions to calculate other loquaciousness metrics
more_thoroughness <- function(chat){
  #link from participant codd to color code
  out <- chat$messages %>%
    group_by(participantCode) %>%
    summarize(message_count = length(message),
              unique_word_count = length(unique(unlist(strsplit(message, " ")))),
              median_word_length = median(nchar(unlist(strsplit(message, " ")))),
              mean_word_length = mean(nchar(unlist(strsplit(message, " "))))
              )
  out <- chat$participants %>% dplyr::select(receiptCode, colorCode) %>% left_join(., out, by = c("colorCode" = "participantCode"))
  out <- out %>% mutate(room_id = chat$room_id, started_at = chat$started_at + 12*60*60) #manually adjust for 24-hour format
  return(out)
}

more_thoroughnesses <- function(chats, simplify = TRUE){
  out <- lapply(chats, more_thoroughness)
  if (simplify) {out <- data.table::rbindlist(l = out)}
  return(out)
}

# Implement above functions and append results to survey data
more_thoroughnesses <- more_thoroughnesses(chat_data)

s2_survey_data <- left_join(s2_survey_data %>% filter(!is.na(ego_code)), more_thoroughnesses %>% filter(!is.na(receiptCode)) %>% dplyr::select(-c("colorCode", "started_at")), by = c("ego_code" = "receiptCode"))

# 
# 
#### Do Sentiment Analysis
# 
# Custom wrapper for sentimentr sentiment analysis function
mean_sentiment <- function(message){
  this_sentiment <- sentimentr::sentiment(message) %>% group_by(element_id) %>% summarise(sentiment = mean(sentiment)) %>% pull(sentiment)
  return(this_sentiment)
}

# Apply sentiment analysis to chat data
chat_data <- featurizeChat(chat_data, feature_name = "sentiment", featurization_function = mean_sentiment)
s2_survey_data <- summarizeChat(s2_survey_data, chat_data, chat_feature_name = "sentiment", summary_function = mean, na.rm = T, summary_feature_name = "ego_mean_sentiment")
# 
# 
#### Match to Alters
# 
# For each participant, find their partner, and extract partner's variables of interest
s2_survey_data <- matchAlters(s2_survey_data, chat_data)
s2_survey_data <- getAlterVars(s2_survey_data, var_names = c("t", "ideo_7", "affpol_pre", "identity_scale", "male", "PID_6", "message_count", "char_count", "ego_mean_sentiment"))
# 
# 
### DV: Loquaciousness
# 
# Generate and export summary statistics for character count
s2_char_count_mean <- s2_survey_data %>% filter(!is.na(room_id)) %>% pull(char_count) %>% mean %>% round
s2_char_count_sd <- s2_survey_data %>% filter(!is.na(room_id)) %>% pull(char_count) %>% sd %>% round

paste0("(Mean $=",s2_char_count_mean,"$, SD $=", s2_char_count_sd,"$)") %>% writeLines(., "../results/s2_charcount_meansd.txt")

# 
# 
# Estimate pre-registered model, and variations with different operationalizations of loquaciousness
s2_prereg_base <- lm(char_count ~ t + ideo_ext + male, data = s2_survey_data %>% filter(!is.na(room_id)))
s2_prereg_mcount <- lm(message_count ~ t + ideo_ext + male, data = s2_survey_data %>% filter(!is.na(room_id)))
s2_prereg_wcount <- lm(word_count ~ t + ideo_ext + male, data = s2_survey_data %>% filter(!is.na(room_id)))
s2_prereg_uwcount <- lm(unique_word_count ~ t + ideo_ext + male, data = s2_survey_data %>% filter(!is.na(room_id)))
# 
# 
# Export regression tables
resizebox.stargazer(s2_prereg_base, s2_prereg_mcount, s2_prereg_wcount, s2_prereg_uwcount, type = "latex", title = "Loquaciousness (Study 2)",
                     se = list(make_robust_se(s2_prereg_base),make_robust_se(s2_prereg_mcount),make_robust_se(s2_prereg_wcount),make_robust_se(s2_prereg_uwcount)), label = "tab:s2_prereg_loquaciousness", out = "../results/s2_prereg_loquaciousness.tex", table.placement = "H", single.row = TRUE,
                    covariate.labels = c("Treatment\\\\(Neutral Accountability)", "Ideological Extremity", "Male"),
                    dep.var.labels = c("\\# Characters (Pre-Reg)", "\\# Messages", "\\# Words", "\\# Unique Words")
                    )
# 
# 
# Export analysis N to include in paper
length(s2_prereg_base$residuals) %>% as.character() %>% writeLines(., "../results/s2_analysis_n.txt")
# 
# 
#### Export P Values
# 
# export p-values to put in paper
summary(s2_prereg_base)$coefficients[3,4] %>% round(digits = 3) %>% str_replace(., "0.", "($p=.") %>% paste0(.,"$)") %>% writeLines(., "../results/s2_extro_p.txt")

# 
# 
# Extract p-values (based on robust SEs)

pvals_raw_s2_prereg_base <- stargazer(s2_prereg_base, type = "text", title = "Loquaciousness (Study 2)",
                     se = list(make_robust_se(s2_prereg_base)), report = "vp",
                    covariate.labels = c("Treatment\\\\(Neutral Accountability)", "Ideological Extremity", "Male"),
                    dep.var.labels = c("\\# Characters (Pre-Reg)")
                    )


pvals_raw_s2_prereg_mcount <- stargazer(s2_prereg_mcount, type = "text", title = "Loquaciousness (Study 2)",
                     se = list(make_robust_se(s2_prereg_mcount)), report = "vp",
                    covariate.labels = c("Treatment\\\\(Neutral Accountability)", "Ideological Extremity", "Male"),
                    dep.var.labels = c("\\# Messages")
                    )


# 
# 
# Export char count p-values to put in paper
pvals_raw_s2_prereg_base[which(str_detect(pvals_raw_s2_prereg_base, "TreatmentAccountability"))] %>% gsub(x = ., pattern = ".*p", replacement="p") %>% gsub(x = ., pattern = " ", replacement = "") %>% paste0("($", . , "$)") %>% writeLines(., "../results/s2_char_treatment_p.txt")
pvals_raw_s2_prereg_base[which(str_detect(pvals_raw_s2_prereg_base, "Ideological Extremity"))] %>% gsub(x = ., pattern = ".*p", replacement="p") %>% gsub(x = ., pattern = " ", replacement = "") %>% paste0("($", . , "$)") %>% writeLines(., "../results/s2_char_ideoext_p.txt")
pvals_raw_s2_prereg_base[which(str_detect(pvals_raw_s2_prereg_base, "Male"))] %>% gsub(x = ., pattern = ".*p", replacement="p") %>% gsub(x = ., pattern = " ", replacement = "") %>% paste0("($", . , "$)") %>% writeLines(., "../results/s2_char_male_p.txt")
# 
# 
# Export message count p-values to put in paper
pvals_raw_s2_prereg_mcount[which(str_detect(pvals_raw_s2_prereg_mcount, "TreatmentAccountability"))] %>% gsub(x = ., pattern = ".*p", replacement="p") %>% gsub(x = ., pattern = " ", replacement = "") %>% paste0("($", . , "$)") %>% writeLines(., "../results/s2_mess_treatment_p.txt")
pvals_raw_s2_prereg_mcount[which(str_detect(pvals_raw_s2_prereg_mcount, "Ideological Extremity"))] %>% gsub(x = ., pattern = ".*p", replacement="p") %>% gsub(x = ., pattern = " ", replacement = "") %>% paste0("($", . , "$)") %>% writeLines(., "../results/s2_mess_ideoext_p.txt")
pvals_raw_s2_prereg_mcount[which(str_detect(pvals_raw_s2_prereg_mcount, "Male"))] %>% gsub(x = ., pattern = ".*p", replacement="p") %>% gsub(x = ., pattern = " ", replacement = "") %>% paste0("($", . , "$)") %>% writeLines(., "../results/s2_mess_male_p.txt")
# 
# 
#### Coefplots
# 
# Export coefplots for appendix
pdf(file = "../results/s2_coefplot_char.pdf", width = 4, height = 4)
var_lab_cex = .8

#coefplot for pre-reg char count model
these_coefs <- summary(s2_prereg_base)$coefficients
these_ses <- make_robust_se(s2_prereg_base)

vars <- match(c("t", "ideo_ext", "maleTRUE"), rownames(these_coefs)) %>% rev
clean_names <- c("Neutral Accountability Treatment", "Ideological Extremity", "Male") %>% rev

par(bty = "n", pty = "s")
plot(x = these_coefs[vars,1], y = 1:length(vars), xlab = "Coefficient", yaxt = "n", ylab = "", type = "n",
     xlim = c(-150, 100),
     ylim = c(0,length(vars)+1),
     main = "Character Count")
abline(v = 0, lty = 3)
points(x = these_coefs[vars,1], y = 1:length(vars), pch = 16)
segments(x0 = these_coefs[vars,1]-1.96*these_ses[vars], x1 = these_coefs[vars,1]+1.96*these_ses[vars], y0 = 1:length(vars))
par(xpd = T)
text(x = these_coefs[vars,1], y = 1:length(vars), labels = clean_names, pos = 3, cex = var_lab_cex)
par(xpd = F)
dev.off()

pdf(file = "../results/s2_coefplot_mess.pdf", width = 4, height = 4)

#coefplot for exploratory message-count model
these_coefs <- summary(s2_prereg_mcount)$coefficients
these_ses <- make_robust_se(s2_prereg_mcount)

vars <- match(c("t", "ideo_ext", "maleTRUE"), rownames(these_coefs)) %>% rev
clean_names <- c("Neutral Accountability Treatment", "Ideological Extremity", "Male") %>% rev

par(bty = "n", pty = "s")
plot(x = these_coefs[vars,1], y = 1:length(vars), xlab = "Coefficient", yaxt = "n", ylab = "", type = "n",
     xlim = c(-5, 3),
     ylim = c(0,length(vars)+1),
     main = "Message Count")
abline(v = 0, lty = 3)
points(x = these_coefs[vars,1], y = 1:length(vars), pch = 16)
segments(x0 = these_coefs[vars,1]-1.96*these_ses[vars], x1 = these_coefs[vars,1]+1.96*these_ses[vars], y0 = 1:length(vars))
par(xpd = T)
text(x = these_coefs[vars,1], y = 1:length(vars), labels = clean_names, pos = 3, cex = var_lab_cex)
par(xpd = F)
dev.off()
# 
### Check Randomization
# 
# Export percentage of participants who completed the study who received the neutral treatment
s2_survey_data %>% filter(consent_response=="Yes, I want to participate") %>% filter(chat_consent) %>% filter(!is.na(ego_code)) %>% pull(t) %>% mean(., na.rm = T) %>% {.*100} %>% round %>% paste0(., "\\%") %>% writeLines(., con = "../results/s2_completed_treatment_percentage.txt")
# 
# 
# Estimate association between treatment and dropout

s2_lm_treatment_dropout <- lm(!is.na(room_id) ~ t, data = left_join(s2_raw_survey_data %>% select(ResponseId, t), s2_survey_data))

# Export p-value to put in paper
summary(s2_lm_treatment_dropout)$coefficients[2,4] %>% round(digits = 3) %>% str_replace(., "0.", "($p=.") %>% paste0(.,"$)") %>% writeLines(., "../results/s2_treatmentattrition_p.txt")

#summary(log_ex_base)$coefficients[2,4]

# 
# 
# Export Ns for treatment and control, to include in paper
s2_assignment_table <- s2_survey_data %>% filter(!is.na(room_id)) %>% group_by(t) %>% summarise(total = n())

s2_assignment_table %>% filter(t==0) %>% pull(total) %>% as.character %>% writeLines(., "../results/s2_t0_n.txt")
s2_assignment_table %>% filter(t==1) %>% pull(total) %>% as.character %>% writeLines(., "../results/s2_t1_n.txt")

# 
# 
# 
### Lee (2009) Bounds
# 
# In light of concerns about randomization, we estimated Lee (2009) bounds.  We adapted the Stata procedure recommended in the resource the reviewer provided, implementing a version of it in R, and estimating bounds for the accountability treatment effect on all four outcomes.
# 
# 
# Prep data for Lee bounds analysis
lee_bounds_data <- s2_survey_data %>% filter(chat_consent)
lee_bounds_data_2 <- left_join(lee_bounds_data, s2_survey_data) %>% mutate(observed = (!is.na(room_id)), treated = t)
# 
# 
# Analyze attrition rates
attrition_rates <- lee_bounds_data_2 %>%
  group_by(treated) %>%
  summarise(
    AttritionRate = mean(!observed),
    Count = n()
  )

# Custom function to trim data for estimating Lee bounds
trim_data <- function(data, outcome_col, group_to_trim, trim_top) {
  # Arrange and rank data based on the outcome column, conditionally on trimming direction
  group_not_to_trim_df = data %>% filter(observed == 1) %>% filter(treated != group_to_trim)
  group_to_trim_df = data %>% filter(observed == 1) %>% filter(treated == group_to_trim)
  
  if (trim_top){
    group_to_trim_df <- group_to_trim_df %>% arrange(desc(!!sym(outcome_col)))
  } else {
    group_to_trim_df <- group_to_trim_df %>% arrange((!!sym(outcome_col)))
  }
  
  group_to_trim_df_trimmed <- group_to_trim_df[1:nrow(group_not_to_trim_df),]
  
  trimmed_data <- rbind(group_not_to_trim_df, group_to_trim_df_trimmed)
  return(trimmed_data)
}

group_to_trim = ifelse(attrition_rates$AttritionRate[2] > attrition_rates$AttritionRate[1], 0, 1)

# Implementing function on data

trimmed_char_lower <- trim_data(lee_bounds_data_2, "char_count", group_to_trim = group_to_trim, trim_top = FALSE)
trimmed_char_upper <- trim_data(lee_bounds_data_2, "char_count", group_to_trim = group_to_trim, trim_top = TRUE)

trimmed_mcount_lower <- trim_data(lee_bounds_data_2, "message_count", group_to_trim = group_to_trim, trim_top = FALSE)
trimmed_mcount_upper <- trim_data(lee_bounds_data_2, "message_count", group_to_trim = group_to_trim, trim_top = TRUE)

trimmed_wcount_lower <- trim_data(lee_bounds_data_2, "word_count", group_to_trim = group_to_trim, trim_top = FALSE)
trimmed_wcount_upper <- trim_data(lee_bounds_data_2, "word_count", group_to_trim = group_to_trim, trim_top = TRUE)

trimmed_uwcount_lower <- trim_data(lee_bounds_data_2, "unique_word_count", group_to_trim = group_to_trim, trim_top = FALSE)
trimmed_uwcount_upper <- trim_data(lee_bounds_data_2, "unique_word_count", group_to_trim = group_to_trim, trim_top = TRUE)
# 
# 
# Estimating models on trimmed datasets

s2_prereg_base_lee_lower <- lm(char_count ~ t + ideo_7 + male, data = trimmed_char_lower %>% filter(!is.na(room_id)))
s2_prereg_base_lee_upper <- lm(char_count ~ t + ideo_7 + male, data = trimmed_char_upper %>% filter(!is.na(room_id)))

s2_prereg_mcount_lee_lower <- lm(message_count ~ t + ideo_7 + male, data = trimmed_mcount_lower %>% filter(!is.na(room_id)))
s2_prereg_mcount_lee_upper <- lm(message_count ~ t + ideo_7 + male, data = trimmed_mcount_upper %>% filter(!is.na(room_id)))

s2_prereg_wcount_lee_lower <- lm(word_count ~ t + ideo_7 + male, data = trimmed_wcount_lower %>% filter(!is.na(room_id)))
s2_prereg_wcount_lee_upper <- lm(word_count ~ t + ideo_7 + male, data = trimmed_wcount_upper %>% filter(!is.na(room_id)))

s2_prereg_uwcount_lee_lower <- lm(unique_word_count ~ t + ideo_7 + male, data = trimmed_uwcount_lower %>% filter(!is.na(room_id)))
s2_prereg_uwcount_lee_upper <- lm(unique_word_count ~ t + ideo_7 + male, data = trimmed_uwcount_upper %>% filter(!is.na(room_id)))

# Extracting summary info from estimated models
s2_prereg_base_lee_lower_est <- summary(s2_prereg_base_lee_lower)$coefficients %>% {.[which(rownames(.) == "t"), which(colnames(.) == "Estimate")]}
s2_prereg_base_lee_upper_est <- summary(s2_prereg_base_lee_upper)$coefficients %>% {.[which(rownames(.) == "t"), which(colnames(.) == "Estimate")]}

s2_prereg_mcount_lee_lower_est <- summary(s2_prereg_mcount_lee_lower)$coefficients %>% {.[which(rownames(.) == "t"), which(colnames(.) == "Estimate")]}
s2_prereg_mcount_lee_upper_est <- summary(s2_prereg_mcount_lee_upper)$coefficients %>% {.[which(rownames(.) == "t"), which(colnames(.) == "Estimate")]}

s2_prereg_wcount_lee_lower_est <- summary(s2_prereg_wcount_lee_lower)$coefficients %>% {.[which(rownames(.) == "t"), which(colnames(.) == "Estimate")]}
s2_prereg_wcount_lee_upper_est <- summary(s2_prereg_wcount_lee_upper)$coefficients %>% {.[which(rownames(.) == "t"), which(colnames(.) == "Estimate")]}

s2_prereg_uwcount_lee_lower_est <- summary(s2_prereg_uwcount_lee_lower)$coefficients %>% {.[which(rownames(.) == "t"), which(colnames(.) == "Estimate")]}
s2_prereg_uwcount_lee_upper_est <- summary(s2_prereg_uwcount_lee_upper)$coefficients %>% {.[which(rownames(.) == "t"), which(colnames(.) == "Estimate")]}

# 
# 
# Extracting estimates for plotting

s2_prereg_base_est <- summary(s2_prereg_base)$coefficients %>% {.[which(rownames(.) == "t"), which(colnames(.) == "Estimate")]}
s2_prereg_mcount_est <- summary(s2_prereg_mcount)$coefficients %>% {.[which(rownames(.) == "t"), which(colnames(.) == "Estimate")]}
s2_prereg_wcount_est <- summary(s2_prereg_wcount)$coefficients %>% {.[which(rownames(.) == "t"), which(colnames(.) == "Estimate")]}
s2_prereg_uwcount_est <- summary(s2_prereg_uwcount)$coefficients %>% {.[which(rownames(.) == "t"), which(colnames(.) == "Estimate")]}

# 
# 
# Plot Lee (2009) Bounds

pdf("../results/lee_bounds.pdf", width = 8, height = 4)

par(mar = c(5.1, 4.1, 4.1, 8.1))
#c(5.1, 4.1, 4.1, 2.1)

plot(y = c(1,1,2,2,3,3,4,4), ylim = c(0,5), xlab = "Estimate", ylab = "", yaxt = "n", main = "Lee (2009) Treatment Effect Bounds",
     x = c(s2_prereg_base_lee_lower_est, s2_prereg_base_lee_upper_est,
           s2_prereg_mcount_lee_lower_est, s2_prereg_mcount_lee_upper_est,
           s2_prereg_wcount_lee_lower_est, s2_prereg_wcount_lee_upper_est,
           s2_prereg_uwcount_lee_lower_est, s2_prereg_uwcount_lee_upper_est), pch = 4)
points(y = 1:4, x = c(s2_prereg_base_est, s2_prereg_mcount_est, s2_prereg_wcount_est, s2_prereg_uwcount_est), pch = 16)
segments(y0 = 1:4, x0 = c(s2_prereg_base_lee_lower_est, s2_prereg_mcount_lee_lower_est, s2_prereg_wcount_lee_lower_est, s2_prereg_uwcount_lee_lower_est), x1 = c(s2_prereg_base_lee_upper_est, s2_prereg_mcount_lee_upper_est, s2_prereg_wcount_lee_upper_est, s2_prereg_uwcount_lee_upper_est))

abline(v = 0, lty = 2)

par(xpd=T)

#text(y = 1:4, x = rep(140,4), labels = c("# Characters", "# Messages", "# Words", "# Unique Words"), pos = 4)
text(y = 1:4, x = rep(s2_prereg_base_lee_lower_est+10,4), labels = c("# Characters", "# Messages", "# Words", "# Unique Words"), pos = 4)

dev.off()
# 
# 
### DV: Affective Polarization
# 
# Analyze whether attitudes became more polarized after (vs before) the chat

# Prep dataframe
affpol_data <- s2_survey_data %>% filter(chat_consent) %>% filter(!is.na(affpol_pre) & !is.na(affpol_post))
affpol_data <- affpol_data %>% mutate(message_diff = (alter_message_count - message_count), char_diff = (alter_char_count - char_count)) %>% mutate(message_diff_bin = as.numeric(message_diff>0), char_diff_bin = as.numeric(char_diff>0))

# Estimate T-Test for change in affective polarization
affpol_change_ttest <- t.test(affpol_data$affpol_post, affpol_data$affpol_pre, paired = TRUE, alternative = "two.sided")

# Export mean change for inclusion in paper
affpol_change_ttest$estimate %>% round(., 1) %>% as.character %>% writeLines(., "../results/s2_mean_affpol_change.txt")

# Export significance for inclusion in paper
print_sig_level(x = affpol_change_ttest$p.value, max_places = 4) %>% paste0("($", ., "$, 2-tailed paired t-test)") %>% writeLines("../results/s2_affpol_change_ttest_p.txt")
# 
# 
# 
# Estimate linear regression models to characterize group polarization
lm_affpol_0 <- lm(affpol_post ~ affpol_pre, data = affpol_data %>% filter(!is.na(room_id)))
lm_affpol_1 <- lm(affpol_post ~ affpol_pre + alter_affpol_pre, data = affpol_data %>% filter(!is.na(room_id)))
lm_affpol_1a <- lm(affpol_post ~ affpol_pre + alter_affpol_pre*message_diff_bin, data = affpol_data %>% filter(!is.na(room_id)))
lm_affpol_1b <- lm(affpol_post ~ affpol_pre + alter_affpol_pre*char_diff_bin, data = affpol_data %>% filter(!is.na(room_id)))
lm_affpol_2 <- lm(affpol_post ~ affpol_pre + alter_affpol_pre + t, data = affpol_data %>% filter(!is.na(room_id)))

lm_affpol_3a <- lm(affpol_post ~ affpol_pre + t*alter_affpol_pre, data = affpol_data %>% filter(!is.na(room_id)))
lm_affpol_3b <- lm(affpol_post ~ affpol_pre + t*pid_str, data = affpol_data %>% filter(!is.na(room_id)) %>% mutate(pid_str = (4-PID_6)))

# Export regression tables
stargazer(lm_affpol_0,lm_affpol_1,lm_affpol_2, type = "latex",
          se = list(make_robust_se(lm_affpol_0), make_robust_se(lm_affpol_1), make_robust_se(lm_affpol_2)), title = "Group Polarization", label = "tab:s2_group_polarization_main", out = "../results/s2_group_polarization_main.tex", dep.var.labels = "Post-Chat Partisan Affect (D - R Feeling Thermo Diff)", covariate.labels = c("Pre-Chat Affect", "Partner's Pre-Chat Affect", "Treatment\\\\(Neutral Accountability)"), omit.stat = c("f"), table.placement = "H", single.row = T)


stargazer(lm_affpol_1a,lm_affpol_1b, type = "latex",
         se = list(make_robust_se(lm_affpol_1a), make_robust_se(lm_affpol_1b)), title = "Group Polarization (Interaction Tests of Relative Loquaciousness)", label = "tab:s2_group_polarization_interactions", omit.stat = c("f"), out = "../results/s2_group_polarization_interactions.tex", table.placement = "H", single.row = T,
         dep.var.labels = "Post-Chat Partisan Affect (D - R Feeling Thermo Diff)", covariate.labels = c("Pre-Chat Affect", "Partner's Pre-Chat Affect", "Partner Sent More Messages", "Partner's Pre-Chat Affect\\\\$\\times$ Partner Sent More Messages", "Partner Wrote More Characters", "Partner's Pre-Chat Affect\\\\$\\times$ Partner Wrote More Characters"))


stargazer(lm_affpol_3a, lm_affpol_3b, type = "latex",
          se = list(make_robust_se(lm_affpol_3a), make_robust_se(lm_affpol_3b)),
           title = "Group Polarization (Interaction Tests with Accountability Treatment)",
          out = "../results/s2_group_polarization_interactions_withtreat.tex",
          label = "tab:s2_group_polarization_interactions_withtreat", omit.stat = c("f"),
          table.placement = "H", single.row = T,
          dep.var.labels = "Post-Chat Partisan Affect (D - R Feeling Thermo Diff)",
          covariate.labels = c("Pre-Chat Affect", "Treatment\\\\(Neutral Accountability)", "Partner's Pre-Chat Affect", "Partner's Pre-Chat Affect\\\\$\\times$ Treatment", "Partisanship Strength", "Partisanship Strength\\\\$\\times$ Treatment"))

# 
# 
### Sentiment Analysis
# 
# Estimate linear regression models to characterize predictors of message sentiment
s2_mean_sentiment_base <- lm(ego_mean_sentiment ~ affpol_pre, data = affpol_data %>% filter(!is.na(alter_affpol_pre)))

s2_sentiment_lm_ks <- lm(ego_mean_sentiment ~ affpol_pre + extro_scale + selfmon_scale + as.numeric(polInt) + PID_6 + ideo_7 + identity_scale + media_scale, data = affpol_data %>% filter(!is.na(alter_affpol_pre)))

# Export regression tables
stargazer(s2_mean_sentiment_base, s2_sentiment_lm_ks, type = "latex",
         se = list(make_robust_se(s2_mean_sentiment_base), make_robust_se(s2_sentiment_lm_ks)), title = "Baseline Partisan Affect and Mean Message Sentiment", label = "tab:s2_sentiment", omit.stat = c("f"), out = "../results/s2_sentiment.tex", table.placement = "H", single.row = T,
         dep.var.labels = "Mean Sentiment of Messages", covariate.labels = c("Baseline Partisan Affect", "Extroversion", "Self-Monitoring", "Political Interest", "6-Point PID", "7-Point Ideology", "Ideological Identity", "Media Consumption Scale"))
# 
# 
## Export Chats
# 
# Export chats from Study 2 for online appendices
pb <- progress::progress_bar$new(
  format = "[:bar] :percent ETA: :eta",
  total = length(chat_data),
  clear = FALSE)

for (i in 1:length(chat_data)){
  pb$tick()
  printChat(chat = chat_data[[i]], file = paste0("../results/all_chat_exports/s2_chat_exports_", chat_data[[i]]$room_id))
}
# 
# 
## Rate of Recruitment
# 
# Notes:
# - Recruitment occurred over 3 days: May 26th, May 27th, and May 30th (29-29 was a weekend)
# - For Study 2, consent_response = chat_consent, so the plot looks different.
# 
# Plot recruitment rate for Study 2

s2_ylim <- c(0,15)

s2_recruitment_sub1 <- s2_survey_data %>% filter(StartDate < as.POSIXct("2022-05-27 00:00:00 EDT"))

s2_xlim <- s2_recruitment_sub1 %>% filter(consent_response=="Yes, I want to participate") %>% pull(StartDate) %>% range

s2_recruitment_sub1 %>% filter(Finished) %>% filter(consent_response=="Yes, I want to participate") %>% pull(StartDate) %>% hist(breaks = "min", freq = T, ylim = s2_ylim, xlim = s2_xlim, main = "Study 2 Recruitment", xlab = "Start Time", ylab = "Recruits Per Minute", border = border_color, xaxt = "n")
s2_recruitment_sub1 %>% filter(Finished) %>%
  filter(consent_response=="Yes, I want to participate") %>%
  filter(chat_consent) %>%
  pull(StartDate) %>% hist(breaks = "min", freq = T, col = "forestgreen", ylim = s2_ylim, xlim = s2_xlim, add = T, xaxt = "n", border = border_color)
s2_recruitment_sub1 %>% filter(Finished) %>%
  filter(consent_response=="Yes, I want to participate") %>%
  filter(chat_consent) %>%
  filter(is.na(ego_code) | str_starts(ego_code, "2CC", negate = T)) %>%
  pull(StartDate) %>% hist(breaks = "min", freq = T, col = "red3", ylim = s2_ylim, xlim = s2_xlim, add = T, xaxt = "n", border = border_color)
axis(1, at = seq(s2_xlim[1], s2_xlim[2], by = "10 min"), labels = seq(s2_xlim[1], s2_xlim[2], by = "10 min") %>% strftime(., format="%H:%M"))

s2_recruitment_sub2 <- s2_survey_data %>% filter(StartDate > as.POSIXct("2022-05-27 00:00:00 EDT")) %>% filter(StartDate < as.POSIXct("2022-05-29 00:00:00 EDT"))

s2_xlim <- s2_recruitment_sub2 %>% filter(consent_response=="Yes, I want to participate") %>% pull(StartDate) %>% range

s2_recruitment_sub2 %>% filter(Finished) %>% filter(consent_response=="Yes, I want to participate") %>% pull(StartDate) %>% hist(breaks = "min", freq = T, ylim = s2_ylim, xlim = s2_xlim, main = "May 27th", xlab = "Start Time", ylab = "Recruits Per Minute", border = border_color, xaxt = "n")
s2_recruitment_sub2 %>% filter(Finished) %>%
  filter(consent_response=="Yes, I want to participate") %>%
  filter(chat_consent) %>%
  pull(StartDate) %>% hist(breaks = "min", freq = T, col = "forestgreen", ylim = s2_ylim, xlim = s2_xlim, add = T, xaxt = "n", border = border_color)
s2_recruitment_sub2 %>% filter(Finished) %>%
  filter(consent_response=="Yes, I want to participate") %>%
  filter(chat_consent) %>%
  filter(is.na(ego_code) | str_starts(ego_code, "2CC", negate = T)) %>%
  pull(StartDate) %>% hist(breaks = "min", freq = T, col = "red3", ylim = s2_ylim, xlim = s2_xlim, add = T, xaxt = "n", border = border_color)
axis(1, at = seq(s2_xlim[1], s2_xlim[2], by = "10 min"), labels = seq(s2_xlim[1], s2_xlim[2], by = "10 min") %>% strftime(., format="%H:%M"))


s2_recruitment_sub3 <- s2_survey_data %>% filter(StartDate > as.POSIXct("2022-05-29 00:00:00 EDT"))

s2_xlim <- s2_recruitment_sub3 %>% filter(consent_response=="Yes, I want to participate") %>% pull(StartDate) %>% range

s2_recruitment_sub3 %>% filter(Finished) %>% filter(consent_response=="Yes, I want to participate") %>% pull(StartDate) %>% hist(breaks = "min", freq = T, ylim = s2_ylim, xlim = s2_xlim, main = "May 30th", xlab = "Start Time", ylab = "Recruits Per Minute", border = border_color, xaxt = "n")
s2_recruitment_sub3 %>% filter(Finished) %>%
  filter(consent_response=="Yes, I want to participate") %>%
  filter(chat_consent) %>%
  pull(StartDate) %>% hist(breaks = "min", freq = T, col = "forestgreen", ylim = s2_ylim, xlim = s2_xlim, add = T, xaxt = "n", border = border_color)
s2_recruitment_sub3 %>% filter(Finished) %>%
  filter(consent_response=="Yes, I want to participate") %>%
  filter(chat_consent) %>%
  filter(is.na(ego_code) | str_starts(ego_code, "2CC", negate = T)) %>%
  pull(StartDate) %>% hist(breaks = "min", freq = T, col = "red3", ylim = s2_ylim, xlim = s2_xlim, add = T, xaxt = "n", border = border_color)
axis(1, at = seq(s2_xlim[1], s2_xlim[2], by = "10 min"), labels = seq(s2_xlim[1], s2_xlim[2], by = "10 min") %>% strftime(., format="%H:%M"))

# 
# 
# 
# Export Study 2 recruitment plots
s2_ylim <- c(0,15)

s2_recruitment_sub1 <- s2_survey_data %>% filter(StartDate < as.POSIXct("2022-05-27 00:00:00 EDT"))
s2_recruitment_sub2 <- s2_survey_data %>% filter(StartDate > as.POSIXct("2022-05-27 00:00:00 EDT")) %>% filter(StartDate < as.POSIXct("2022-05-29 00:00:00 EDT"))
s2_recruitment_sub3 <- s2_survey_data %>% filter(StartDate > as.POSIXct("2022-05-29 00:00:00 EDT"))

## Export plots with different dimensions, for layout
day1_range <- s2_recruitment_sub1 %>% filter(consent_response=="Yes, I want to participate") %>% pull(StartDate) %>% range
day2_range <- s2_recruitment_sub2 %>% filter(consent_response=="Yes, I want to participate") %>% pull(StartDate) %>% range
top_left_width_fraction <- as.numeric(difftime(day1_range[2], day1_range[1]))/(as.numeric(difftime(day1_range[2], day1_range[1]))+as.numeric(difftime(day2_range[2], day2_range[1])))

top_left_width_fraction*.99
(1-top_left_width_fraction)*.99

s2_xlim <- s2_recruitment_sub1 %>% filter(consent_response=="Yes, I want to participate") %>% pull(StartDate) %>% range

pdf(file = "../results/s2_recruitment_timehist_may26.pdf", width = 12*top_left_width_fraction, height = 4)
s2_recruitment_sub1 %>% filter(Finished) %>% filter(consent_response=="Yes, I want to participate") %>% pull(StartDate) %>% hist(breaks = "min", freq = T, ylim = s2_ylim, xlim = s2_xlim, main = "May 26th", xlab = "Start Time", ylab = "Recruits Per Minute", border = border_color, xaxt = "n")
s2_recruitment_sub1 %>% filter(Finished) %>%
  filter(consent_response=="Yes, I want to participate") %>%
  filter(chat_consent) %>%
  pull(StartDate) %>% hist(breaks = "min", freq = T, col = "forestgreen", ylim = s2_ylim, xlim = s2_xlim, add = T, xaxt = "n", border = border_color)
s2_recruitment_sub1 %>% filter(Finished) %>%
  filter(consent_response=="Yes, I want to participate") %>%
  filter(chat_consent) %>%
  filter(is.na(ego_code) | str_starts(ego_code, "2CC", negate = T)) %>%
  pull(StartDate) %>% hist(breaks = "min", freq = T, col = "red3", ylim = s2_ylim, xlim = s2_xlim, add = T, xaxt = "n", border = border_color)
axis(1, at = seq(s2_xlim[1], s2_xlim[2], by = "10 min"), labels = seq(s2_xlim[1], s2_xlim[2], by = "10 min") %>% strftime(., format="%H:%M", tz = "EST"))
dev.off()

s2_xlim <- s2_recruitment_sub2 %>% filter(consent_response=="Yes, I want to participate") %>% pull(StartDate) %>% range

pdf(file = "../results/s2_recruitment_timehist_may27.pdf", width = 12*(1-top_left_width_fraction), height = 4)
s2_recruitment_sub2 %>% filter(Finished) %>% filter(consent_response=="Yes, I want to participate") %>% pull(StartDate) %>% hist(breaks = "min", freq = T, ylim = s2_ylim, xlim = s2_xlim, main = "May 27th", xlab = "Start Time", ylab = "Recruits Per Minute", border = border_color, xaxt = "n")
s2_recruitment_sub2 %>% filter(Finished) %>%
  filter(consent_response=="Yes, I want to participate") %>%
  filter(chat_consent) %>%
  pull(StartDate) %>% hist(breaks = "min", freq = T, col = "forestgreen", ylim = s2_ylim, xlim = s2_xlim, add = T, xaxt = "n", border = border_color)
s2_recruitment_sub2 %>% filter(Finished) %>%
  filter(consent_response=="Yes, I want to participate") %>%
  filter(chat_consent) %>%
  filter(is.na(ego_code) | str_starts(ego_code, "2CC", negate = T)) %>%
  pull(StartDate) %>% hist(breaks = "min", freq = T, col = "red3", ylim = s2_ylim, xlim = s2_xlim, add = T, xaxt = "n", border = border_color)
axis(1, at = seq(s2_xlim[1], s2_xlim[2], by = "10 min"), labels = seq(s2_xlim[1], s2_xlim[2], by = "10 min") %>% strftime(., format="%H:%M", tz = "EST"))
dev.off()

s2_xlim <- s2_recruitment_sub3 %>% filter(consent_response=="Yes, I want to participate") %>% pull(StartDate) %>% range

pdf(file = "../results/s2_recruitment_timehist_may30.pdf", width = 12, height = 4)
s2_recruitment_sub3 %>% filter(Finished) %>% filter(consent_response=="Yes, I want to participate") %>% pull(StartDate) %>% hist(breaks = "min", freq = T, ylim = s2_ylim, xlim = s2_xlim, main = "May 30th", xlab = "Start Time", ylab = "Recruits Per Minute", border = border_color, xaxt = "n")
s2_recruitment_sub3 %>% filter(Finished) %>%
  filter(consent_response=="Yes, I want to participate") %>%
  filter(chat_consent) %>%
  pull(StartDate) %>% hist(breaks = "min", freq = T, col = "forestgreen", ylim = s2_ylim, xlim = s2_xlim, add = T, xaxt = "n", border = border_color)
s2_recruitment_sub3 %>% filter(Finished) %>%
  filter(consent_response=="Yes, I want to participate") %>%
  filter(chat_consent) %>%
  filter(is.na(ego_code) | str_starts(ego_code, "2CC", negate = T)) %>%
  pull(StartDate) %>% hist(breaks = "min", freq = T, col = "red3", ylim = s2_ylim, xlim = s2_xlim, add = T, xaxt = "n", border = border_color)
legend("topright", legend = c("Consented to Chat", "Failed to Chat"), col = c("forestgreen", "red3"), pch = 15, bty = "n")
axis(1, at = seq(s2_xlim[1], s2_xlim[2], by = "10 min"), labels = seq(s2_xlim[1], s2_xlim[2], by = "10 min") %>% strftime(., format="%H:%M", tz = "EST"))
dev.off()


## Addressing Skew of Count Outcomes
# We address the issue of skew in the count outcomes by logging the count outcomes and re-estimating the models.

# First, plot histograms of the count outcomes to illustrate skew
pdf(file = "../results/s2_char_hist.pdf", width = 4, height = 4)
par(pty = "s")
s2_survey_data %>% filter(!is.na(room_id)) %>% pull(char_count) %>% hist(., main = "Character Count", xlab = "Character Count")
dev.off()

pdf(file = "../results/s2_mess_hist.pdf", width = 4, height = 4)
par(pty = "s")
s2_survey_data %>% filter(!is.na(room_id)) %>% pull(message_count) %>% hist(., main = "Message Count", xlab = "Message Count")
dev.off()

pdf(file = "../results/s2_word_hist.pdf", width = 4, height = 4)
par(pty = "s")
s2_survey_data %>% filter(!is.na(room_id)) %>% pull(word_count) %>% hist(., main = "Word Count", xlab = "Word Count")
dev.off()

pdf(file = "../results/s2_uniqueword_hist.pdf", width = 4, height = 4)
par(pty = "s")
s2_survey_data %>% filter(!is.na(room_id)) %>% pull(unique_word_count) %>% hist(., main = "Unique Word Count", xlab = "Unique Word Count")
dev.off()

# Next, conduct D'Agostino's test for normality
dagostino_test_charcount <- agostino.test(s2_survey_data %>% filter(!is.na(room_id)) %>% pull(char_count))
dagostino_test_messcount <- agostino.test(s2_survey_data %>% filter(!is.na(room_id)) %>% pull(message_count))
dagostino_test_wordcount <- agostino.test(s2_survey_data %>% filter(!is.na(room_id)) %>% pull(word_count))
dagostino_test_uniquewordcount <- agostino.test(s2_survey_data %>% filter(!is.na(room_id)) %>% pull(unique_word_count))

# Estimate models with logged count outcomes
s2_base_log <- lm(log(char_count) ~ t + ideo_7 + male, data = s2_survey_data %>% filter(!is.na(room_id)))
s2_mcount_log <- lm(log(message_count) ~ t + ideo_7 + male, data = s2_survey_data %>% filter(!is.na(room_id)))
s2_wcount_log <- lm(log(word_count) ~ t + ideo_7 + male, data = s2_survey_data %>% filter(!is.na(room_id)))
s2_uwcount_log <- lm(log(unique_word_count) ~ t + ideo_7 + male, data = s2_survey_data %>% filter(!is.na(room_id)))

# Export regression tables
resizebox.stargazer(s2_base_log, s2_mcount_log, s2_wcount_log, s2_uwcount_log, type = "latex",
                     se = list(make_robust_se(s2_base_log),make_robust_se(s2_mcount_log),make_robust_se(s2_wcount_log),make_robust_se(s2_uwcount_log)), label = "tab:s2_prereg_loquaciousness_logged", title = "Logged Loquaciousness (Study 2)", out = "../results/s2_prereg_loquaciousness_logged.tex", table.placement = "H", tab.width = "\\textwidth",
                    covariate.labels = c("Treatment\\\\(Neutral Accountability)", "7-Point Ideology", "Male"),
                    dep.var.labels = c("Log \\# Characters (Pre-Reg)", "Log \\# Messages", "Log \\# Words", "Log \\# Unique Words"),
                    single.row = T)

# 
# 
## Check and Report Number of Cluster
# 
# Export number of chatrooms
s2_survey_data %>% filter(!is.na(room_id)) %>% pull(room_id) %>% unique %>% length %>% as.character %>% writeLines(., "../results/s2_chatroom_n.txt")
# 
# 
# # Pooled Analyses
# 
# Here we conduct analyses that include data from both Studies 1 and 2.
# 
## PID Strength Distribution
# 
# Plot distribution of PID strength for appendix
pooled_PID_6 <- c(s1_survey_data$PID_6, s2_survey_data$PID_6)

pdf("../results/pooled_pid_dist.pdf", width = 8, height = 4)

hist(4-(pooled_PID_6), main = "Distribution of Partisanship Strength", xaxt = "n", xlab = "")
axis(1, at = c(1.1,1.9,2.9), labels = c("Lean Dem", "Not So Strong Dem", "Strong Dem"), lwd = 0)

dev.off()
# 
# 
## Examples of What People Actually Wrote
# 
# Function to get dimensions of the first page of all PDFs in a directory
get_pdf_dimensions <- function(directory_path, identifier) {
  # Ensure the pdftools package is installed and loaded
  if (!requireNamespace("pdftools", quietly = TRUE)) {
    install.packages("pdftools")
    library(pdftools)
  }
  
  # List all PDF files in the directory
  pdf_files <- list.files(path = directory_path, pattern = "\\.pdf$", full.names = TRUE)
  pdf_files <- pdf_files[which(str_detect(pdf_files, identifier))]
  
  # Initialize an empty data frame to store results
  results <- data.frame("name" = character(), "width" = numeric(), "height" = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each PDF file
  for (pdf_file in pdf_files) {
    # Render the first page of the PDF
    dims <- pdf_render_page(pdf_file, page = 1) %>% dim %>% {.[-1]}
    
    # Add the information to the results data frame
    results <- rbind(results, data.frame("name" = basename(pdf_file), "width" = dims[1], "height" = dims[2], stringsAsFactors = FALSE))
  }
  
  results <- results %>% arrange(name)
  
  return(results)
}

s1_chat_export_dims <- get_pdf_dimensions("../results/all_chat_exports/", identifier = "s1_chat_exports_")
s2_chat_export_dims <- get_pdf_dimensions("../results/all_chat_exports/", identifier = "s2_chat_exports_")
# 
# 
# Select random sample of chats for appendix
set.seed(1) # set seed for reproducibility
save_path <- "../results/"
#dir.create(save_path)

example_chats_for_appendix_s1 <- s1_chat_export_dims %>% slice_sample(n = 6)
example_chats_for_appendix_s2 <- s2_chat_export_dims %>% slice_sample(n = 6)

files_to_copy <- paste0("../results/all_chat_exports/", example_chats_for_appendix_s1$name)
for (i in seq_along(files_to_copy)) {
  file.copy(from = files_to_copy[i], to = paste0(save_path, "s1_chat_example_", i, ".pdf"))
}

files_to_copy <- paste0("../results/all_chat_exports/", example_chats_for_appendix_s2$name)
for (i in seq_along(files_to_copy)) {
  file.copy(from = files_to_copy[i], to = paste0(save_path, "s2_chat_example_", i, ".pdf"))
}

# 
# 
## Diff b/w Mobile vs Desktop
# 
# Merge data for analyzing differences between mobile and desktop users
s12_merged_survey_data <- bind_rows(s1_survey_data_nona %>% mutate(study=1) %>% mutate(gender = case_when(gender == "Male" ~ "Male",
                                                                                 gender == "Female" ~ "Female",
                                                                                 gender == "Other (please specify):" ~ "Other")),
                                    s2_survey_data %>% mutate(study=2))

# Use operating system names to identify mobile vs desktop
os_mobile <- "Android|iPad|iPhone"
os_computer <- "Macintosh|Windows NT|CrOS|Linux|X11|Ubuntu"

s12_merged_survey_data <- s12_merged_survey_data %>% mutate(took_on_mobile = case_when(str_detect(string = `browser_info_Operating System`, pattern = os_mobile) ~ T,
                                                                                       str_detect(string = `browser_info_Operating System`, pattern = os_computer) ~ F))
# 
# 
# Export prevalence of mobile usage
s12_merged_survey_data %>% filter(!is.na(took_on_mobile)) %>% pull(took_on_mobile) %>% mean %>% {.*100} %>% round %>% paste0(., "\\%") %>% writeLines(., con = "../results/s12_mobile_percentage.txt")
# 
# 
# Estimate regressions to characterize association between mobile and loquaciousness
lm_mob_char <- lm(char_count ~ took_on_mobile + as.factor(study), data = s12_merged_survey_data)
lm_mob_mess <- lm(message_count ~ took_on_mobile + as.factor(study), data = s12_merged_survey_data)
lm_mob_word <- lm(word_count ~ took_on_mobile + as.factor(study), data = s12_merged_survey_data)
lm_mob_uniqueword <- lm(unique_word_count ~ took_on_mobile + as.factor(study), data = s12_merged_survey_data)

# 
# 
# Export regression tables
stargazer(lm_mob_char, lm_mob_mess, lm_mob_word, lm_mob_uniqueword, type = "latex", label = "tab:s12_mobile_loq", title = "Mobile Devices and Loquaciousness", out = "../results/s12_mobile_loq.tex", table.placement = "H",
          #tab.width = "textwidth",
          single.row = T)

resizebox.stargazer(lm_mob_char, lm_mob_mess, lm_mob_word, lm_mob_uniqueword, 
                    type = "latex", label = "tab:s12_mobile_loq", title = "Mobile Devices and Loquaciousness",
                    out = "../results/s12_mobile_loq.tex",
                    table.placement = "H",
                    tab.width = "\\textwidth", omit.stat = c("f"), single.row = TRUE)

# 
# 
# Merge data for analyzing technical issues in mobile vs desktop users
s12_merged_survey_data2 <- bind_rows(s1_survey_data %>% mutate(study=1) %>% mutate(gender = case_when(gender == "Male" ~ "Male",
                                                                                 gender == "Female" ~ "Female",
                                                                                 gender == "Other (please specify):" ~ "Other")),
                                    s2_survey_data %>% mutate(study=2))


os_mobile <- "Android|iPad|iPhone"
os_computer <- "Macintosh|Windows NT|CrOS|Linux|X11|Ubuntu"

s12_merged_survey_data2 <- s12_merged_survey_data2 %>% mutate(took_on_mobile = case_when(str_detect(string = `browser_info_Operating System`, pattern = os_mobile) ~ T,
                                                                                         str_detect(string = `browser_info_Operating System`, pattern = os_computer) ~ F))
# Estimate regression to characterize association between mobile and technical problems
lm_mob_techproblem <- lm(!is.na(had_problem_1) ~ took_on_mobile + as.factor(study) + age, data = s12_merged_survey_data)

# Export regression table
stargazer(lm_mob_techproblem, type = "latex", label = "tab:s12_mobile_tech", title = "Mobile Devices and Technical Problems Reported", out = "../results/s12_mobile_tech.tex", table.placement = "H",
          #tab.width = "textwidth",
          single.row = T)
# 
# 
# 
# 