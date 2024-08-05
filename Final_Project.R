install.packages('modelsummary')
library(modelsummary)
library(forecast)
library(tidyverse)
library(caret)
library(ROCR)
library(broom)
library(arules)

## Logistic Regression
Spotify <-read_csv("Downloads/Popular_Spotify_Songs (1).csv")
Spotify <- Spotify[-575,]
Spotify$streams <- as.numeric(Spotify$streams)
Spotify <- na.omit(Spotify)
Spotify$key <- factor(Spotify$key)
Spotify$mode <- factor(Spotify$mode)
Spotify$artist_count <- factor(Spotify$artist_count)
Spotify$released_year <- factor(Spotify$released_year)
Spotify$released_month <- factor(Spotify$released_month)
Spotify$released_day <- factor(Spotify$released_day)
mean(Spotify$streams)

# 468985764 = average number of streams

AvgStreams = 468985764
Spotify <- Spotify %>% mutate(IsPopular = ifelse(streams >= AvgStreams,1,0))

Spotify$Obs <- seq_len(nrow(Spotify))

set.seed(1)
Train <- Spotify %>% sample_frac(0.7)
Validation <- Spotify %>% anti_join(Train, by = "Obs")

names(Train)[names(Train) == 'artist(s)_name'] <- 'artist_s_name'
names(Spotify)[names(Spotify) == 'artist(s)_name'] <- 'artist_s_name'
names(Validation)[names(Validation) == 'artist(s)_name'] <- 'artist_s_name'

lr1 <- glm(IsPopular ~ . -Obs -streams -track_name -artist_s_name, data = Train, family = "binomial", maxit = 100)
lr1 <- glm(IsPopular ~ . - Obs - streams - track_name - artist_s_name, 
           data = Train, family = "binomial", maxit = 100)

lr1.step <-step(lr1, direction = "backward")
summary(lr1.step)

initial_vars <- names(Train)
final_vars <- names(lr1.step$model)
removed_vars <- setdiff(initial_vars, final_vars)

names(Train)[names(Train) == 'acousticness_%'] <- 'acousticness'
names(Spotify)[names(Spotify) == 'acousticness_%'] <- 'acousticness'
names(Validation)[names(Validation) == 'acousticness_%'] <- 'acousticness'

lr2 <- glm(IsPopular ~ in_spotify_playlists + in_deezer_playlists + in_shazam_charts + released_year + in_spotify_charts + in_deezer_charts + acousticness, 
           data = Train,
           family = "binomial",
           maxit = 100)

## Association Rules

transactions <- read_csv("Downloads/Popular_Spotify_Songs (1).csv")
transactions$in_spotify_playlists <- factor(transactions$in_spotify_playlists)
transactions$in_deezer_playlists <- factor(transactions$in_deezer_playlists)
transactions$in_shazam_charts <- factor(transactions$in_shazam_charts)
transactions$released_year <- factor(transactions$released_year)
transactions$in_spotify_charts <- factor(transactions$in_spotify_charts)
transactions$in_deezer_charts <- factor(transactions$in_deezer_charts)
transactions$`acousticness_%` <- as.numeric(transactions$`acousticness_%`)
transactions$IsPopular <- ifelse(transactions$streams >= AvgStreams, 1, 0)

subset_transactions <- transactions[, c("in_spotify_playlists", "in_deezer_playlists", "in_shazam_charts", "released_year", "in_spotify_charts", "in_deezer_charts", "acousticness_%", "IsPopular")]

subset_transactions <- subset_transactions[, c("in_spotify_playlists", "in_deezer_playlists", "in_shazam_charts", "released_year", "in_spotify_charts", "in_deezer_charts", "acousticness_%", "IsPopular")]

subset_transactions$IsPopular <- factor(subset_transactions$IsPopular)
transactions <- as(subset_transactions, "transactions")
inspect(transactions)
rules <- apriori(transactions, parameter = list(support = 0.1, confidence = 0.5))
inspect(rules)




transactions <- as(subset_transactions, "transactions")

Transactions <-as.matrix(subset_transactions, "transactions")
fp.trans <- as(Transactions, "transactions")

inspect(Transactions)
rules <- apriori(subset_transactions, parameter = list(support = 0.1, confidence = 0.5, target = "rules"), 
                 appearance = list(rhs= MyList))




Validation$track_name <- factor(Validation$track_name, levels = levels(Train$track_name))
Train$track_name <- factor(Train$track_name, levels = union(levels(Train$track_name), levels(Validation$track_name)))


predTrain <- predict(lr1, type="response")
predVal <- predict(lr1, type="response",newdata = Validation)
confusionMatrix(as.factor(ifelse(predVal > 0.5, 1, 0)), as.factor(Validation$IsPopular), positive ="1")





