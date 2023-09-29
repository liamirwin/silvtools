# Building and Selecting Models

# How can we predict differences in observed basal area increment pattern

# BAI - summarized as mean and sum of last five and ten years

# Tree Size -> Height, Crown Area, Crown Volume
# Species -> Douglas-fir, Lodgepole Pine, Interior Spruce
# ^ Requires prediction for wall-to-wall area
# Soil Moisture -> Topographic Wetness Index (mean)
# Solar Irradiance -> Rayshader

library(easystats)
library(tidyverse)
library(tidymodels)
library(jtools)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(ggpmisc)
library(patchwork)
library(ggrepel)
library(modelsummary)
library(gt)

model_df <- read.csv('G:/Quesnel_2022/Modelling/core_trees.csv')
match_df <- matches_df %>% filter(status == 'TP')
model_df_filt <- model_df %>%
  filter(!(tree_id.x %in% c('CT5P3-Pl040','CT5P2-Fd027', 'CT1P2-Sx044', 'CT5P2-Fd047')) & !is.na(vol_concave))
model_df <- model_df_filt


# Convert sum_bai_5 to m^2
model_df <- model_df %>%
  mutate(sum_bai_5_m2 = sum_bai_5 / 1000000,
         mean_bai_m2 = sum_bai_5_m2/5,
         ba = (pi * (Diametr/200)^2),
         growth_percent = (mean_bai_m2 / ba * 100 ),
         growth_percent_5y = sum_bai_5_m2/ba * 100)

# Calculate the annual growth percentage
model_df <- model_df %>%
  mutate(growth_percent = (sum_bai_5_m2 / (pi * (Diametr/200)^2)) * 100 / 5)










plot_relationship <- function(df, xvar, yvar, groupvar=NULL, xlab=NULL, ylab=NULL, label = F) {
  # Create scatter plot with line of fit for each group in groupvar
  xlab <- ifelse(is.null(xlab), xvar, xlab)
  ylab <- ifelse(is.null(ylab), yvar, ylab)

  if(is.null(groupvar)){

    p <- ggplot(df, aes_string(x=xvar, y=yvar)) +
      geom_point(size = 3, alpha = 0.6) +
      stat_poly_line(method='lm', se = FALSE, linetype = 'dashed') +
      stat_poly_eq((use_label(c('eq','R2')))) +
      labs(x=xlab, y=ylab) +
      theme_classic() +
      theme(legend.position="bottom",
            text = element_text(size=12), # Change global text size
            axis.title = element_text(size=14), # Change axis title text size
            legend.title = element_text(size=12), # Change legend title text size
            legend.text = element_text(size=12) # Change legend text size
      )
    if(label == T){
             p <- p +
                 geom_text_repel(aes(label=tree_id.x))
    }

  } else {
    p <- ggplot(df, aes_string(x=xvar, y=yvar, color=groupvar)) +
      geom_point(size = 3, alpha = 0.6) +
      stat_poly_line(method='lm', se = FALSE, linetype = 'dashed') +
      stat_poly_eq((use_label(c('eq','R2')))) +
      scale_color_brewer(palette = 'Set1') +
      labs(x=xlab, y=ylab, color=groupvar) +
      theme_classic() +
      theme(legend.position="bottom",
            text = element_text(size=14), # Change global text size
            axis.title = element_text(size=14), # Change axis title text size
            legend.title = element_text(size=12), # Change legend title text size
            legend.text = element_text(size=12) # Change legend text size
      )
    if(label == T){
      p <- p +
        geom_text_repel(aes(label=tree_id.x))
    }
  }

  return(p)
}

p1 <- plot_relationship(model_df,
                        'log(vol_concave)',
                        "log(sum_bai_5)",
                        xlab = 'Crown Volume (Log)',
                        ylab = 'Cumulative BAI (5 Year) (Log)', label = F)

p2 <- plot_relationship(model_df,
                  'log(vol_concave)',
                  "log(sum_bai_5)",
                  "Species",
                  xlab = 'Crown Volume (Log)',
                  ylab = 'Cumulative BAI (5 Year) (Log)')

p3 <- plot_relationship(model_df_filt,
                        'log(vol_concave)',
                        "log(sum_bai_5)",
                        xlab = 'Crown Volume (Log)',
                        ylab = 'Cumulative BAI (5 Year) (Log)', label = F)

p4 <- plot_relationship(model_df_filt,
                        'log(vol_concave)',
                        "log(sum_bai_5)",
                        "Species",
                        xlab = 'Crown Volume (Log)',
                        ylab = 'Cumulative BAI (5 Year) (Log)')

p5 <- plot_relationship(model_df,
                        'vol_concave',
                        "sum_bai_5",
                        "PlotID",
                        xlab = 'Crown Volume (Log)',
                        ylab = 'Cumulative BAI (5 Year) (Log)')





p1 + p2 + p3 + p4

model_df <- model_df_filt

model_df <- model_df %>% rename(PlotID = PlotID.x, Diameter = Diametr)


matches_df %>% filter(Height > 0 ) %>% plot_relationship(., 'Height', 'Zmax') + geom_abline(slope = 1, intercept = 0, color = "red")

# WITHOUT RANDOM EFFECTS

# Create log variables
model_df <- model_df %>%
  mutate(log_vol_concave = log(vol_concave),
         log_vol_convex = log(vol_convex),
         log_irr_mean = log(irr_mean),
         log_afree = log(afree))

# Model 1: bai vs concave vol without plot random effect
model_1 <- lm(log_sum_bai_5 ~ log_vol_concave, data = model_df)

# Model 2: bai vs convex vol without plot random effect
model_2 <- lm(log_sum_bai_5 ~ log_vol_convex, data = model_df)

# Model 3: Solar Irradiance
model_3 <- lm(log_sum_bai_5 ~ log_irr_mean, data = model_df)

# Model 4: Competition Index (Heygi)
model_4 <- lm(log_sum_bai_5 ~ cindex, data = model_df)

# Model 5: Freegrowing Area
model_5 <- lm(log_sum_bai_5 ~ log_afree, data = model_df)

# Model 6: Topographic Wetness
model_6 <- model_df %>% filter(!is.na(twi_mean)) %>% lm(log_sum_bai_5 ~ twi_mean, data = .)

# Model 7: bai vs tree height without plot random effect
model_7 <- lm(log_sum_bai_5 ~ Zmax, data = model_df)

# Model 8: Vol Concave plus Vol Convex
model_8 <- lm(log_sum_bai_5 ~ log_vol_convex + log_vol_concave, data = model_df)

# Model 9: Vol concave plus irradiance
model_9 <- lm(log_sum_bai_5 ~ log_vol_concave + log_irr_mean, data = model_df)

# Model 10: Vol Concave plus competition
model_10 <- lm(log_sum_bai_5 ~ log_vol_concave + cindex, data = model_df)

# Model 11: Vol concave plus height
model_11 <- lm(log_sum_bai_5 ~ log_vol_concave + Zmax, data = model_df)

# Model 12: Vol Concave plus free growing area
model_12 <- lm(log_sum_bai_5 ~ log_vol_concave + log_afree, data = model_df)

# Model 13:
model_13 <- lm(log_sum_bai_5 ~ log_vol_concave + log_irr_mean + log_afree + cindex, data = model_df)

# Model 14 - Field Variables
model_14 <- lm(log_sum_bai_5 ~ Diameter, data = model_df)


models <- list(
  "Model 1"  = model_1,
  "Model 2"  = model_2,
  "Model 3"  = model_3,
  "Model 4"  = model_4,
  "Model 5"  = model_5,
  "Model 6" = model_6,
  "Model 7" = model_7,
  "Model 8" = model_8,
  "Model 9" = model_9,
  "Model 10" = model_10,
  "Model 11" = model_11,
  "Model 12" = model_12,
  "Model 13" = model_13,
  "Model 14" = model_14
)

summary_lm <- modelsummary(models,
                        fmt = fmt_decimal(2,3),
                        statistic = c("p.value"),
                        output = 'gt',
                        stars = TRUE)

summary_lm <- summary_lm %>%
  tab_spanner(label = md('**Single Variable Models**'), columns = 2:8) %>%
  tab_spanner(label = md('**Multi-variable Models**'), columns = 9:14) %>%
  tab_spanner(label = md('**Field Measured Diameter**'), columns = 15) %>%
  tab_options(
    table.font.size = px(14),
    table.font.names = "Times New Roman",
    table.border.bottom.style = "none",
    table.border.top.style = 'none'
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = c('Model 1'),
      rows = c(24,25,26)
    )
  ) %>%
  tab_style(
    style = cell_fill(color = "gray", alpha = 0.3),
    locations = cells_body(
      columns = c('Model 1'),
      rows = c(24,25,26)
    )
  )

summary_lm

# WITH RANDOM EFFECTS

# Pre-create log variables
model_df <- model_df %>%
  mutate(log_sum_bai_5 = log(sum_bai_5),
         log_vol_concave = log(vol_concave),
         log_vol_convex = log(vol_convex),
         log_irr_mean = log(irr_mean),
         log_afree = log(afree))

# Model 1: bai vs concave vol with plot random effect
model_1 <- lmerTest::lmer(log_sum_bai_5 ~ log_vol_concave + (1|PlotID), data = model_df)

# Model 2: bai vs convex vol with plot random effect
model_2 <- lmerTest::lmer(log_sum_bai_5 ~ log_vol_convex + (1|PlotID), data = model_df)

# Model 3: Solar Irradiance
model_3 <- lmerTest::lmer(log_sum_bai_5 ~ log_irr_mean + (1|PlotID), data = model_df)

# Model 4: Competition Index (Heygi)
model_4 <- lmerTest::lmer(log_sum_bai_5 ~ cindex + (1|PlotID), data = model_df)

# Model 5: Freegrowing Area
model_5 <- lmerTest::lmer(log_sum_bai_5 ~ log_afree + (1|PlotID), data = model_df)

# Model 6: Topographic Wetness
model_6 <- model_df %>% filter(!is.na(twi_mean)) %>% lmerTest::lmer(log_sum_bai_5 ~ twi_mean + (1|PlotID), data = .)

# Model 7: bai vs tree height with plot random effect
model_7 <- lmerTest::lmer(log_sum_bai_5 ~ Zmax + (1|PlotID), data = model_df)

# Model 8: Vol Concave plus Vol Convex
model_8 <- lmerTest::lmer(log_sum_bai_5 ~ log_vol_convex + log_vol_concave + (1|PlotID), data = model_df)

# Model 9: Vol concave plus irradiance
model_9 <- lmerTest::lmer(log_sum_bai_5 ~ log_vol_concave + log_irr_mean + (1|PlotID), data = model_df)

# Model 10: Vol Concave plus competition
model_10 <- lmerTest::lmer(log_sum_bai_5 ~ log_vol_concave + cindex + (1|PlotID), data = model_df)

# Model 11: Vol concave plus height
model_11 <- lmerTest::lmer(log_sum_bai_5 ~ log_vol_concave + Zmax + (1|PlotID), data = model_df)

# Model 12: Vol Concave plus free growing area
model_12 <- lmerTest::lmer(log_sum_bai_5 ~ log_vol_concave + log_afree + (1|PlotID), data = model_df)

# Model 13:
model_13 <- lmerTest::lmer(log_sum_bai_5 ~ log_vol_concave + irr_mean + log_afree + cindex + (1|PlotID), data = model_df)

# Model 14 - Field Variables
model_14 <- lmerTest::lmer(log_sum_bai_5 ~ Diameter + (1|PlotID), data = model_df)


models <- list(
  "Model 1"  = model_1,
  "Model 2"  = model_2,
  "Model 3"  = model_3,
  "Model 4"  = model_4,
  "Model 5"  = model_5,
  "Model 6" = model_6,
  "Model 7" = model_7,
  "Model 8" = model_8,
  "Model 9" = model_9,
  "Model 10" = model_10,
  "Model 11" = model_11,
  "Model 12" = model_12,
  "Model 13" = model_13,
  "Model 14" = model_14
)

summary_lme <- modelsummary(models,
                        fmt = fmt_decimal(2,3),
                        statistic = c("p.value"),
                        output = 'gt',
                        stars = TRUE)


summary_lme <- summary_lme %>%
  tab_spanner(label = md('**Single Variable Models**'), columns = 2:8) %>%
  tab_spanner(label = md('**Multi-variable Models**'), columns = 9:14) %>%
  tab_spanner(label = md('**Field Measured Diameter**'), columns = 15) %>%
  tab_options(
    table.font.size = px(14),
    table.font.names = "Times New Roman",
    table.border.bottom.style = "none",
    table.border.top.style = 'none'
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = c('Model 1'),
      rows = c(25,26,27)
    )
  ) %>%
  tab_style(
    style = cell_fill(color = "gray", alpha = 0.3),
    locations = cells_body(
      columns = c('Model 1'),
      rows = c(25,26,27)
    )
  )

summary_lme

file = "D:/Proposal_2022/Thinning Paper/Figures/model_summary_table_may29.png"

gtsave(summary, filename = file)

# CROSS VALIDATION

crossv_strat =
  function (data, var, id = '.id') {
    n <- nrow(data)
    folds = data[var][[1]]
    idx <- seq_len(n)
    fold_idx <- split(idx, folds)
    fold <- function(test) {
      list(train = modelr::resample(data, setdiff(idx, test)),
           test = modelr::resample(data, test))
    }
    cols <- purrr::transpose(purrr::map(fold_idx, fold))
    tibble::as_tibble(cols)
  }

candidate_list <- list(
  'log(vol_concave)',
  'log(vol_convex)',
  'log(irr_mean)',
  'cindex',
  'log(afree)',
  'twi_mean',
  'Zmax',
  'log(vol_convex) + log(vol_concave)',
  'log(vol_concave) + log(irr_mean)',
  'log(vol_concave) + cindex',
  'log(vol_concave) + Zmax',
  'log(vol_concave) + log(afree)',
  'log(vol_concave) + irr_mean + log(afree) + cindex',
  'Diameter'
)

model_df <- model_df %>% rename(TreeID = tree_id.x)

cross_mod = model_df %>%
  crossv_strat(var = 'TreeID') %>%
  mutate(model_1 = purrr::map(train, ~lm(noquote(paste0('sum_bai_5 ~ ', candidate_list[[1]])), data=.))) %>%
  mutate(model_2 = purrr::map(train, ~lm(noquote(paste0('sum_bai_5 ~ ', candidate_list[[2]])), data=.))) %>%
  mutate(model_3 = purrr::map(train, ~lm(noquote(paste0('sum_bai_5 ~ ', candidate_list[[3]])), data=.))) %>%
  mutate(model_4 = purrr::map(train, ~lm(noquote(paste0('sum_bai_5 ~ ', candidate_list[[4]])), data=.))) %>%
  mutate(model_5 = purrr::map(train, ~lm(noquote(paste0('sum_bai_5 ~ ', candidate_list[[5]])), data=.))) %>%
  mutate(model_6 = purrr::map(train, ~lm(noquote(paste0('sum_bai_5 ~ ', candidate_list[[6]])), data=.))) %>%
  mutate(model_7 = purrr::map(train, ~lm(noquote(paste0('sum_bai_5 ~ ', candidate_list[[7]])), data=.))) %>%
  mutate(model_8 = purrr::map(train, ~lm(noquote(paste0('sum_bai_5 ~ ', candidate_list[[8]])), data=.))) %>%
  mutate(model_9 = purrr::map(train, ~lm(noquote(paste0('sum_bai_5 ~ ', candidate_list[[9]])), data=.))) %>%
  mutate(model_10 = purrr::map(train, ~lm(noquote(paste0('sum_bai_5 ~ ', candidate_list[[10]])), data=.))) %>%
  mutate(model_11 = purrr::map(train, ~lm(noquote(paste0('sum_bai_5 ~ ', candidate_list[[11]])), data=.))) %>%
  mutate(model_12 = purrr::map(train, ~lm(noquote(paste0('sum_bai_5 ~ ', candidate_list[[12]])), data=.))) %>%
  mutate(model_13 = purrr::map(train, ~lm(noquote(paste0('sum_bai_5 ~ ', candidate_list[[13]])), data=.))) %>%
  mutate(model_14 = purrr::map(train, ~lm(noquote(paste0('sum_bai_5 ~ ', candidate_list[[14]])), data=.))) %>%
  pivot_longer(
    cols = starts_with("model"), # Only pivot the model columns
    names_to = 'model_name',
    values_to = 'model'
  ) %>%
  mutate(predicted = map2(model, test, ~ augment(.x, newdata = .y))) %>%
  unnest(predicted) %>%
  dplyr::select(-train, -test, -model) %>%
  mutate(cross_v = 'TreeID')

cross_val <- cross_mod %>%
  group_by(model_name, cross_v) %>%
  mutate(R2 = yardstick::rsq_trad_vec(sum_bai_5, .fitted),
         rmse = caret::RMSE(sum_bai_5, .fitted)) %>% dplyr::filter(!is.na(.fitted)) %>%
  dplyr::select(model_name, cross_v, R2, rmse) %>%
  distinct()

df <- cross_val %>% select(model_name, R2) %>% ungroup()

# Create the table

df %>%
  pivot_wider(names_from = model_name, values_from = R2) %>%
  select(!cross_v) %>%
  rename_with(~str_replace(., "model_", "Model "), everything()) %>%
  gt() %>%
  tab_header(
    title = "Model Performance with Plot Level Cross Validation",
    subtitle = "R-squared values for each model"
  ) %>%
  tab_options(
    table.font.size = px(18),
    table.font.names = "Times New Roman",
    table.border.bottom.style = "none",
    table.border.top.style = 'none'
  ) %>%
  data_color(
    columns = everything(),
    colors = scales::col_numeric(
      palette = c("red", "green"),
      domain = c(0, 0.4)
    )
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 2
  ) %>%
  fmt_missing(
    columns = everything(),
    missing_text = "NA"
  )


# Cross Validation Attempt 2 - GPT

library(caret)
library(e1071)

library(lme4)
library(lattice)
library(Matrix)

# Install required packages
if(!require(lme4)) install.packages("lme4")
if(!require(lattice)) install.packages("lattice")
if(!require(Matrix)) install.packages("Matrix")

library(lme4)
library(lattice)
library(Matrix)


cross_validate_model <- function(formula, data, n_folds = 10){

  formula <- as.formula(formula)

  # Set seed for reproducibility
  set.seed(123)

  # Create vectors to hold the results of each fold
  rmse_results <- c()
  rsquared_results <- c()
  mae_results <- c()

  # Perform n_folds-fold cross-validation manually
  folds <- sample(1:n_folds, size = nrow(data), replace = TRUE)

  for(i in 1:n_folds){
    # Split the data into a training set and a test set
    trainSet <- data[folds != i,]
    testSet  <- data[folds == i,]

    # Train the mixed-effects model
    model <- lmer(formula, data = trainSet, REML = FALSE)

    # Test the model on the test set
    predictions <- predict(model, newdata = testSet, re.form=NA)

    # Evaluate its performance
    performance <- postResample(pred = predictions, obs = testSet[[as.character(formula[[2]])]])

    # Store the results
    rmse_results <- c(rmse_results, performance["RMSE"])
    rsquared_results <- c(rsquared_results, performance["Rsquared"])
    mae_results <- c(mae_results, performance["MAE"])
  }

  # Average the results
  avg_rmse <- mean(rmse_results)
  avg_rsquared <- mean(rsquared_results)
  avg_mae <- mean(mae_results)

  results_df <- data.frame(formula = deparse(formula), avg_rmse = avg_rmse, avg_rsquared = avg_rsquared, avg_mae = avg_mae)

  # Print the average performance
  print(c("Average RMSE" = avg_rmse, "Average Rsquared" = avg_rsquared, "Average MAE" = avg_mae))

  # Return the average performance
  return(results_df)
}

results <- cross_validate_model("log_sum_bai_5 ~ log_vol_concave + (1|PlotID)", model_df, n_folds = 10)

# Apply the cross validation function to each model
results_list <- lapply(models, function(model) {
  formula <- deparse(formula(model))
  cross_validate_model(formula, model_df, n_folds = 10)
})

# Combine the results into one data frame
results_df <- do.call(rbind, results_list) %>%
  mutate(avg_rmse_exp = exp(avg_rmse), avg_mae_exp = exp(avg_mae))

# Load the gt package
library(gt)

# Convert the data frame to a gt table
results_table <- gt(results_df)

# Print the table
print(results_table)





# Correlation Plots

library(correlation)
library(see)

corr_df <- model_df %>% select('vol_concave','vol_convex','Zmax','n_points','cindex','afree','CRR','Diameter','sum_bai_5')

results <- summary(correlation(corr_df)) + theme_classic()

plot(results, show_data = "points") +




# Dredging (Test all possible combinations)

# Give lm with all variables; after filtering to colinearity

library(MuMIn)

mod_all = lm( formula = log(sum_bai_5) ~
               vol_concave +
               log(vol_concave) +
               vol_convex +
               log(vol_convex) +
               CRR +
               n_points +
               irr_mean +
               twi_mean +
               cindex +
               Zmax,
              data = model_df,
             na.action = 'na.fail')


dd <- dredge(mod_all)
# Filter the table
dd_nv <- dd[!grepl("vol_concave|log\\(vol_concave\\)|vol_convex|log\\(vol_convex\\)", rownames(dd)), ] %>% as.data.frame() %>%
  tibble::rownames_to_column()

dd <- dd %>% as.data.frame() %>%
  tibble::rownames_to_column()

View(dd_nv)
# Filter out rows where 'vol_concave', 'vol_convex', 'log(vol_concave)', and 'log(vol_convex)' are NA

# Forward Stepwise Regression

# Create a null model with only the intercept
null_model <- lm(log_sum_bai_5 ~ 1, data = model_df)

# Define the full model with all predictors
full_model <- lm(log_sum_bai_5 ~ log_vol_concave + log_vol_convex + log_irr_mean + cindex +
                   log_afree + twi_mean + Zmax, data = model_df)

# Perform forward selection starting from null model and moving towards the full model
# Set direction = "forward"
stepwise_model <- step(null_model, scope = list(lower = null_model, upper = full_model), direction = "forward")

# Print summary of the final stepwise model
summary(stepwise_model)



















# Three Dimensional Scatter Plot

library(dplyr)

color_set <- c("Fd" = "#BF382A", "Sx" = "#0C4B8E", "Pl" = "#31a354")


p <- model_df %>%
  plot_ly(x = ~RCC_mean, y = ~GCC_mean, z = ~BCC_mean, color = ~Species, colors = color_set) %>%
  add_markers() %>%
  layout(title = list(text = "Three Dimensional Scatter Plot Between Chromatic Coordinates",
                      font = list(size = 24, family = "Arial")),
         legend = list(orientation = "h", xanchor = "center", x = 0.5, y = -0.1,
                       font = list(size = 18, family = "Arial")),
         scene = list(xaxis = list(title = 'RCC'),
                      yaxis = list(title = 'GCC'),
                      zaxis = list(title = 'BCC')))

p


orca(p, file="D:/Proposal_2022/Thinning Paper/Figures/species_cc.pdf")








library(merDeriv)
# GPT Model selection
model_df <- model_df %>% filter(!is.na(vol_concave)) %>% rename(PlotID = PlotID.x)

# Fit a linear mixed model
model <- lmer(mean_bai_5 ~ vol_concave + (1|PlotID), data = model_df)

# You can change the response variable and the predictor based on your research question
# Extract and print model parameters
parameters <- model_parameters(model)
print(parameters)

# Check model assumptions
check <- check_model(model)
print(check)

# Check model performance
performance <- model_performance(model)
print(performance)

# Estimate effects
effects <- estimate_means(model)
print(effects)

# Calculate effect sizes
effect_sizes <- effectsize(model)
print(effect_sizes)

# Generate a report
report <- report(model)
print(report)

# Visualize the estimated effects
plot(effects)

# Visualize the model predictions
predicted_values <- estimate_expectation(model) # obtain the predicted values

# Add predicted values to the original dataframe
model_df$predicted_values <- predicted_values$Predicted

# Plot the predicted vs. actual values
ggplot(model_df, aes(x = mean_bai_5, y = predicted_values)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(x = "Actual Values", y = "Predicted Values", title = "Model Predictions vs Actual Values")




analyze_model <- function(data, response, predictor, random_effect) {

  # Fit a linear mixed model
  formula <- as.formula(paste0(response, " ~ ", predictor, " + (1|", random_effect, ")"))
  model <- lme4::lmer(formula, data = data)

  # Extract and print model parameters
  parameters <- parameters::model_parameters(model)
  print(parameters)

  # Check model assumptions
  check <- performance::check_model(model)
  print(check)

  # Check model performance
  performance <- performance::model_performance(model)
  print(performance)


  # Check if model contains at least one categorical factor
  fixed_effect_vars <- names(fixef(model))
  fixed_effect_vars <- fixed_effect_vars[fixed_effect_vars != "(Intercept)"]  # remove "(Intercept)"
  if (any(sapply(data[, fixed_effect_vars], is.factor))) {
    # Estimate effects
    effects <- modelbased::estimate_means(model)
    print(effects)

    # Visualize the estimated effects
    plot(effects)
  } else {
    message("Model contains no categorical factor. Skipping estimate_means and its plot.")
  }


  # Calculate effect sizes
  effect_sizes <- effectsize::effectsize(model)
  print(effect_sizes)

  # Generate a report
  report <- report::report(model)
  print(report)

  # Estimate model predictions
  predicted_values <- modelbased::estimate_expectation(model)

  # Ensure data has the same number of rows as predicted_values
  if (nrow(predicted_values) != nrow(data)) {
    stop("Number of predicted values does not match number of rows in data frame.")
  }

  # Add predicted values to the original dataframe
  data$predicted_values <- predicted_values$Predicted

  # Plot the predicted vs. actual values
  predicted_plot <- ggplot(data, aes_string(x = response, y = "predicted_values")) +
    geom_point() +
    geom_smooth(method = "lm", color = "red") +
    theme_classic() +
    labs(x = "Actual Values", y = "Predicted Values", title = "Model Predictions vs Actual Values")

  return()
}

# Example usage
analyze_model(data = model_df, response = "sum_bai_5", predictor = "vol_concave", random_effect = "PlotID")



