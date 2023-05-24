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



model_df <- read.csv('G:/Quesnel_2022/Modelling/core_trees.csv')


plot_relationship <- function(df, xvar, yvar, groupvar=NULL, xlab=NULL, ylab=NULL) {
  # Create scatter plot with line of fit for each group in groupvar
  xlab <- ifelse(is.null(xlab), xvar, xlab)
  ylab <- ifelse(is.null(ylab), yvar, ylab)

  if(is.null(groupvar)){

    p <- ggplot(df, aes_string(x=xvar, y=yvar)) +
      geom_point(size = 3) +
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
  } else {
    p <- ggplot(df, aes_string(x=xvar, y=yvar, color=groupvar)) +
      geom_point(size = 3) +
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
  }

  return(p)
}

p1 <- plot_relationship(model_df,
                        'log(vol_concave)',
                        "log(sum_bai_5)",
                        xlab = 'Crown Volume (Log)',
                        ylab = 'BAI (5 Year) (Log)')

p2 <- plot_relationship(model_df,
                  'log(vol_concave)',
                  "log(sum_bai_5)",
                  "Species",
                  xlab = 'Crown Volume (Log)',
                  ylab = 'BAI (5 Year) (Log)')


p1 + p2


model_df <- model_df %>% filter()


# Model 1: simple no plot effect - BAI vs volume

model_1 <- lm(log(sum_bai_5) ~ log(vol_concave), data = model_df)

# Model 2: bai vs vol with plot random effect

model_2 <- lmerTest::lmer(log(sum_bai_5) ~ log(vol_concave) + (1|PlotID), data = model_df)

# Model 3: simple with species included as a fixed effect

model_3 <- lm(log(sum_bai_5) ~ log(vol_concave) + Species, data = model_df)

# Model 4: species included as a fixed effect plus random effect of plotID

model_4 <- lmerTest::lmer(log(sum_bai_5) ~ log(vol_concave) + Species + (1|PlotID), data = model_df)

# Model 5: simple no plot effect competition vs BAI

model_5 <- lm(log(sum_bai_5) ~ cindex, data = model_df)

# Model 6: Competition plus plot effect

model_6 <- lmerTest::lmer(log(sum_bai_5) ~ cindex + (1|PlotID), data = model_df)

# Model 7: Competiton plus species

model_7 <- lm(log(sum_bai_5) ~ cindex + Species, data = model_df)

# Model 8: Competiton plus species plus random effect of plot

model_8 <- lmerTest::lmer(log(sum_bai_5) ~ cindex + Species + (1|PlotID), data = model_df)

# Model 9: Competition plus Volume

model_9 <- lm(log(sum_bai_5) ~ cindex + log(vol_concave), data = model_df)

# Model 10: Competition plus volume with random effect of plot

model_10 <- lmerTest::lmer(log(sum_bai_5) ~ cindex + log(vol_concave) + (1|PlotID), data = model_df)

# Model 11: Competition plus volume with fixed effect of species

model_11 <- lm(log(sum_bai_5) ~ cindex + log(vol_concave) + Species, data = model_df)

# Model 12: Competition plus volume with fixed effect of species + random effect of plot

model_12 <- lmerTest::lmer(log(sum_bai_5) ~ cindex + log(vol_concave) + Species + (1|PlotID), data = model_df)

# Model 13: Solar Irradiance

model_13 <- lm(log(sum_bai_5) ~ irr_mean, data = model_df)

# Model 14: Solar Irradiance plus volume

model_14 <- lm(log(sum_bai_5) ~ irr_mean + log(vol_concave), data = model_df)

# Model 15: Solar Irradiance plus volume + random effect of plot

model_15 <- lmerTest::lmer(log(sum_bai_5) ~ irr_mean + log(vol_concave) + (1|PlotID), data = model_df)

# Model 16: Solar Irradiance plus volume + random effect of plot + plus competition

model_16 <- lmerTest::lmer(log(sum_bai_5) ~ irr_mean + log(vol_concave) + cindex + (1|PlotID), data = model_df)

models <- list(
  "Model 1"  = model_1,
  "Model 2"  = model_2,
  "Model 3"  = model_3,
  "Model 4"  = model_4,
  "Model 5"  = model_5,
  "Model 6"  = model_6,
  "Model 7"  = model_7,
  "Model 8"  = model_8,
  "Model 9"  = model_9,
  "Model 10" = model_10,
  "Model 11" = model_11,
  "Model 12" = model_12,
  "Model 13" = model_13,
  "Model 14" = model_14,
  "Model 15" = model_15,
  "Model 16" = model_16
)

library(modelsummary)

modelsummary(models, fmt = fmt_decimal(1,3), statistic = c('std.error', "p.value"))

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

# Conv to dataframe

dd = dredge(mod_all) %>%
  tibble::rownames_to_column()

lmerTest::lmer(log(sum_bai_5) ~ log(vol_concave) + Zmax + (1|PlotID), data = model_df) %>% summ()

summ(model)

model <- lm(log(sum_bai_5) ~ cindex, data = model_df)

check_model(model)

# Model 1 - Crown Volume



model_1 <- lmerTest::lmer(log(sum_bai_5) ~ afree + aindex + cindex + (1|PlotID), data = model_df)

summ(model_1)
report(model_1)
check_model(model_1)
# Model 2 - Competition Indicies

model_2 <- lmerTest::lmer(log(sum_bai_5) ~ log(cindex) + (1|PlotID), data = model_df)

summ(model_2)


model_df %>% ggplot(aes(x = cindex, y = vol_concave)) + geom_point()




summary_stats <- function(vector) {
  stats <- list(
    mean = mean(vector),
    sd = sd(vector),
    min = min(vector),
    max = max(vector)
  )
  return(stats)
}


x <- model_df$sum_bai_10


summary_stats(x)
