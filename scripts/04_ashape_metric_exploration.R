# Read in the data
df <- read.csv('H:/Quesnel_2022/process/CT1/output/crowns/ashapes/ULS22_CT1_chunk_5cmvoxel_ashapes.csv')
x <- df
#x <- df %>% select(-treeID)
# Get unique column names
unique_colnames <- unique(names(x))

# Remove duplicated rows of column names
x <- x[, unique_colnames]

# Select all character columns
char_cols <- sapply(x, is.character)

# Convert selected columns to numeric
x <- mutate_if(x, char_cols, as.numeric) %>% filter(!is.na(vol_convex))



library(purrr)
library(ggplot2)
library(tidyverse)
library(ggplot2)

# Create a list to store the histograms
histograms <- list()

# Loop through each column in the data frame and create a histogram
library(ggplot2)

x <- x %>% select(-c(treeID,X,Y))

# Histograms of each column


histograms <- map(names(x), ~ ggplot(x, aes(x = .data[[.x]])) +
                    geom_histogram(binwidth = 0.5, color = "black", fill = "white") +
                    labs(x = .x)) + theme_light()


scatter_plots <- list()



for(i in 2:ncol(x)) {

  for(j in 1:(i-1)) {
    scatter_plots[[i]] <-
      ggplot(x, aes(x = .data[[names(x)[j]]], y = .data[[names(x)[i]]])) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = "lm", se = FALSE) +
      labs(x = names(x)[j], y = names(x)[i]) +
      geom_text(x = min(x[,j])+0.1*diff(range(x[,j])), y = max(x[,i])-0.1*diff(range(x[,i])),
                label = paste0("Corr = ", round(cor(x[,j], x[,i]), 2))) + theme_light()
  }
}
library(gridExtra)
# Combine plots into a grid
grid.arrange(histograms[[1]], histograms[[2]], histograms[[3]], histograms[[4]],
             histograms[[5]], histograms[[6]], histograms[[7]], histograms[[8]],
             scatter_plots[[7]], ncol = 3, nrow = 3, top = "Distribution of Metrics and Correlations")


# Display the histograms in a grid
gridExtra::grid.arrange(grobs = histograms, ncol = 3)
