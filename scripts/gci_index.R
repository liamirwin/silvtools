# Changes to NDGCI

# Original Index

calculate_ndgci <- function(norm_growth, norm_competition) {
  ndgci = (norm_growth - norm_competition) / (norm_growth + norm_competition)
  return(ndgci)
}


# New Index


calculate_gci <- function(norm_growth, norm_competition){
  # New index should range from 0 to 1
  gci = ((1 - norm_growth) + norm_competition)/2
  # Tree with value of 1 - slow growing, high competition
  # Tree with value of 0 - fast growing, low competition
  # Middle range of 0.5 - both low growth/comp and high growth/comp
}

((1-bai_norm)+cindex_norm)/2


# Calculate Linear BAI Models
calculate_linear_bai_models <- function(treetops) {
  # Compute several new variables related to Basal Area Increment (BAI)
  treetops <- treetops %>% mutate(log_sum_bai_5 = (4.75962 + 0.75351 * log(vol_concave)), sum_bai_5 = exp(log_sum_bai_5),
                                  norm_cindex = ((cindex - min(cindex))/(max(cindex) - min(cindex))),
                                  norm_sum_bai_5 = ((sum_bai_5 - min(sum_bai_5))/(max(sum_bai_5) - min(sum_bai_5))),
                                  ndgci = (norm_sum_bai_5 - norm_cindex)/(norm_sum_bai_5 + norm_cindex),
                                  gci = ((1 - norm_sum_bai_5) + norm_cindex)/2)

  # Return the updated treetops data
  return(treetops)
}

ttops <- ttops %>%
    calculate_linear_bai_models() %>%
     calculate_gc_ratio()


p1 <- ggplot(ttops, aes(x = norm_cindex, y = norm_sum_bai_5, color = gci)) +
  geom_point(alpha = 0.6) +
  viridis::scale_color_viridis(name = "GCI") +
  labs(x = "Normalized Competition", y = "Normalized Growth", title = 'New Index') +
  theme_bw() +
  theme(legend.position = 'bottom')

p2 <- ggplot(ttops, aes(x = norm_cindex, y = norm_sum_bai_5, color = ndgci)) +
  geom_point(alpha = 0.6) +
  viridis::scale_color_viridis(name = "NDGCI") +
  labs(x = "Normalized Competition", y = "Normalized Growth", title = 'Old Index') +
  theme_bw() +
  theme(legend.position = 'bottom')


p1 + p2


