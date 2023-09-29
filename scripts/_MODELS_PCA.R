

# Give lm with all variables; after filtering to colinearity

mod_all = lm(blups_mod,
             formula = DBH16 ~
               vol_convex +
               lag:vol_convex +
               apex_angle +
               CRR +
               mean.NDVI +
               mean.SIPI +
               mean.NDRE1 +
               mean.NIRvNDRE +
               Zq999 +
               mean.CCI,
             na.action = 'na.fail')

# Conv to dataframe

dd = dredge(mod_all) %>%
  tibble::rownames_to_column()

# Gives AIC orders models

### best performing model
# dd_complex = dd %>%
#   subset(delta == 0)
# dd_best_complex = MuMIn::get.models(dd, subset = 1)[[1]]

### simple
dd_simple = dd %>%
  subset((has('vol_convex', 'lag:vol_convex'))) %>%
  subset(df == 4)

# Best model with 3, 4, 5; compare models with similar levels of complexity
