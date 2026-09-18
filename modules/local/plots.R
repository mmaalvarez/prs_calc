library(tidyverse)
library(pROC)
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")

if(interactive()){setwd("/g/strcombio/fsupek_data/users/malvarez/projects/lucia/prs/external/hartwig/plots/")}


prs_model = "young_hg37"

predictions_table = read_tsv(paste0("../calculate_PRS/res/predictions_", prs_model))


## Stratify individuals into deciles based on their PRS
## Calculate the Odds Ratio (OR) for each decile compared to a reference group (usually the middle decile or the lowest decile)
## Plot the ORs across deciles

## Also do the ROC AUC


### OR deciles

reference_decile = 1

OR_predictions_table = predictions_table %>% 
  # Assign deciles to PRS
  mutate(decile = ntile(prs, 10)) %>% 
  # Count cases and controls per decile
  group_by(decile) %>%
  summarise(cases    = sum(type == "Lung_tumor"),
            controls = sum(type == "Control"))
  
# Get reference counts
ref = OR_predictions_table %>% 
  filter(decile == reference_decile)

# Calculate OR and 95% CI relative to reference
OR_predictions_table = OR_predictions_table %>% 
  mutate(or    = (cases / controls) / (ref$cases / ref$controls),
         se    = sqrt(1/cases + 1/controls + 1/ref$cases + 1/ref$controls),
         lower = exp(log(or) - 1.96 * se),
         upper = exp(log(or) + 1.96 * se))

OR_deciles_plot = ggplot(OR_predictions_table, 
                         aes(x = factor(decile),
                             y = or)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "PRS decile",
       y = paste0("Odds Ratio: (cases / controls) / (decile ", reference_decile, "'s cases / decile ", reference_decile, "'s controls)"),
       title = paste0("PRS model from ", gsub("_hg.*", " et al.", str_to_title(prs_model)))) +
  theme_bw() +
  theme(text = element_text(size = 15))
ggsave(filename = paste0("OR_deciles_plot_", prs_model, ".jpg"),
       plot = OR_deciles_plot,
       width = 12.5,
       height = 7,
       dpi = 300,
       bg = "white")


### ROC AUC

roc = roc(response = predictions_table$y, predictor = predictions_table$prs)

auc_val = as.numeric(auc(roc))

roc_plot = ggroc(roc) +
  geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1), 
               linetype = "dashed", color = "grey50") +
  ggtitle(sprintf("PRS model from %s - Hartwig's lung vs ovary+uterus+prostate+melanoma\nPRS ROC (AUC=%.3f)", 
                  gsub("_hg.*", " et al.", str_to_title(prs_model)), 
                  auc_val)) +
  theme_bw()

ggsave(filename = paste0("roc_", prs_model, ".jpg"),
       plot = roc_plot,
       width = 12.5,
       height = 7,
       dpi = 300,
       bg = "white")
