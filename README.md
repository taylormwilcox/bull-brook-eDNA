# bull-brook-eDNA
Scripts for Wilcox et al. (2018) in Ecosphere (DOI: 10.1002/ecs2.2500)

Taylor Wilcox
18 September 2018

All scripts use the dataset "FlintRock_UClarkFork_Bull_Brook.csv" which is composed of 630 sites across both HUCs,
except for BullBrook_electro_eDNA_concordance, which draws on "UpperClarkFork_eDNA-electrofishing-concordance.csv"

BullBrook_data_summaries
* Summary stats for covariates 
* Correlation matrix of covariates and bull trout eDNA presence/absence
* Various other summaries in Results of ms (# positive sites, flow in positive sites, etc)

BullBrook_model_selection
* Hierarchical model selection (RE then FE) for bull trout models
* Includes user-defined functions to estimate AUC and confusion matrices (model accuracy)

BullBrook_influence_and_CI
* Fitting of patch influence object
* Calculating/plotting Cook's Distance/patch
* Parametric bootstrap of fixed effect confidence intervals

BullBrook_kernal_density_delta
* Script to calculate and plot density of sites by TEMP and FLOW

BullBrook_temp_flow_response_curves
* Figures of Pr(bull trout) ~ FLOW and TEMP

BullBrook_brk_response_curves_fig
* Figure showing modeled effected of brook trout relative abundance on probability of bull trout

BullBrook_random_patch_response_cuves
* Figure showing response curves by patch with random intercepts and slopes

BullBrook_brook_eDNA_quants
* Figure showing relationship between brook trout eDNA concentrations and habitat covariates

BullBrook_BRK_projections
* Brook trout eDNA concentration model and assessment
* Brook trout future simulations
* Bull trout climate projections
* Associated figures and summary statistics

BullBrook_electro_eDNA_concordance
* concordance between historic electrofishing and eDNA sampling

Except for Figure 1, all figures should be completely reproduceable from the R scripts. Figure 1 was constructed from a composite of ArcGIS maps and the underlying data for this visualization is contained within the dataset (locations and detection results).
