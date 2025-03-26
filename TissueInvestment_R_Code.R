#libraries for packages used
library(emmeans)
library(lme4)
library(dplyr) #formatting
library(gt) #creating supplementary tables

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#import data in
df1 <- read.csv("~/Data/TissueInvestmentData.csv") 
df1$Morph <- gsub(" ", "", df1$Morph)
df1$Species <- gsub(". ", "oecilia_", df1$Species)

#brain residuals
lm_brain <- lm(log10(Brain.Weight..mg.) ~ log10(df1$Length..mm.), data = df1)
summary(lm_brain)
br_in <- coef(lm_brain)[1]
br_s <- coef(lm_brain)[2]
df1$ResBL <-log10(df1$Brain.Weight..mg.) - (log10(df1$Length..mm.)*br_s+br_in)

#testis residuals
lm_gonad <- lm(log10(Gonad.Weight..mg.) ~ log10(df1$Length..mm.), data = df1)
summary(lm_gonad)
g_in <- coef(lm_gonad)[1]
g_s <- coef(lm_gonad)[2]
df1$ResGL <-log10(df1$Gonad.Weight..mg.) - (log10(df1$Length..mm.)*g_s+g_in)

#subgroups
nfdf <- df1[df1$Sex != "Female",] #no females
nmdf <- df1[df1$Sex != "Male",]#no males

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ANOVA Comparisons
#Male Gonad
#raw gonad
##all males
raw_g <- aov(nfdf$Gonad.Weight..mg.~nfdf$Species+nfdf$Species:nfdf$Morph) #ANOVA
summary(raw_g)
raw_g_comp <-pairs(emmeans(raw_g, "Morph", adjust = "tukey")) #pairwise comparisons
summary(raw_g_comp)

#residual gonad weight (males only)
##all males
res_g <-aov(nfdf$ResGL~nfdf$Species+nfdf$Species:nfdf$Morph) #ANOVA
summary(res_g)
res_g_comp <-pairs(emmeans(res_g, "Morph", adjust = "tukey")) #pairwise comparisons
summary(res_g_comp)

#-------
#Brain

#raw brain weight (males and females)
raw_br_all <- aov(df1$Brain.Weight..mg.~df1$Sex*df1$Species+df1$Species:df1$Morph) #ANOVA
summary(raw_br_all)

##just males
raw_br_m <- aov(nfdf$Brain.Weight..mg.~nfdf$Species+nfdf$Species:nfdf$Morph) #ANOVA
summary(raw_br_m)
raw_br_m_comp <-pairs(emmeans(raw_br_m, "Morph", adjust = "tukey")) #pairwise comparisons
summary(raw_br_m_comp)

##just females
raw_br_f <- aov(nmdf$Brain.Weight..mg.~nmdf$Species) #ANOVA
summary(raw_br_f)
raw_br_f_comp <-pairs(emmeans(aov(nmdf$Brain.Weight..mg.~nmdf$Species), "Species", adjust = "tukey")) #pairwise comparisons
summary(raw_br_f_comp)


#residual brain weight (males and females)
res_br_all <- aov(df1$ResBL~df1$Sex*df1$Species+df1$Species:df1$Morph) #ANOVA
summary(res_br_all)
##just males
res_br_m <- aov(nfdf$ResBL~nfdf$Species+nfdf$Species:nfdf$Morph) #ANOVA
summary(res_br_m)
res_br_m_comp <- pairs(emmeans(res_br_m, "Morph", adjust = "tukey")) #pairwise comparisons
summary(res_br_m_comp)

##just females
res_br_f <-aov(nmdf$ResBL~nmdf$Species) #ANOVA
summary(res_br_f)
res_br_f_comp <-pairs(emmeans(aov(nmdf$ResBL~nmdf$Species), "Species", adjust = "tukey")) #pairwise comparisons
summary(res_br_f_comp)

#-------

#Neuron/glia ratio
ngr_all <- aov(df1$Neuron.Ratio~df1$Sex*df1$Species+df1$Species:df1$Morph) #ANOVA
summary(ngr_all)
##just males
ngr_m <-aov(nfdf$Neuron.Ratio~nfdf$Species+nfdf$Species:nfdf$Morph) #ANOVA
summary(ngr_m)
ngr_sp_comp <- pairs(emmeans(ngr_m, "Species", adjust = "tukey")) #pairwise comparisons with P. parae morphs combined
summary(ngr_sp_comp)
ngr_m_comp <-pairs(emmeans(ngr_m, "Morph", adjust = "tukey")) #pairwise comparisons separated morphs of P. parae
summary(ngr_m_comp)

##just females
ngr_f <-aov(nmdf$Neuron.Ratio~nmdf$Species) #ANOVA
summary(ngr_f)
ngr_f_comp <- pairs(emmeans(aov(nmdf$Neuron.Ratio~nmdf$Species), "Species", adjust = "tukey")) #pairwise comparisons
summary(ngr_f_comp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Linear models comparing the relationship between gonad size and brain traits across males (either combining or separating morphs)

#residual brain by residual gonad across males
##by species (P. pare morphs combined)
bvg <- lm(ResBL ~ ResGL*Species,data=nfdf)
summary(bvg)
bvg_s0 <- test(emtrends(bvg, ~Species, var="ResGL"), adjust="none") #slope to 0
bvg_s_comp <- emtrends(bvg, pairwise~Species, var="ResGL") #slope significant from each other
bvg_i_comp <- emmeans(bvg, pairwise~Species, adjust = "tukey") #intercepts different from each other


## by morph (P. pare morphs separated)
bvg_m <- lm(ResBL ~ ResGL*Morph,data=nfdf)
summary(bvg_m)
bvg_m_s0 <- test(emtrends(bvg_m, ~Morph, var="ResGL"), adjust="none") #slope to 0
bvg_m_s_comp <- emtrends(bvg_m, pairwise~Morph, var="ResGL", adjust="tukey") #slope significant from each other
bvg_m_i_comp <- emmeans(bvg_m, pairwise~Morph, adjust = "tukey") #intercepts different from each other


#-------

#neuron/glia ratio by residual gonad across males
#Species (P. pare morphs combined)
ngrvg<- lm(Neuron.Ratio ~ ResGL*Species,data=nfdf)
summary(ngrvg)
ngrvg_s0 <-test(emtrends(ngrvg, ~Species, var="ResGL")) #slope to 0
ngrvg_s_comp <-emtrends(ngrvg, pairwise~Species, var="ResGL") #slope significant from each other
ngrvg_i_comp <-emmeans(ngrvg, pairwise~Species, adjust = "tukey") #intercepts different from each other

#morph (P. pare morphs separated)
ngrvg_m<- lm(Neuron.Ratio ~ ResGL*Morph,data=nfdf)
summary(ngrvg_m)
ngrvg_m_s0 <-test(emtrends(ngrvg_m, ~Morph, var="ResGL")) #slope to 0
ngrvg_m_s_comp <-emtrends(ngrvg_m, pairwise~Morph, var="ResGL") #slope significant from each other
ngrvg_m_i_comp <-emmeans(ngrvg_m, pairwise~Morph, adjust = "tukey") #intercepts different from each other


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Create supplementary tables
#functions
#create a table of Effects, degrees of freedom, Sum square, Mean square, F value, and P value pased on anova object input and title
create_anova_table <- function(aov_df, title) {
  # Summarize the ANOVA output
  df <- summary(aov_df)[[1]]

  # Add effect names as a column
  Effects <- rownames(df)
  df <- cbind(Effects, data.frame(df, row.names = NULL))

  # Rename the columns
  colnames(df) <- c("Effect", "df", "Sum Sq", "Mean Sq", "F value", "Pr (>F)")

  # Format the dataframe and convert the columns to numeric
  df <- format.data.frame(df, digits = 3)
  df[, 2:6] <- lapply(df[, 2:6], as.numeric)

  # Create the gt table
  df %>%
    gt() %>%
    tab_header(title = title) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        columns = vars(`Pr (>F)`),
        rows = `Pr (>F)` < 0.05
      )
    )
}
#create a table of comparison results inclduing the comparison, estimate, standard error, degrees of freedom, t ratio, and p value when given pairwise comparison of emmeans
create_comp_table <- function(comp, title) {
  # Create a dataframe from the summary of the pairwise comparison
  comp_df <- as.data.frame(cbind(
    summary(comp)$contrast,
    summary(comp)$estimate,
    summary(comp)$SE,
    summary(comp)$df,
    summary(comp)$t.ratio,
    summary(comp)$p.value
  ))

  # Convert columns 2 to 6 to numeric
  comp_df[, 2:6] <- lapply(comp_df[, 2:6], as.numeric)

  # Rename the columns
  colnames(comp_df) <- c("Contrast", "Estimate", "SE", "df", "t ratio", "P value")

  # Format the data to show 3 digits
  comp_df <- format.data.frame(comp_df, digits = 3)
  comp_df[, 2:6] <- lapply(comp_df[, 2:6], as.numeric)

  # Create the gt table and apply styling
  comp_df %>%
    gt() %>%
    tab_header(title = title) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        columns = vars(`P value`),
        rows = `P value` < 0.05
      )
    )
}
#create a table inclduing a linear model and its comparisons of slope to 0, pairwise comparisons of slope, and pairwise comparisons of intercept
create_regression_table <- function(lm_model, emtrends_s0, emtrends_comp, emmeans_comp, title) {
  lm_summary <- summary(lm_model)$coefficients
  lm_df <- data.frame(Effect = rownames(lm_summary), lm_summary)
  rownames(lm_df) <- NULL

  # Rename the columns for clarity
  colnames(lm_df) <- c("Effect", "Estimate", "Std. Error", "t value", "P value")
  # Format the lm_df with 3 digits
  lm_df <- format.data.frame(lm_df, digits = 3)

  # Get the emtrends (Slope to 0 and Pairwise Comparisons)
  emtrends_s0_df <- data.frame(emtrends_s0)
  emtrends_s0_df <- subset(emtrends_s0_df, select =-df)
  colnames(emtrends_s0_df) <- c("Effect", "Estimate", "Std. Error", "t value", "P value")
  emtrends_s0_df[, 2:5] <- lapply(emtrends_s0_df[, 2:5], as.numeric)
  emtrends_s0_df <- format.data.frame(emtrends_s0_df, digits = 3)

  #slope comparisons
  emtrends_comp_df <- as.data.frame(summary(emtrends_comp)$contrast)
  emtrends_comp_df <- subset(emtrends_comp_df, select =-df)
  colnames(emtrends_comp_df) <- c("Effect", "Estimate", "Std. Error", "t value", "P value")
  emtrends_comp_df[, 2:5] <- lapply(emtrends_comp_df[, 2:5], as.numeric)
  emtrends_comp_df <- format.data.frame(emtrends_comp_df, digits = 3)

  # Get the emmeans (Intercepts Pairwise Comparisons)
 emmeans_comp_df <- as.data.frame(summary(emmeans_comp)$contrast)
  emmeans_comp_df <- subset(emmeans_comp_df, select =-df)
  colnames(emmeans_comp_df) <- c("Effect", "Estimate", "Std. Error", "t value", "P value")
  emmeans_comp_df[, 2:5] <- lapply(emmeans_comp_df[, 2:5], as.numeric)
  emmeans_comp_df <- format.data.frame(emmeans_comp_df, digits = 3)

  all_dfs <- list(lm_df, emtrends_s0_df, emtrends_comp_df, emmeans_comp_df)
  results_df <- do.call(rbind, all_dfs)
  # Combine all results into a single data frame for display
  title_reg <- data.frame(Effect = "Linear Model", Estimate = NA, `Std. Error` = NA, `t value` = NA, `P value` = NA)
  colnames(title_reg) <- c("Effect", "Estimate", "Std. Error", "t value", "P value")
  title_slope <- data.frame(Effect = "Slope to 0", Estimate = NA, `Std. Error` = NA, `t value` = NA, `P value` = NA)
    colnames(title_slope) <- c("Effect", "Estimate", "Std. Error", "t value", "P value")
  title_slope_comp <- data.frame(Effect = "Pairwise Slope Comparisons", Estimate = NA, `Std. Error` = NA, `t value` = NA, `P value` = NA)
  colnames(title_slope_comp) <- c("Effect", "Estimate", "Std. Error", "t value", "P value")
title_inter_comp <- data.frame(Effect = "Pairwise Intercept Comparisons", Estimate = NA, `Std. Error` = NA, `t value` = NA, `P value` = NA)
  colnames(title_inter_comp) <- c("Effect", "Estimate", "Std. Error", "t value", "P value")


  results_df <- rbind(title_reg, lm_df, title_slope, emtrends_s0_df , title_slope_comp, emtrends_comp_df, title_inter_comp, emmeans_comp_df)

  # Format the results table using `gt`
  results_df %>%
    gt() %>%
    tab_header(title = title) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        columns = vars(`P value`),
        rows = `P value` < 0.05
      )
    ) %>%
    tab_spanner(
      label = "Model Results",
      columns = c("Effect", "Estimate", "Std. Error", "t value", "P value")
    )
}

#-------

#raw gonads
create_anova_table(raw_g, "Raw Gonads Across Male Morph and Species ANOVA")
create_comp_table(raw_g_comp, "Male Raw Gonad Pairwise Comparisons")

#residual gonad weight (males only)
create_anova_table(res_g, "Residual Gonads Across Male Morph and Species ANOVA")
create_comp_table(res_g_comp, "Male Residual Gonad Pairwise Comparisons")

#raw brain weight
create_anova_table(raw_br_all, "Raw Brains Across Sex and Species ANOVA")
##just males
create_anova_table(raw_br_m, "Raw Brains Across Male Morph and Species ANOVA")
create_comp_table(raw_br_m_comp, "Male Raw Brain Pairwise Comparisons")
##just females
create_anova_table(raw_br_f, "raw brain female")
create_comp_table(raw_br_f_comp, "raw brain female comp")

#residual brain weight
create_anova_table(res_br_all, "Residual Brains Across Sex and Species ANOVA")
##just males
create_anova_table(res_br_m, "Residual Brains Across Male Morph and Species ANOVA")
create_comp_table(res_br_m_comp, "Male Residual Brain Pairwise Comparisons")
##just females
create_anova_table(res_br_f, "res brain female")
create_comp_table(res_br_f_comp, "Female Residual Brain Pairwise Comparisons")


#Neuron/glia ratio
create_anova_table(ngr_all, "Neuron/Glia Ratio Across Sex and Species ANOVA")
##just males
create_anova_table(ngr_m, "Neuron/Glia Ratio Across Male Morph and Species ANOVA")
create_comp_table(ngr_m_comp, "Male Neuron/Glia Ratio Pairwise Comparisons")
##just females
create_anova_table(ngr_f, "Neuron/Glia Ratio Across Female Species ANOVA")
create_comp_table(ngr_f_comp, "Female Neuron/Glia Ratio Pairwise Comparisons")

#-------

#male brain x gonad
## species (P. pare morphs combined)
create_regression_table(bvg, bvg_s0, bvg_s_comp, bvg_i_comp, "Male Brain Weight by Gonad Weight Species")
## morphs (P. pare morphs separated)
create_regression_table(bvg_m, bvg_m_s0, bvg_m_s_comp, bvg_m_i_comp, "Male Brain Weight by Gonad Weight Morph")

#male neuron/glia ratio x gonad
## species (P. pare morphs combined)
create_regression_table(ngrvg, ngrvg_s0, ngrvg_s_comp, ngrvg_i_comp, "Male Neuron/Glia Ratio by Gonad Weight Species")
## morphs(P. pare morphs separated)
create_regression_table(ngrvg_m, ngrvg_m_s0, ngrvg_m_s_comp, ngrvg_m_i_comp, "Male Neuron/Glia Ratio by Gonad Weight Morph")

