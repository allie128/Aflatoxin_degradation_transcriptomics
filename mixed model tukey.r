# ============================================================
# Mixed-Effects Model (Treatment x Time, Experiment as Random)
# with Tukey-Adjusted Post-Hoc Comparisons
#
# Design: 2 experiments (random) x 2 treatments (fixed) x
#         3 time points (fixed) x 3 replicates = 36 obs
# ============================================================

library(lme4)
library(lmerTest)   # adds p-values to lmer via Satterthwaite df
library(emmeans)    # post-hoc comparisons for mixed models
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load your data
# ------------------------------------------------------------
# Expects a long-format CSV with one row per sample:
#   experiment,treatment,time,rep,response
# See sample_data.csv for the expected format (2 experiments x
# 2 treatments x 3 time points x 3 replicates = 36 rows).
mneg <- read.csv("aflatoxin_m_neg.csv")
mpos <- read.csv("aflatoxin_m_pos.csv")

mneg$experiment <- factor(mneg$experiment)
mneg$treatment  <- factor(mneg$treatment)
mneg$time       <- factor(mneg$time, levels = c("D0", "D4", "D7"))
mpos$experiment <- factor(mpos$experiment)
mpos$treatment  <- factor(mpos$treatment)
mpos$time       <- factor(mpos$time, levels = c("D0", "D4", "D7"))


str(mneg)
str(mpos)

# ------------------------------------------------------------
# 2. Mixed-effects model
# ------------------------------------------------------------
# Treatment and time (and their interaction) are fixed effects.
# Experiment is a random effect -- this treats the 2 experiments
# as a random sample from a broader population of possible
# experimental runs, and accounts for the fact that samples
# within the same experiment aren't fully independent.
model_neg <- lmer(response ~ treatment * time + (1 | experiment), data = mneg)
model_pos <- lmer(response ~ treatment * time + (1 | experiment), data = mpos)

summary(model_neg)
summary(model_pos)
anova(model_neg)   # Type III tests (Satterthwaite df) for fixed effects
anova(model_pos)

# ------------------------------------------------------------
# 3. Check assumptions
# ------------------------------------------------------------
# Normality of residuals
shapiro.test(resid(model_neg))
shapiro.test(resid(model_pos))
qqnorm(resid(model_neg)); qqline(resid(model_neg))
qqnorm(resid(model_pos)); qqline(resid(model_pos))

# Residuals vs fitted (check homogeneity of variance)
plot(fitted(model_neg), resid(model_neg),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, lty = 2)
plot(fitted(model_pos), resid(model_pos),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, lty = 2)


# Random effect (experiment) variance -- if this is ~0, the
# random effect may not be doing much and a simple two-way
# ANOVA might suffice. If it's sizable, keep the mixed model.
VarCorr(model_neg)
VarCorr(model_pos)
print(VarCorr(model_neg), comp = "Var")
print(VarCorr(model_pos), comp = "Var")
# ------------------------------------------------------------
# 4. Tukey-adjusted post-hoc comparisons (via emmeans)
# ------------------------------------------------------------
# TukeyHSD() does not work on mixed models (lmer objects), so
# emmeans is the standard replacement -- it computes estimated
# marginal means and applies the same Tukey adjustment for
# multiple comparisons.

# a) All pairwise comparisons across every treatment x time
#    combination (best choice if the interaction is significant)
emm_int_neg <- emmeans(model_neg, ~ treatment * time)
tukey_int_neg <- pairs(emm_int_neg, adjust = "tukey")
tukey_int_df_neg <- as.data.frame(tukey_int_neg)
tukey_int_df_neg <- tukey_int_df_neg[order(tukey_int_df_neg$p.value), ]
print(tukey_int_df_neg, digits = 3)
emm_int_pos <- emmeans(model_pos, ~ treatment * time)
tukey_int_pos <- pairs(emm_int_pos, adjust = "tukey")
tukey_int_df_pos <- as.data.frame(tukey_int_pos)
tukey_int_df_pos <- tukey_int_df_pos[order(tukey_int_df_pos$p.value), ]
print(tukey_int_df_pos, digits = 3)

# b) Main effect comparisons (use if interaction is NOT significant)
pairs(emmeans(model, ~ treatment), adjust = "tukey")
pairs(emmeans(model, ~ time), adjust = "tukey")

# ------------------------------------------------------------
# 5. Pull out specific comparisons of interest
#    e.g., compare time points WITHIN each treatment
#    (this is what "compare 2 time points" usually means)
# ------------------------------------------------------------
# emmeans lets you do this directly and cleanly using "by":
time_within_treatment_neg <- emmeans(model_neg, ~ time | treatment)
pairs(time_within_treatment_neg, adjust = "tukey")
time_within_treatment_pos <- emmeans(model_pos, ~ time | treatment)
pairs(time_within_treatment_pos, adjust = "tukey")


# Or treatment comparisons within each time point:
treatment_within_time_neg <- emmeans(model_neg, ~ treatment | time)
pairs(treatment_within_time_neg, adjust = "tukey")
treatment_within_time_pos <- emmeans(model_pos, ~ treatment | time)
pairs(treatment_within_time_pos, adjust = "tukey")

# ------------------------------------------------------------
# 6. Visualize
# ------------------------------------------------------------
ggplot(mneg, aes(x = time, y = response, color = treatment, group = treatment)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  facet_wrap(~experiment) +
  labs(title = "Treatment x Time interaction, by experiment",
       y = "Response", x = "Time point") +
  theme_minimal()
