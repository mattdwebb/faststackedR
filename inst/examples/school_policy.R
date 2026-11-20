library(faststackedR)
school_data <- read.csv(system.file("extdata", "school_policy.csv",
                                    package = "faststackedR"))
# This standalone API uses Inf, not zero, for never-treated schools.
school_data$first_treat[school_data$first_treat == 0] <- Inf
fit <- fit(school_data, "test_score", "year", gvar = "first_treat",
           cluster = "school", data_mode = "rcs", event_window = c(-1, 1),
           covariates = "baseline_score", inference = c("native", "jackknife"))
example_result <- data.frame(estimate = fit$estimate, cluster_se = fit$native_se,
                             jackknife_se = fit$jk_se, jackknife_p = fit$jk_p)
print(example_result, digits = 7, row.names = FALSE)
