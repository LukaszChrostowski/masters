library(dplyr)

# Bias, SE, RMSE for the given result
# filter out CI/variance related results
df_filtered <- dr1 %>%
  filter(!grepl("ci", variable))

result <- df_filtered %>%
  group_by(y, variable) %>%
  summarise(
    trues = mean(trues),
    mean_estimate = mean(value),
    Bias = mean(value) - mean(trues),
    SE = sqrt(sum((value - mean(value))^2) / (n() - 1)),
    RMSE = sqrt((mean(value) - mean(trues))^2 + (sqrt(sum((value - mean(value))^2) / (n() - 1)))^2)
  ) %>%
  ungroup()


result <- variance %>%
  filter(!grepl("ci", variable)) %>%
  group_by(y, variable, type) %>%
  summarise(mean_value = mean(value)) %>%
  ungroup()

# Display results
print(result, n = nrow(result))

# Coverage rate
df_ci <- variance %>%
  filter(grepl("ci", variable))

coverage_rate <- df_ci %>%
  group_by(y, type, variable) %>%
  summarise(coverage_rate = mean(value)) %>%
  ungroup()

# Display results
print(coverage_rate, n = nrow(coverage_rate))


# The same for ipw 
ipw_cloglog_new <- ipw_cloglog_new[!is.nan(ipw_cloglog_new$value), ] # may be replace by probit/logit(ipw.rds) results 
df <- ipw_cloglog_new
ipw <- df %>%
  filter(!grepl("ci", variable))

true_values <- ipw %>%
  filter(variable == "trues") %>%
  select(sample, y, true_value = value)


results <- ipw %>%
  filter(variable != "trues") %>%
  left_join(true_values, by = c("sample", "y")) %>%
  group_by(sample, y, variable) %>%
  summarise(
    Bias = mean(value) - mean(true_value), 
    SE = sd(value),                         
    RMSE = sqrt(mean((value - true_value)^2))  
  )

print(results, n = 24)
