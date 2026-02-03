### Estimate Wealth Effect of Stock Crash on COnsumption

# --- 1. Load Required Libraries ---
library(readxl)    # Excel import
library(rugarch)   # GARCH model for shocks
library(seasonal)  # X-13 Seasonal Adjustment
library(urca)      # Cointegration tests
library(tseries)   # ADF tests
library(dynlm)     # ECM estimation
library(vars)      # VECM and IRF
library(lmtest)    # For Durbin-Watson test
library(tidyverse)     # Data wrangling
library(zoo)       # Date handling

# --- 2. Data Import ---
file_path <- "data_WE.xlsx" # Ensure this file is in your working directory
daily_raw <- read_excel(file_path, sheet = "Daily")
monthly_raw <- read_excel(file_path, sheet = "Monthly")

# --- 3. First Step: Daily Stock Shocks (GARCH) ---
daily_raw <- daily_raw %>%
  arrange(Date) %>%
  mutate(returns = c(NA, diff(log(ihsg)))) %>%
  filter(!is.na(returns))

# Estimate GARCH(1,1)
garch_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0))
)
garch_fit <- ugarchfit(spec = garch_spec, data = daily_raw$returns)
daily_raw$shock <- as.numeric(residuals(garch_fit, standardize = TRUE))

# Aggregate to Monthly Cumulative Shock
monthly_shocks <- daily_raw %>%
  mutate(YearMonth = format(as.Date(Date), "%Y-%m")) %>%
  group_by(YearMonth) %>%
  summarise(stock_shock = sum(shock, na.rm = TRUE))

# Merge with Monthly data
monthly_raw$YearMonth <- format(as.Date(monthly_raw$Date), "%Y-%m")
final_df <- inner_join(monthly_raw, monthly_shocks, by = "YearMonth")

# --- 4. Seasonal Adjustment & TS Conversion ---
start_p <- c(2019, 1)

# Convert to TS objects
rsi_raw <- ts(log(final_df$rsi), start = start_p, frequency = 12)
inc_raw <- ts(log(final_df$Income), start = start_p, frequency = 12)
cpi_raw <- ts(log(final_df$ihk), start = start_p, frequency = 12)
infl <- ts(final_df$infl, start = start_p, frequency = 12)
ir_ts <- ts(final_df$ir, start = start_p, frequency = 12) 
shock_ts <- ts(final_df$stock_shock, start = start_p, frequency = 12)
data_raw <- cbind(rsi_raw, inc_raw, ir_ts, infl, shock_ts)
plot(data_raw)

# Plot shock_ts
library(ggplot2)
library(scales)

# Convert ts object to data frame for ggplot
shock_df <- data.frame(
  Date = as.Date(time(shock_ts)),
  Shock = as.numeric(shock_ts)
)

# Create the plot
p_shock <- ggplot(shock_df, aes(x = Date, y = Shock)) +
  geom_col(fill = "#A0522D", alpha = 0.8) + # Terracotta color to match theme
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  # Add a threshold line for "Panic" (1.5 SD)
  geom_hline(yintercept = -1.5 * sd(shock_df$Shock), linetype = "dashed", color = "red") +
  annotate("text", x = min(shock_df$Date), y = -1.5 * sd(shock_df$Shock), 
           label = "Panic Threshold (1.5 SD=6.28)", vjust = 1.5, hjust = 0.2, size = 3, color = "red") +
  labs(title = "Monthly Cumulative Stock Market Shocks (GARCH-derived)",
       subtitle = "Aggregated daily standardized residuals from IHSG",
       y = "Cumulative Shock Value (%)", x = "Year") +
  theme_minimal() +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")

print(p_shock)
#ggsave("stock_shock_plot.png", width = 8, height = 4)

# Seasonal Adjustment (X-13ARIMA-SEATS) for RSI, Income, and CPI
rsi_sa <- final(seas(rsi_raw))
inc_sa <- final(seas(inc_raw))
cpi_sa <- final(seas(cpi_raw))

# Create a combined data matrix for testing
data_matrix <- cbind(rsi_raw, inc_raw, ir_ts, cpi_sa)
colnames(data_matrix) <- c("RSI", "Income", "IR", "CPI")
plot(data_matrix)

# --- 5. Stationarity & Cointegration ---
cat("\n--- ADF Test on Seasonally Adjusted Series (First Difference) ---\n")
print(apply(data_matrix, 2, function(x) adf.test(diff(x))$p.value))

# Johansen Test (Trace Statistic)
jo_test <- ca.jo(data_matrix, type = "trace", ecdet = "const", K = 2)
summary(jo_test)

# --- 6. Step Two: Error Correction Model (ECM) ---
# Long-Run Equilibrium
long_run <- lm(rsi_raw ~ inc_raw + shock_ts + ir_ts + cpi_sa)
summary(long_run)
ect <- ts(residuals(long_run), start = start_p, frequency = 12)

# Short-Run Dynamics (The ECM)
ecm_model <- dynlm(d(log(rsi_raw)) ~ L(ect, 1) + d(log(inc_raw)) + + shock_ts + L(shock_ts,1) + L(shock_ts,2) + ir_ts)

# --- 7. Diagnostics & Half-Life ---
cat("\n--- ECM Model Summary ---\n")
print(summary(ecm_model))

# Durbin-Watson Test for Autocorrelation
dw_result <- dwtest(ecm_model)
cat("\n--- Durbin-Watson Test ---\n")
print(dw_result)

# Speed of Adjustment & Half-Life
rho <- coef(ecm_model)["L(ect, 1)"]
half_life <- log(0.5) / log(1 + rho)
cat("\nSpeed of Adjustment (ECT):", round(rho, 4))
cat("\nHalf-life of shock impact:", round(half_life, 2), "months\n")


##Non Linear Analysis

# --- 1. Create the Panic Dummy ---
# We define a "Panic" as a monthly shock that is 1.5 standard deviations below mean
threshold <- -1.5 * sd(shock_ts, na.rm = TRUE)
panic_dummy <- ifelse(shock_ts < threshold, 1, 0)

# --- 2. Create Asymmetric Shocks ---
# This separates the 'good news' from the 'bad news'
neg_shock <- ifelse(shock_ts < 0, shock_ts, 0) # Only contains negative values
pos_shock <- ifelse(shock_ts > 0, shock_ts, 0) # Only contains positive values

# --- 3. Re-estimate the ECM with these new variables ---

# Model A: Do only big crashes matter?
ecm_panic <- dynlm(d(rsi_raw) ~ L(ect, 1) + d(inc_raw) + 
                     panic_dummy + L(panic_dummy, 1) + L(panic_dummy, 2) + ir_ts)
dwtest(ecm_panic)

# Model B: Do consumers react differently to crashes vs. booms?
ecm_asymmetric <- dynlm(d(rsi_raw) ~ L(ect, 1) + d(inc_raw) + 
                          neg_shock + L(neg_shock, 2) + 
                          pos_shock + L(pos_shock, 2) + ir_ts)
dwtest(ecm_asymmetric)

# --- 4. Compare Results ---
summary(ecm_panic)
summary(ecm_asymmetric)

#Model diagnostics
# Load the library
library(FinTS)

# Run the ARCH-LM test on your ECM residuals
# 'lags = 12' is typical for monthly data
arch_test <- ArchTest(residuals(ecm_panic), lags = 12)

print(arch_test)

library(sandwich)
library(lmtest)

# Recalculate the ECM Summary using Newey-West Robust Standard Errors
# This "punishes" the t-stats for the ARCH effects you found
robust_results <- coeftest(ecm_panic, vcov = vcovHAC(ecm_panic))

cat("\n--- ECM Results with Robust Standard Errors (HAC) ---\n")
print(robust_results)

# Re-estimate with seasonal dummies
ecm_seasonal <- dynlm(d(rsi_raw) ~ L(ect, 1) + 
                        d(inc_raw) + 
                        L(panic_dummy, 2) + 
                        ir_ts + 
                        season(rsi_raw))

summary(ecm_seasonal)
dwtest(ecm_seasonal)

arch_test <- ArchTest(residuals(ecm_seasonal), lags = 12)
print(arch_test)

library(sandwich)
library(lmtest)

# Recalculate the ECM Summary using Newey-West Robust Standard Errors
# This "punishes" the t-stats for the ARCH effects you found
robust_results <- coeftest(ecm_seasonal, vcov = vcovHAC(ecm_seasonal))

cat("\n--- ECM Results with Robust Standard Errors (HAC) ---\n")
print(robust_results)


# --- 8. VECM Framework & Impulse Response Functions (IRF) ---
# Add stock_shock to the system
library(vars)

# 1. Prepare data (Variables are already logged)
data_vecm <- cbind(rsi_raw, inc_raw, ir_ts, shock_ts)
colnames(data_vecm) <- c("RSI", "Income", "IR", "Stock_Shock")

# 2. Run Johansen to get the VECM
# We add 'season = 12' to control for those monthly spikes
jo_test <- ca.jo(data_vecm, type="trace", ecdet="const", K=2, spec="longrun", season=12)
vecm_obj <- vec2var(jo_test, r = 1)

# 3. Generate the IRF for the Wealth Effect
# We look at how a Shock to the Market affects RSI over 12 months
irf_wealth <- irf(vecm_obj, impulse = "Stock_Shock", response = "RSI", 
                  n.ahead = 12, boot = TRUE, runs = 1000)

# 4. Plot the result (uncomment to plot a positive shock IRF)
#plot(irf_wealth, main = "Impulse Response: Market Crash -> Retail Sales",
#     ylab = "Response of log(RSI)", xlab = "Months After Crash")

# Generate the IRF for all variables
# Remove 'response' to get everything at once
irf_crash <- irf(vecm_obj, impulse = "Stock_Shock", 
                 n.ahead = 12, boot = TRUE, runs = 1000)

# Function to invert and scale the IRF to a -1% (0.01) shock
# By default, IRF is +1 Standard Deviation. We flip it to negative.
invert_irf <- function(irf_obj) {
  irf_obj$irf$Stock_Shock <- irf_obj$irf$Stock_Shock * -1
  irf_obj$Lower$Stock_Shock <- irf_obj$Lower$Stock_Shock * -1
  irf_obj$Upper$Stock_Shock <- irf_obj$Upper$Stock_Shock * -1
  return(irf_obj)
}

irf_negative <- invert_irf(irf_crash)

# Plotting the results
library(ggplot2)
library(gridExtra)
library(tidyr)

# Function to extract IRF data into a clean data frame
get_irf_df <- function(irf_obj, resp_name) {
  data.frame(
    Month = 0:12,
    # Scale by -1 to show a CRASH instead of a boom
    Effect = irf_obj$irf$Stock_Shock[, resp_name] * -1,
    Lower  = irf_obj$Lower$Stock_Shock[, resp_name] * -1,
    Upper  = irf_obj$Upper$Stock_Shock[, resp_name] * -1,
    Variable = resp_name
  )
}

# 1. Extract data for each variable
df_rsi    <- get_irf_df(irf_crash, "RSI")
df_income <- get_irf_df(irf_crash, "Income")
df_ir     <- get_irf_df(irf_crash, "IR")

# 2. Plotting Function for consistency
plot_irf_pretty <- function(df, title, color) {
  ggplot(df, aes(x = Month, y = Effect)) +
    geom_hline(yintercept = 0, color = "black", linetype = "solid") +
    # The Shaded Confidence Interval
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = color, alpha = 0.2) +
    # The Impulse Response Line
    geom_line(color = color, size = 1) +
    geom_point(color = color) +
    scale_x_continuous(breaks = 0:12) +
    labs(title = title, y = "Response (%)", x = "Months after Crash") +
    theme_minimal()
}

# 3. Create the two plots
p1 <- plot_irf_pretty(df_rsi, "Response of Retail Sales (RSI)", "#c0392b")
p2 <- plot_irf_pretty(df_income, "Response of Income", "#2980b9")
#p3 <- plot_irf_pretty(df_ir, "Response of Interest Rate (IR)", "#27ae60")

# 4. Arrange in a 1x3 grid
#grid.arrange(p1, p2, p3, ncol = 3)
grid.arrange(p1, p2, ncol = 2)

############### Alternative Version of drawing IRFs ######################

library(ggplot2)
library(gridExtra)
library(scales)

# 1. Calculate Standard Deviation and 8% Scaling Factor
shock_sd <- sd(shock_ts, na.rm = TRUE)
shock_sd
scaling_factor_8pct <- -8 / shock_sd

# 2. Extract and Scale function (Now multiplying by 100 for percentage points)
get_fan_data_8pct <- function(irf_obj, resp_name, scale_val) {
  # Scaling to percentage points: Log-change * 100
  percent_scale <- scale_val * 100
  
  est <- irf_obj$irf$Stock_Shock[, resp_name] * percent_scale
  
  # SE calculation scaled to percent
  se <- (irf_obj$Upper$Stock_Shock[, resp_name] - irf_obj$irf$Stock_Shock[, resp_name]) / 1.96
  se_scaled <- abs(se * percent_scale) 
  
  data.frame(
    Month = 0:12,
    Effect = est,
    up95 = est + (1.96 * se_scaled), lo95 = est - (1.96 * se_scaled),
    up90 = est + (1.645 * se_scaled), lo90 = est - (1.645 * se_scaled),
    up68 = est + (1 * se_scaled), lo68 = est - (1 * se_scaled)
  )
}

# 3. Create DataFrames
df_rsi_8 <- get_fan_data_8pct(irf_main, "RSI", scaling_factor_8pct)
df_inc_8 <- get_fan_data_8pct ## ada bug di sini
(irf_main, "Income", scaling_factor_8pct)

# 4. Plotting Function
plot_fan_8pct <- function(df, title, base_color) {
  ggplot(df, aes(x = Month, y = Effect)) +
    # Horizontal baseline
    geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
    # Vertical line at Month 2 to highlight the significant lag
    geom_vline(xintercept = 2, linetype = "dotted", color = "grey40") +
    # Fan Chart Intervals
    geom_ribbon(aes(ymin = lo95, ymax = up95), fill = base_color, alpha = 0.1) +
    geom_ribbon(aes(ymin = lo90, ymax = up90), fill = base_color, alpha = 0.15) +
    geom_ribbon(aes(ymin = lo68, ymax = up68), fill = base_color, alpha = 0.25) +
    # Main IRF Line
    geom_line(color = base_color, linewidth = 1.2) +
    # Highlight point at Month 2
    geom_point(data = subset(df, Month == 2), color = base_color, size = 3) +
    # X-axis months
    scale_x_continuous(breaks = 0:12) +
    # CLARITY UPDATE: Fixed Y-axis labels to show 2 decimal places and %
    scale_y_continuous(labels = function(x) sprintf("%.1f%%", x)) + 
    labs(
      title = title, 
      subtitle = "Dynamic response to a 8% market crash (Shaded: 68%, 90%, 95% CI)", 
      y = "Percentage Point Change from Baseline", 
      x = "Months after Stock Market Crash"
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 12, face = "bold")
    )
}

# 5. Render
p1_final <- plot_fan_8pct(df_rsi_8, "Response of Retail Sales (RSI)", "#c0392b")
p2_final <- plot_fan_8pct(df_inc_8, "Response of Household Income", "#2980b9")

grid.arrange(p1_final, p2_final, ncol = 2)

# View the exact percentage impacts for RSI
print(df_rsi_8[, c("Month", "Effect")])

# View the exact percentage impacts for Income
print(df_inc_8[, c("Month", "Effect")])
