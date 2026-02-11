# Data

Two datasets:
AVAGLIO and DFCI. Everyone received the same treatment.
Columns: age, sex, kps, mgmt, eor, os, os_status
We make the `eor` column as binary, named `isGTR`
With this change, all columns now are numeric (no strings, no factors)

# Main Method

## Input:

PLD from treatment group. Columns: age, sex, kps, mgmt, eor, os, os_status
These columns are not pre-processed (i.e., mean and variance may not be standardized)

Summary from external data predictors: One dataframe with three columns
xname, type, value
xname can be one of the predictor nmes: {age, sex, kps, mgmt, isGTR}
type: a number, typically 1 or 2.
value: mean(xname^type) in the control data

Summary from external data outcome: Unweighted RMST defined as E(min{T, tau}).
Here `tau` is a pre-defined cut-off, typically 24 months.

## Procedure:
We find the normalization constants for each predictor column.
We normalize each predictor column.
We adjust the target data frame accordingly.
We fit a Weibull regression model on PLD.
For patients with censored outcome, we generate E(min{T, tau} | X, delta, Y).
We setup and solve the linear program in terms of alpha.

## Output:
We return alpha and the weights and theta_min, theta_max.

# Required Helper Functions:
Given a dataframe with numeric predictors, fit a Weibull regression and
estimate E(min{T, tau} | X, delta, Y).
> [!NOTE]
> See impute_censored

Given vector of (Y, delta) and weights and tau, find RMST.
> [!NOTE]
> See weighted_RMST
