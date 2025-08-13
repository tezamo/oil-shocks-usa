# The Macroeconomic Effects of Oil Supply News Shocks in the USA

This repository contains the code, data, and report for a research project analyzing the macroeconomic impacts of oil supply news shocks on key US economic indicators, such as the Consumer Price Index (CPI) and Industrial Production Index (IPI). The study covers the period from January 1975 to December 2022 and employs Vector Autoregression (VAR) models, including reduced-form and triangular VAR with Cholesky decomposition, to examine dynamic interrelationships and impulse responses.

The analysis draws inspiration from related work (e.g., Känzig's findings) and focuses on how informational shocks in oil supply propagate through the US economy, affecting inflation and industrial output.

## Repository Structure

- **Report**: 
  - [`report`](https://tezamo.github.io/oil-shocks-usa/report.html):  It includes the introduction, methodology (data aggregation, lag length determination, VAR estimation, stability assessments), results (optimal lag lengths, coefficient matrices, eigenvalues, impulse response functions), and conclusion. The report embeds R code snippets for visualization and analysis setup.

- **Data Files**:
  - `cpi.xlsx`: Monthly Consumer Price Index (CPI) data for the US from January 1975 to December 2022.
  - `oilSupplyNewsShocks.xlsx`: Oil supply news shock series, capturing anticipated changes and informational variances in oil supply.
  - `ipi.xlsx`: Monthly Industrial Production Index (IPI) data for the US from January 1975 to December 2022.
  - (Note: An additional oil supply surprise series is referenced in the report but may be derived or included within the oilSupplyNewsShocks file.)

- **Code**:
  - `analysis.m`: MATLAB script for performing the core analysis, including lag length selection (using AIC/BIC), estimation of reduced-form and triangular VAR models, Cholesky decomposition, eigenvalue stability checks, and generation of Impulse Response Functions (IRFs). This script loads the Excel data files and replicates the empirical results described in the report.

  ```
  For privacy reasons, the full code  are not included in this public version. Please contact me to request access to the complete project.
  ```


## Requirements

- **Software**: MATLAB (version R2020a or later recommended for compatibility with econometric toolboxes like the Econometrics Toolbox).
- **Libraries/Toolboxes**: MATLAB's Econometrics Toolbox for VAR modeling, impulse responses, and Cholesky decomposition. No additional installations are needed beyond standard MATLAB.
- **Data Handling**: Ensure Excel files are readable; MATLAB's `readtable` or `xlsread` functions can import the data.

## How to Run the Analysis

1. Clone this repository:
   ```
   git clone https://github.com/tezamo/oil-shocks-usa.git
   ```
   
2. Open MATLAB and navigate to the repository directory.

3. Load the data:
   - Use commands like `cpi_data = readtable('cpi.xlsx');` to import the Excel files.

4. Run the MATLAB script:
   ```
   analysis
   ```
   - This will:
     - Determine optimal lag lengths (AIC suggests 20 lags; BIC suggests 17).
     - Estimate the reduced-form VAR and triangular VAR models.
     - Compute coefficient matrices, residual covariances, and eigenvalues for stability.
     - Generate Impulse Response Functions (IRFs) for responses of CPI and IPI to oil supply surprises and news shocks.
     - Output results to console or save figures/plots (e.g., IRFs as described in Section 5.2 of the report).

5. View the report: Open [`report`](https://tezamo.github.io/oil-shocks-usa/report.html)  in a web browser for a detailed narrative, including embedded figures (e.g., IRFs) and methodological explanations.

## Key Findings (from the Report)

- **Optimal Lag Length**: AIC = 20 lags; BIC = 17 lags, balancing model complexity and fit.
- **VAR Results**: Negative autoregressive tendencies in oil supply news shocks; inverse relationships between CPI and IPI residuals.
- **Stability**: Eigenvalues indicate system stability with some oscillatory behavior; Matrix A (reduced-form) shows decreasing magnitudes over iterations.
- **IRFs**: Oil supply news shocks have immediate, intense effects on CPI (fluctuating then dissipating) and short-lived impacts on IPI.
- **Implications**: Highlights the need for policymakers to monitor oil supply information dynamics for broader economic stability.


## References

- The analysis references Känzig's work on oil supply shocks and macroeconomic dynamics.
- Data sourced from standard economic databases (e.g., US CPI and IPI from official sources like the Bureau of Labor Statistics or Federal Reserve).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details (if not present, assume open for academic use).

For questions or contributions, open an issue or contact the author.