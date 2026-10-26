# Stability Analysis of Canonical Correlation Analysis (CCA)

> ** Read the Full Report:** You can view the complete research paper [here (PDF)](STA496_FinalPaper.pdf).

An end-to-end generative simulation framework in **R** for evaluating the statistical stability, parameter recovery, and empirical power of **Canonical Correlation Analysis (CCA)** across both univariate ($k = 1$) and multivariate ($k > 1$) canonical settings.

Building upon the benchmark methodology of [Helmer et al. (2024)](https://doi.org/10.1038/s42003-024-05869-4), this project introduces a power-law covariance generative model with orthonormal QR-factorized weight structures and trace feasibility constraints to systematically evaluate how **samples-per-feature ratios ($n/p$)**, **between-set correlation magnitude ($\rho$)**, **feature dimensionality ($p_x, p_y$)**, and **dimensionality asymmetry** govern the reliability of canonical estimates.

---

## Overview & Methodology

For Canonical Correlation Analysis to yield reproducible insights, estimated canonical correlations, weights, scores, and loadings must remain stable across independent samples drawn from the same population. This repository implements a controlled Monte Carlo and Wishart simulation pipeline to quantify parameter recovery under diverse data-generating regimes:

1. **Generative Population Modeling:**
   * Within-set covariance matrices $\Sigma_X \in \mathbb{R}^{p \times p}$ and $\Sigma_Y \in \mathbb{R}^{q \times q}$ are constructed with power-law eigenvalue decay: $\Sigma_X = c_x \cdot \text{diag}(1^{-\alpha_x}, \dots, p^{-\alpha_x})$.
   * Population weight matrices $W_X = \Sigma_X^{-1/2} U_k$ and $W_Y = \Sigma_Y^{-1/2} V_k$ are generated via QR decomposition of standard normal matrices, filtered through a signal-concentration feasibility constraint:
     $$\text{tr}(\tilde{W}_X^T \Sigma_X \tilde{W}_X) \ge \frac{k \cdot \text{tr}(\Sigma_X)}{2p}$$
   * Synthetic observations $(x_i, y_i)^T \sim \mathcal{N}_{p+q}(0, \Sigma)$ and sample covariance realizations (via multivariate normal and Wishart draws) are generated across logarithmic grids of sample-to-feature ratios ($n/p \in [3, 1000]$).

2. **Component-Wise Stability & Recovery Metrics:**
   Evaluated individually across each canonical pair $i \in \{1, \dots, k\}$:
   * **Relative Association Error ($\Delta \rho_i$):** Discrepancy between empirical and population canonical correlations, $\frac{\hat{\rho}_i - \rho_i}{\rho_i}$.
   * **Weight Error ($\Delta W_i$):** Sign-invariant maximum complement of absolute cosine similarity across sets $s \in \{X, Y\}$, $\max_s (1 - \vert{}\text{cossim}(\hat{W}_{s,i}, W_{s,i})\vert{})$.
   * **Out-of-Sample Score & Loading Errors ($\Delta T_i, \Delta L_i$):** Evaluated on independent holdout test datasets $(X^{\text{test}}, Y^{\text{test}})$ to prevent in-sample overfitting bias.
   * **Inter-Sample Weight Stability ($s_i$):** Average pairwise absolute cosine similarity $\vert{}\text{cossim}(V_{1,i}, V_{2,i})\vert{}$ across repeated independent draws.
   * **Empirical Statistical Power:** Permutation-based hypothesis testing ($M = 100$ row permutations of $Y$) with Bonferroni family-wise error correction ($\alpha / k$) in the multicanonical setting.

---

## Key Empirical Findings

* **Samples-per-Feature ($n/p$) & Signal Strength Hierarchy:** While $n/p$ is the primary driver of asymptotic convergence, true correlation magnitude ($\rho_{\text{true}}$) dictates the required sample threshold. In univariate CCA, moderate correlations ($\rho = 0.5$) achieve near-perfect power and low weight error by $n/p \in [32, 64]$, whereas weak signals ($\rho = 0.1$) exhibit substantial instability even at $n/p = 300$.
* **Error Propagation in Multicanonical CCA ($k = 3$):** Higher-order canonical components stabilize significantly slower than first canonical pairs, even when holding true correlation magnitude constant ($\rho = 0.5$). Because successive canonical directions are estimated conditionally under orthogonality constraints ($a_i^T \hat{\Sigma}_X a_j = 0$), estimation variance in leading components propagates into subsequent pairs.
* **Feature Dimensionality & Asymmetry Penalty:** While balanced increases in feature dimensionality ($p_x = p_y \in \{4, 8, 16, 32\}$) have a secondary impact in univariate CCA, they induce substantial spread in multicanonical estimation. Furthermore, **dimensionality asymmetry** ($p_x \neq p_y$ with $p_x + p_y = 64$) severely degrades stability across both univariate and multivariate settings, requiring significantly higher sample densities to match balanced configurations.
* **Absence of Systematic PC1 Bias:** Evaluations of relative Principal Component 1 (PC1) similarity across $p_x = p_y = 100$ confirm that CCA weight vectors do not exhibit strong systematic bias toward the leading within-set variance axis, contrasting with behavior typically observed in Partial Least Squares (PLS).

---

## Repository Structure

```text
├── Reproduce_Results.Rmd   # Master R Markdown pipeline to simulate data, compute metrics, and generate all figures
└── README.md               # Project documentation and methodology overview
```

---

## Reproducing the Simulation & Figures

To run the simulations and regenerate all figures from the paper (`Figures 1–10`):

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/phamali/](https://github.com/phamali/)<your-repo-name>.git
   cd <your-repo-name>
   ```

2. **Knit or execute `Reproduce_Results.Rmd`:**
   Open `Reproduce_Results.Rmd` in RStudio and click **Knit**, or run it directly from the terminal:
   ```bash
   Rscript -e "rmarkdown::render('Reproduce_Results.Rmd')"
   ```

*(Optional: If you share the library calls at the top of your `Reproduce_Results.Rmd`—such as `MASS`, `ggplot2`, `dplyr`, or `parallel`—I can add the exact `install.packages(...)` one-liner right above step 2!)*

---

## Authors, Supervision & Reference

* **Authors:** Alisha Pham and Jianyang Xiao
* **Faculty Supervisor:** Professor Elena Tuzhilina (Department of Statistical Sciences, University of Toronto)
* **Primary Reference:** Helmer, M., Warrington, S., Mohammadi-Nejad, A.-R., et al. (2024). *On the stability of canonical correlation analysis and partial least squares with application to brain-behavior associations.* Communications Biology, 7(1):217.
