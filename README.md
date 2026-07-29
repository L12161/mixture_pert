# Mixture LDP — Frequency Estimation under Heterogeneous Local Differential Privacy

Code accompanying the paper on **mixture local-differential-privacy (LDP) frequency
estimation**, where users hold different privacy budgets and a *mixture* of GRR
(Generalized Randomized Response) and UE/BUE (Unary/Binary Unary Encoding) is
optimized per privacy level.

The code covers three things:

1. **Synthetic experiments** — mean-squared-error (MSE) and L1 comparisons of the
   mixture mechanism against GRR, UE, and GRR baselines on synthetic (uniform single attribute and multi-attribute dataset) data.
2. **Real-data experiments** — the same comparisons driven by an empirical and theoreticl
   per-user, per-attribute privacy-budget distribution (multi-attribute dataset : UCI Adult Dataset, Single Attribute Dataset:Kosarak and Retail ).
3. **Heatmap Visualizations** — RDP-to-DP conversion analysis used for the budget-accounting heatmap. THis helps to show how the privacy budget was managed before and after deploying moments accountant. 

There is separate codes for the baseline LDP mechaanisms and the adaptive LDP mechanisms. The baaseline mechanism involves uniform (strongest privacy budget) allocation of privacy budget across all the users, and aadaptive mechanism involves comparison with ID-LDP. 

---

## Repository layout

```
.
├── setup.m                     Run once per MATLAB session; puts src/ + experiments/ on the path.
│
├── src/                        Reusable functions (called by the experiment scripts).
│   ├── generators/             Synthetic data generators.
│   │   ├── generate_uniform.m
│   │   ├── generate_gaussian.m
│   │   ├── generate_powlaw.m
│   │   ├── generate_zipf.m        (uses sample_discrete)
│   │   ├── generate_exponential.m
│   │   └── sample_discrete.m
│   ├── optimization/           Mechanism-parameter solvers.
│   │   ├── min_opt0.m             baseline / GRR objective
│   │   ├── min_opt1.m             UE objective
│   │   ├── min_opt2.m             symmetric (a = 0.5) objective
│   │   ├── min_opt3.m             mixture objective (takes domain size d)
│   │   └── privacy_parameters.m   solves for [a, b, alpha, p, q] per epsilon level
│   ├── estimation/
│   │   ├── Est_mixture.m          unbiased frequency estimator for the mixture mechanism
│   │   └── actual_MSE.m           empirical MSE / relative error
│   ├── io/
│   │   └── reading_distribution.m reads unique epsilon values + proportions from a dataset
│   └── thirdparty/
│       └── TightPlots.m           figure-layout utility (external; keep attribution header)
│
├── experiments/                Top-level scripts — these are the ones you run.
│   ├── adaptive_LDP_mechanism_comparison.m   MSE/L1 comparison for the adaptive (mixture) mechanism — synthetic or real data, set via the domain size matrix
│   ├── baseline_LDP_mechanism_comparison.m   MSE/L1 comparison for the baseline (uniform-budget) mechanism — synthetic or real data, set via the domain size matrix
│   └── plot_error_gap.m              plots a pairwise-gap CSV
│
├── notebooks/                  Python / Jupyter analyses.
│   ├── Heatmap_DP_RDP_DP_conversion.ipynb   RDP-derived vs. direct epsilon heatmap (Colab)
│   └── Mixture_Motivation_Graph_.ipynb      motivating GRR-vs-UE comparison figure
│
├── data/                       Input CSVs (not tracked). Put PubMed files here:
│                               pubmed_probabilistic_data.csv, pubmed_unique_counts.csv
│
└── results/
    └── performance_gaps/       pairwise_gap_table_n_*_d_*.csv  + generated figures
```

---

## Requirements

- **MATLAB** R2021a or newer, with the **Optimization Toolbox** (`fmincon`, used by
  the `min_opt*` solvers).
- **Python 3.10+** for the notebooks: `numpy`, `pandas`, `matplotlib`, `seaborn`.
  The heatmap notebook is written for Google Colab (`/content/...` paths); adjust
  the `FILE_*` constants if running locally.

---

## Getting started (MATLAB)

From the repository root, once per MATLAB session:

```matlab
setup        % adds src/ and experiments/ to the path
```

After that, every function (e.g. `generate_uniform`, `min_opt3`, `Est_mixture`) is
callable by name from any script, and you can run any experiment directly:

```matlab
adaptive_LDP_mechanism_comparison    % adaptive (mixture) mechanism, synthetic or real data
baseline_LDP_mechanism_comparison    % baseline (uniform-budget) mechanism, synthetic or real data
```

> The experiment scripts begin with `clear`/`close all`. That clears the workspace
> but **not** the path, so you only need to run `setup` once after starting MATLAB.

---

## Typical workflows

### A. Synthetic or real-data comparison

```matlab
setup
adaptive_LDP_mechanism_comparison     % or baseline_LDP_mechanism_comparison
```

Each script runs the same MSE/L1 comparison pipeline for either synthetic or real
data — just set the domain size matrix inside the script to switch between the two.
It generates or loads the data, optimizes the mixture parameters via the `min_opt*`
solvers, estimates frequencies with `Est_mixture`, and plots MSE/L1 against the
baselines using `TightPlots`.

### B. RDP → DP budget heatmap (Python)

Open `notebooks/Heatmap_DP_RDP_DP_conversion.ipynb`. It consumes the three CSVs
(`pubmed_probabilistic_data.csv`, `privacy_parameters_output.csv`,
`pubmed_unique_counts.csv`) and produces the difference heatmap
(direct epsilon − RDP-derived epsilon) over feature counts × Rényi orders alpha.

---

