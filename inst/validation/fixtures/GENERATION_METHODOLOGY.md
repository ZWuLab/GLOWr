# Test Data Generation Methodology for GLOWr Validation

**File Log (reverse chronological order):**
- 2025-10-19: Created by Claude Code - Comprehensive documentation of test data generation methodology

## Related Code Files

This document describes test data generation for GLOWr validation. The code is located in:

**Data Generation Functions:**
- `tests/testthat/helper-validation.R` (lines 393-501): `generate_test_data()` function
- `tests/testthat/generate_fixtures.R`: Wrapper script that generates all 6 fixtures

**Generated Data:**
- `inst/validation/fixtures/*.rds`: The 6 test data fixture files

**Documentation:**
- `inst/validation/fixtures/README.md`: Quick reference guide
- `inst/validation/fixtures/GENERATION_METHODOLOGY.md`: This file (comprehensive methodology)

All paths are relative to the GLOWr package root.

---

**Purpose:** This document provides detailed mathematical and computational specifications for the test data fixtures used in GLOWr validation. It ensures full reproducibility, understanding, and record-keeping for all validation datasets.

---

## Table of Contents

1. [Overview](#overview)
2. [Statistical Framework](#statistical-framework)
3. [Implementation Details](#implementation-details)
4. [Fixture Specifications](#fixture-specifications)
5. [Reproducibility Protocol](#reproducibility-protocol)
6. [Quality Checks](#quality-checks)

---

## Overview

### Purpose of Test Fixtures

The test data fixtures serve three critical purposes in the GLOWr validation framework:

1. **Algorithm Validation**: Compare GLOWr implementations against legacy GLOW code with known ground truth
2. **Numerical Precision**: Verify numerical agreement within tolerance (< 1e-10)
3. **Edge Case Testing**: Cover diverse genetic architectures (LD structure, rare variants, binary/continuous traits)

### Design Principles

All test data are simulated with:
- **Known ground truth**: True causal variants, effect sizes, and causal probabilities are recorded
- **Controlled parameters**: MAF, LD structure, sample size, and genetic architecture are precisely specified
- **Reproducibility**: Fixed random seeds ensure identical data regeneration
- **Realism**: Data mimic real whole-genome sequencing studies (MAF distributions, LD patterns)

---

## Statistical Framework

### Data Generating Model

#### Notation

- **n**: Sample size (number of individuals)
- **p**: Number of genetic variants in the region
- **k**: Number of covariates (excluding intercept)
- **G**: Genotype matrix (n $\times$ p) with entries in {0, 1, 2} (additive genetic model)
- **X**: Covariate matrix (n $\times$ (k+1)) including intercept
- **Y**: Outcome vector (length n), continuous or binary
- **$\beta_G$**: Genetic effect vector (length p), with $\beta_G[j] = 0$ for null variants
- **$\beta_X$**: Covariate effect vector (length k+1)

#### Continuous Outcome Model

For continuous outcomes, we use a linear model:

```
Y = X beta_X + G beta_G + epsilon
```

where:
- **$\varepsilon \sim N(0, \sigma^2)$**: Independent error terms with $\sigma^2 = 1$
- **$\beta_X \sim N(0, 0.3^2)$**: Covariate effects (including intercept)
- **$\beta_G[j]$**:
  - = 0 for null variants (with probability $1 - \pi_{causal}$)
  - $\sim N(\mu_{effect}, (\mu_{effect}/2)^2)$ for causal variants

#### Binary Outcome Model

For binary outcomes, we use a logistic regression model:

```
logit(P(Y_i = 1)) = X_i beta_X + G_i beta_G
```

where:
- **P(Y = 1)**: Disease probability
- The intercept $\beta_X[1]$ is calibrated to achieve target prevalence
- Calibration: $\beta_X[1] \leftarrow \beta_X[1] + \log(\rho/(1-\rho)) - \text{mean}(\log(\hat{p}/(1-\hat{p})))$
  - $\rho$: Target prevalence
  - $\hat{p}$: Predicted probabilities before calibration

### Genotype Generation

#### Independent Variants (ld_structure = FALSE)

For independent variants:

```
G[i,j] ~ Binomial(2, maf[j])
```

where **maf[j]** is the minor allele frequency for variant j.

**MAF Distribution:**
- **Common variants** (rare = FALSE): maf ~ Uniform(0.05, 0.5)
- **Rare variants** (rare = TRUE): maf ~ Uniform(0.001, 0.01)

#### Correlated Variants (ld_structure = TRUE)

For variants in linkage disequilibrium, we use a Gaussian copula approach:

**Step 1: Generate latent multivariate normal**
```
Z ~ MVN(0, Sigma)
```

where $\Sigma$ is a p $\times$ p correlation matrix with AR(1) structure:
```
Sigma[i,j] = rho^|i-j|
```
with $\rho$ = ld_strength (default 0.3 for moderate LD, 0.2 for weak LD)

**Step 2: Transform to uniform**
```
U[i,j] = Phi(Z[i,j])
```
where $\Phi$ is the standard normal CDF

**Step 3: Transform to genotypes**

Using Hardy-Weinberg equilibrium thresholds:
```
G[i,j] = 2  if U[i,j] < maf[j]^2                              (AA)
         1  if maf[j]^2 <= U[i,j] < maf[j]^2 + 2*maf[j]*(1-maf[j])  (Aa)
         0  otherwise                                         (aa)
```

**Expected LD:**
For AR(1) structure with $\rho = 0.3$:
- Adjacent variants (|i-j| = 1): r $\approx$ 0.3
- Variants 2 apart (|i-j| = 2): r $\approx$ 0.09
- Variants 5 apart (|i-j| = 5): r $\approx$ 0.002

### Causal Variant Selection

**Causal fraction:** By default, $\pi_{causal}$ = n_causal / p

For each variant j:
```
true_PI[j] = 1  if j ∈ causal_idx (selected at random)
             0  otherwise
```

**Effect size distribution:**
```
true_B[j] = 0                                   if true_PI[j] = 0
            ~ N(mu_effect, (mu_effect/2)^2)    if true_PI[j] = 1
```

Typical effect sizes:
- Common variants: $\mu_{effect}$ = 0.4 - 0.5 (explains ~5-10% variance per variant)
- Rare variants: $\mu_{effect}$ = 0.7 - 0.8 (larger effects, consistent with rare variant association)

---

## Implementation Details

### Code Location

**Primary Functions:**
1. `generate_test_data()`: Main data generation function
   - Location: `tests/testthat/helper-validation.R` (lines 393-501)
   - Inputs: n, p, k, binary, maf, rare, ld_structure, ld_strength, n_causal, effect_size, seed, prevalence
   - Outputs: List with G, X, Y, true_B, true_PI, causal_idx, maf, description

2. `generate_fixtures.R`: Wrapper script for all fixtures
   - Location: `tests/testthat/generate_fixtures.R`
   - Calls `generate_test_data()` with specific parameters for each fixture
   - Saves results as .rds files in `inst/validation/fixtures/`

### Key Implementation Choices

#### 1. Genotype Encoding
- **0**: Homozygous major (aa)
- **1**: Heterozygous (Aa)
- **2**: Homozygous minor (AA)

This follows the standard **additive genetic model** where each copy of the minor allele contributes equally.

#### 2. Covariate Generation
```r
X <- cbind(1, matrix(rnorm(n * k), n, k))
colnames(X) <- c("Intercept", paste0("Cov", 1:k))
```

Covariates are **independent N(0,1)** random variables, representing:
- Age, sex, principal components in real GWAS
- Purposely uncorrelated with genotypes (conservative for Type I error)

#### 3. Binary Outcome Calibration

The intercept calibration ensures target prevalence:

```r
# Initial probabilities
prob_init <- 1 / (1 + exp(-(X %*% beta_X + G %*% beta_G)))

# Calibrate intercept
beta_X[1] <- beta_X[1] + log(rho/(1-rho)) - mean(log(prob_init/(1-prob_init)))

# Final probabilities
prob_final <- 1 / (1 + exp(-(X %*% beta_X + G %*% beta_G)))

# Generate binary outcome
Y <- rbinom(n, 1, prob_final)
```

**Note:** Actual prevalence may differ slightly from target due to binomial sampling.

#### 4. LD Structure via AR(1)

The AR(1) correlation structure:
```r
rho_mat <- matrix(0, p, p)
for (i in 1:p) {
  for (j in 1:p) {
    rho_mat[i, j] <- ld_strength^abs(i - j)
  }
}
```

This creates **exponentially decaying LD** with distance, mimicking real genomic LD patterns where nearby variants are more correlated.

---

## Fixture Specifications

### 1. test_data_simple.rds

**Purpose:** Basic sanity check with minimal complexity

**Parameters:**
```r
n = 100
p = 10
k = 2
binary = FALSE
rare = FALSE
ld_structure = FALSE
n_causal = 2
effect_size = 0.5
seed = 12345
```

**Characteristics:**
- Continuous outcome
- Independent SNPs (no LD)
- Common variants (MAF: 0.05-0.5)
- 2 causal variants out of 10 (20% causal)
- Moderate effect sizes

**Use Case:**
- Quick validation tests
- Debugging
- Verify basic functionality without confounding from LD or rare variants

**Expected Behavior:**
- Burden, SKAT, Fisher tests should all detect association
- Power depends on which variants are causal (randomly selected)

---

### 2. test_data_correlated.rds

**Purpose:** Test handling of LD structure

**Parameters:**
```r
n = 200
p = 20
k = 2
binary = FALSE
rare = FALSE
ld_structure = TRUE
ld_strength = 0.3
n_causal = 4
effect_size = 0.4
seed = 23456
```

**Characteristics:**
- Continuous outcome
- AR(1) LD structure with $\rho = 0.3$
- Common variants
- 4 causal variants out of 20 (20% causal)

**LD Pattern:**
- Adjacent variants: E[$r^2$] $\approx$ 0.09 (r $\approx$ 0.3)
- 2 SNPs apart: E[$r^2$] $\approx$ 0.008 (r $\approx$ 0.09)
- Gradually decaying correlation

**Use Case:**
- Validate correlation matrix M_Z calculation in `getZ_marg_score()`
- Test optimal weight calculation under LD
- Verify GFisher handles correlated test statistics correctly

**Expected Behavior:**
- Tests should account for correlation in M
- Ignoring correlation (using identity matrix) would inflate Type I error

---

### 3. test_data_binary.rds

**Purpose:** Binary outcome with balanced case-control ratio

**Parameters:**
```r
n = 500
p = 15
k = 2
binary = TRUE
rare = FALSE
ld_structure = FALSE
n_causal = 3
effect_size = 0.6
prevalence = 0.5  # Balanced
seed = 34567
```

**Characteristics:**
- Binary outcome (case-control)
- Target prevalence = 50% (balanced design)
- Independent SNPs
- Common variants
- 3 causal variants out of 15

**Actual Prevalence:** ~50% $\pm$ 2% (binomial variability)

**Use Case:**
- Test binary trait handling in `getZ_marg_score()`
- Validate score calculation: score = G' (Y - Ŷ) for logistic regression
- Compare with SPA-corrected scores (though not needed for balanced design)

**Expected Behavior:**
- Standard score test should perform well (balanced design, common variants)
- SPA correction should give nearly identical results

---

### 4. test_data_continuous.rds

**Purpose:** Test with additional covariates (k=3 instead of 2)

**Parameters:**
```r
n = 300
p = 12
k = 3  # More covariates
binary = FALSE
rare = FALSE
ld_structure = FALSE
n_causal = 2
effect_size = 0.5
seed = 45678
```

**Characteristics:**
- Continuous outcome
- 3 covariates (plus intercept) instead of usual 2
- Independent SNPs
- Common variants

**Use Case:**
- Verify proper handling of covariate matrix dimension
- Test residualization: G̃ = G - X(X'X)^(-1)X'G
- Ensure X matrix dimension flexibility

**Expected Behavior:**
- Results should be robust to number of covariates
- More covariates = more residualization, potentially lower power

---

### 5. test_data_rare.rds

**Purpose:** Rare variant association testing

**Parameters:**
```r
n = 1000
p = 10
k = 2
binary = TRUE
rare = TRUE
ld_structure = FALSE
n_causal = 2
effect_size = 0.8  # Larger effects
prevalence = 0.3
seed = 56789
```

**Characteristics:**
- Binary outcome (case-control)
- **Rare variants:** MAF < 0.01 (0.001 to 0.01)
- Larger sample size (n=1000) for adequate power
- Larger effect sizes (consistent with rare variant hypothesis)
- Prevalence = 30% (unbalanced, more realistic)

**MAF Distribution:** Mean MAF $\approx$ 0.005

**Use Case:**
- **Critical test for SPA correction** in `getZ_marg_score_binary_SPA()`
- Test Burden test performance (rare variants favor Burden)
- Validate MAC (minor allele count) filters

**Expected Behavior:**
- SPA correction should provide more accurate p-values than standard score test
- Burden test should outperform SKAT (rare variants likely same direction)
- Some variants may have MAC < 10, triggering SPA

---

### 6. test_data_sparse.rds

**Purpose:** Sparse causal signal (needle in haystack)

**Parameters:**
```r
n = 400
p = 50
k = 2
binary = FALSE
rare = FALSE
ld_structure = TRUE
ld_strength = 0.2  # Weak LD
n_causal = 2  # Only 2 out of 50!
effect_size = 0.7
seed = 67890
```

**Characteristics:**
- Continuous outcome
- **Sparse causal structure:** Only 4% causal (2/50)
- Weak LD structure ($\rho = 0.2$)
- Common variants
- Larger effect sizes to maintain power

**Use Case:**
- Test optimal weights with sparse true_PI
- Validate BE (Best Estimator) vs APE (Asymptotically Powerful Estimator)
- Challenge for SKAT (many null variants dilute signal)

**Expected Behavior:**
- Optimal weights should downweight null variants
- BE weights (using sparse PI) should outperform equal weights
- APE weights should provide robust performance

---

## Reproducibility Protocol

### Generating All Fixtures

**Method 1: From package root**
```bash
cd /path/to/GLOWr
Rscript tests/testthat/generate_fixtures.R
```

**Method 2: From tests/testthat directory**
```bash
cd /path/to/GLOWr/tests/testthat
Rscript generate_fixtures.R
```

**Expected Output:**
```
Generating test data fixtures...
Fixture directory: <GLOWr>/inst/validation/fixtures

[1/6] Generating test_data_simple.rds...
  n=100, p=10, k=2, continuous outcome, common variants (MAF: 0.051-0.488), independent, 2 causal
  Saved to: test_data_simple.rds

[2/6] Generating test_data_correlated.rds...
  ...

[6/6] Generating test_data_sparse.rds...
  ...

======================================================================
Fixture generation completed!
======================================================================

Generated fixtures:
  - test_data_simple.rds (10.5 KB)
  - test_data_correlated.rds (14.2 KB)
  - test_data_binary.rds (12.8 KB)
  - test_data_continuous.rds (11.7 KB)
  - test_data_rare.rds (22.1 KB)
  - test_data_sparse.rds (18.6 KB)

Fixtures saved to: .../inst/validation/fixtures
Total size: 89.9 KB
```

### Verifying Reproducibility

To verify that fixtures are identical to original:

```r
# Load original fixture
original <- readRDS("inst/validation/fixtures/test_data_simple.rds")

# Regenerate with same seed
regenerated <- generate_test_data(
  n = 100, p = 10, k = 2, binary = FALSE, rare = FALSE,
  ld_structure = FALSE, n_causal = 2, effect_size = 0.5,
  seed = 12345
)

# Verify identity
all.equal(original$G, regenerated$G)           # Should be TRUE
all.equal(original$X, regenerated$X)           # Should be TRUE
all.equal(original$Y, regenerated$Y)           # Should be TRUE
all.equal(original$true_B, regenerated$true_B) # Should be TRUE
```

### Random Seeds

| Fixture                     | Seed  | Purpose |
|-----------------------------|-------|---------|
| test_data_simple.rds        | 12345 | Base seed |
| test_data_correlated.rds    | 23456 | Sequential |
| test_data_binary.rds        | 34567 | Sequential |
| test_data_continuous.rds    | 45678 | Sequential |
| test_data_rare.rds          | 56789 | Sequential |
| test_data_sparse.rds        | 67890 | Sequential |

**Rationale:** Sequential seeds (12345, 23456, ...) are easy to remember and ensure no overlap in random number streams.

---

## Quality Checks

### Post-Generation Validation

After generating fixtures, verify:

#### 1. MAF Distribution

```r
data <- readRDS("inst/validation/fixtures/test_data_simple.rds")

# Check MAF range
range(data$maf)  # Should be in [0.05, 0.5] for common variants

# Check no monomorphic variants (MAF > 0)
all(data$maf > 0)  # Should be TRUE
```

#### 2. Genotype Encoding

```r
# Check values are in {0, 1, 2}
all(data$G %in% c(0, 1, 2))  # Should be TRUE

# Check column MAFs match specified
observed_maf <- colMeans(data$G) / 2
cor(observed_maf, data$maf)  # Should be > 0.99 for large n
```

#### 3. Outcome Distribution

**For continuous:**
```r
# Should be approximately normal (if n large)
shapiro.test(data$Y)  # p-value > 0.05 expected
```

**For binary:**
```r
# Check prevalence
mean(data$Y)  # Should be close to target prevalence +/- binomial SE

# For test_data_binary (prevalence = 0.5, n = 500):
# Expected: 0.5 +/- 1.96 * sqrt(0.5*0.5/500) = 0.5 +/- 0.044
```

#### 4. Causal Variant Consistency

```r
# true_PI should match causal_idx
which(data$true_PI == 1)  # Should equal data$causal_idx
sum(data$true_PI)         # Should equal n_causal
```

#### 5. LD Structure (for correlated fixtures)

```r
data <- readRDS("inst/validation/fixtures/test_data_correlated.rds")

# Compute observed LD
observed_cor <- cor(data$G)

# For AR(1) with rho=0.3, check diagonal elements
observed_cor[1, 2]  # Should be approx 0.3 +/- sampling error
observed_cor[1, 3]  # Should be approx 0.09 +/- sampling error
```

---

## Maintenance and Updates

### When to Regenerate Fixtures

Regenerate fixtures if:
1. **Bug found in `generate_test_data()`**: Fix bug, regenerate, re-validate all functions
2. **Need new test scenario**: Add new fixture with different parameters
3. **R version update**: Ensure consistency across R versions (RNG may differ)

### Version Control

- **Fixtures (.rds files)**: Tracked in git (small file size, ~90 KB total)
- **Generation scripts**: Tracked in git (tests/testthat/generate_fixtures.R)
- **Documentation**: This file, tracked in git

### Adding New Fixtures

To add a new fixture:

1. **Edit `generate_fixtures.R`:**
```r
# Add new fixture
test_data_new <- generate_test_data(
  n = ...,
  p = ...,
  # ... other parameters ...
  seed = XXXXX  # Choose unique seed
)
saveRDS(test_data_new, file.path(fixture_dir, "test_data_new.rds"))
```

2. **Update `inst/validation/fixtures/README.md`:**
   - Add description of new fixture
   - Document use case

3. **Update this file (GENERATION_METHODOLOGY.md):**
   - Add detailed specification
   - Explain rationale for new fixture

4. **Regenerate and commit:**
```bash
Rscript tests/testthat/generate_fixtures.R
git add inst/validation/fixtures/test_data_new.rds
git add tests/testthat/generate_fixtures.R
git add inst/validation/fixtures/README.md
git add inst/validation/fixtures/GENERATION_METHODOLOGY.md
git commit -m "Add new test fixture: test_data_new.rds"
```

---

## References

### Statistical Methods

1. **Gaussian Copula for LD:**
   - Cario, M. C., & Nelson, B. L. (1997). Modeling and generating random vectors with arbitrary marginal distributions and correlation matrix. Technical Report, Department of Industrial Engineering and Management Sciences, Northwestern University.

2. **Logistic Regression Outcome Generation:**
   - Austin, P. C., et al. (2021). Generating survival times to simulate Cox proportional hazards models with time-varying covariates. Statistics in Medicine.

3. **Hardy-Weinberg Equilibrium:**
   - Wigginton, J. E., Cutler, D. J., & Abecasis, G. R. (2005). A note on exact tests of Hardy-Weinberg equilibrium. The American Journal of Human Genetics, 76(5), 887-893.

### Code Dependencies

- **mvtnorm**: Multivariate normal distribution for correlated genotypes
- **stats**: Standard R statistical functions (rnorm, rbinom, pnorm, qnorm)

---

## Contact and Questions

For questions about test data generation:
1. Review this document and `tests/testthat/helper-validation.R`
2. Check fixture generation script: `tests/testthat/generate_fixtures.R`
3. Contact: Zheyang Wu (project PI)

---

**Document Version:** 1.0
**Last Updated:** 2025-10-19
**Author:** GLOWr Development Team
