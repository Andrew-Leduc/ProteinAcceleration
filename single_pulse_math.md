# Single Pulse + Pseudotime — Mathematical Framework

## The Problem

With a **single metabolic pulse** of duration $T$, each cell yields two measurements per protein:

$$L_p(i), \quad M_p(i)$$

But the ODE system has **three unknowns** per protein per cell:

$$\frac{dL}{dt} = -\alpha L, \qquad \frac{dM}{dt} = \beta - \alpha M$$

$$\Rightarrow \quad L(T) = L_0\,e^{-\alpha T}, \qquad M(T) = \frac{\beta}{\alpha}(1-e^{-\alpha T})$$

**Unknowns:** $\alpha$, $\beta$, $L_0$ — **under-determined** at the single-cell level.

---

## The Solution: Pseudotime as a Temporal Axis

Consider a **population** of $N$ cells distributed across a biological trajectory (e.g., cells adapting to drug, differentiating, activating). Assign each cell a pseudotime $\tau_i \in [0,1]$ derived from UMAP + diffusion pseudotime (DPT).

**Key assumption:** A cell at pseudotime $\tau$ was at pseudotime $\tau - \delta\tau$ when the pulse began, where $\delta\tau$ corresponds to $T$ hours of actual time. Therefore:

$$L_0^{(p)}(i) \approx \underbrace{\text{Total}_p\bigl(\tau_i - \delta\tau\bigr)}_{\text{past state}}= \mathbb{E}\!\left[L_p(j) + M_p(j) \;\middle|\; \tau_j \approx \tau_i - \delta\tau\right]$$

This is estimated empirically as a Gaussian-weighted mean over cells at the target pseudotime:

$$\hat{L}_0^{(p)}(i) = \sum_j w_{ij} \cdot \bigl(L_p(j) + M_p(j)\bigr), \qquad w_{ij} \propto \exp\!\left(-\frac{(\tau_j - (\tau_i - \delta\tau))^2}{2\sigma^2}\right)$$

---

## Parameter Recovery

With $\hat{L}_0^{(p)}(i)$ in hand, the system is **exactly determined** — 2 measurements, 2 remaining unknowns ($\alpha$, $\beta$):

### Step 1 — Degradation rate from Light

$$L_p(i) = \hat{L}_0^{(p)}(i)\,e^{-\alpha T}$$

$$\boxed{\hat{\alpha}_p(i) = -\frac{\ln\!\left(L_p(i)\,/\,\hat{L}_0^{(p)}(i)\right)}{T}}$$

### Step 2 — Synthesis rate from Medium

$$M_p(i) = \frac{\beta_p}{\alpha_p}\!\left(1 - e^{-\alpha_p T}\right)$$

$$\boxed{\hat{\beta}_p(i) = \frac{\hat{\alpha}_p(i)\cdot M_p(i)}{1 - e^{-\hat{\alpha}_p(i)\,T}}}$$

### Step 3 — Steady-state abundance

$$\boxed{\hat{P}_{ss}^{(p)}(i) = \frac{\hat{\beta}_p(i)}{\hat{\alpha}_p(i)}}$$

### Step 4 — Classify dynamics

Compare the **recovered past state** $\hat{L}_0$ to the **current total** $L(T)+M(T)$:

| Comparison | Conclusion |
|---|---|
| $L_p(T) + M_p(T) > \hat{L}_0^{(p)}$ | Protein **increasing** along trajectory |
| $L_p(T) + M_p(T) = \hat{L}_0^{(p)}$ | **Steady state** |
| $L_p(T) + M_p(T) < \hat{L}_0^{(p)}$ | Protein **decreasing** along trajectory |

---

## Connection to RNA Velocity

This framework is the **proteomics analogue of RNA velocity** (La Manno et al. 2018, Bergen et al. 2020):

| RNA velocity | This method |
|---|---|
| Unspliced mRNA $u$ | Light protein $L$ (old, pre-existing) |
| Spliced mRNA $s$ | Medium protein $M$ (newly synthesised) |
| $u \to s$ dynamics encode direction | $L \to M$ ratio encodes direction |
| Splicing rate $\beta$ | Degradation rate $\alpha$ |
| Transcription rate $\alpha$ | Synthesis rate $\beta$ |
| RNA velocity vector in PCA/UMAP | Protein velocity vector along pseudotime |

The key difference: the metabolic pulse provides an absolute timescale $T$, so **absolute** values of $\alpha$ and $\beta$ are recoverable, not just relative directions. RNA velocity gives only the direction of change; this method gives the rate.

---

## The Observable Signal: M/(L+M) Fraction

Define the **new-protein fraction**:

$$f_p(i) = \frac{M_p(i)}{L_p(i) + M_p(i)}$$

At steady state ($L_0 = P_{ss}$), this fraction is constant along pseudotime:

$$f_p^{ss} = 1 - e^{-\alpha T}$$

Deviations from $f_p^{ss}$ along pseudotime reveal dynamics:

- $f_p > f_p^{ss}$: protein is **increasing** (more new protein than expected at steady state)
- $f_p < f_p^{ss}$: protein is **decreasing** (more old protein than expected)

This fraction is directly observable from mass spec data without any model fitting — a model-free diagnostic of protein acceleration.

---

## The $\delta\tau$ Calibration Problem

The pseudotime interval $\delta\tau$ that corresponds to the pulse duration $T$ is **unknown** and must be calibrated. Options:

1. **Steady-state proteins as internal calibration:** Proteins known to be at steady state across the trajectory should give $\hat{\alpha} = \alpha$ regardless of $\delta\tau$. The $\delta\tau$ that minimises variance in $\hat{\alpha}$ across steady-state proteins is the optimal calibration.

2. **Prior knowledge of $\alpha$:** If $\alpha$ is known for a subset of proteins (from independent turnover experiments), $\delta\tau$ can be inferred by inverting $\hat{\alpha} = -\ln(L/\hat{L}_0)/T$ with the known $\alpha$.

3. **Two pulses on a subset:** Running the full two-pulse experiment on a small subset of cells provides a direct $\alpha$ estimate, which calibrates $\delta\tau$ for the single-pulse population.

---

## Assumptions and Limitations

| Assumption | Implication if violated |
|---|---|
| $\alpha$, $\beta$ constant within pulse window | Edge effects at trajectory boundaries; $\hat{\alpha}$ biased if rates change within $T$ hours |
| Trajectory is one-dimensional | Multi-branching trajectories require branch-aware pseudotime; $\hat{L}_0$ will be averaged across branches |
| Cells evolve deterministically along trajectory | Stochastic cell fate (e.g., Raj lab EGFR fluctuations) violates the past-state assumption for rare cell states |
| $\delta\tau$ is known | Miscalibrated $\delta\tau$ biases $\hat{\alpha}$; see calibration options above |
| Enough cells at every pseudotime point | Sparse sampling at trajectory endpoints degrades $\hat{L}_0$ estimation |

---

## Why This Works Better Than RNA Velocity for This Question

1. **Metabolic labeling provides ground truth timescale** — $T$ is precisely known; RNA splicing rates vary across transcripts and cell types.

2. **Protein stability is the biological quantity of interest** — for drug resistance, differentiation, and memory formation, protein abundance (not mRNA) is the functional readout.

3. **Post-transcriptional dynamics are captured** — proteins whose regulation is primarily translational or post-translational (e.g., EGFR stabilisation by HSP90, ubiquitin-mediated degradation) are invisible to RNA velocity but directly quantified here.

4. **Absolute rate estimates** — $\hat{\alpha}$ and $\hat{\beta}$ in h$^{-1}$ and AU h$^{-1}$ respectively, not dimensionless ratios.
