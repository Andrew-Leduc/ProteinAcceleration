# Dual Metabolic Pulse — Steady State Derivation

## Setup

Two sequential 8-hour pulses, measurement at $t = 16$ h.

| Symbol | Meaning |
|--------|---------|
| $\alpha$ | Degradation rate constant (h$^{-1}$) |
| $\beta$ | Synthesis rate (AU h$^{-1}$) |
| $P_{ss} = \beta/\alpha$ | Steady-state protein abundance |
| $L_0$ | Initial protein abundance (all Light at $t=0$) |
| $T_1 = 8$ h | End of pulse 1 (Medium) |
| $T_2 = 16$ h | End of pulse 2 (Heavy), sampling time |

**Steady state condition:** $L_0 = P_{ss} = \beta / \alpha$

---

## Phase 1: Medium pulse ($0 \leq t \leq T_1$)

$$\frac{dL}{dt} = -\alpha L, \qquad L(0) = L_0$$

$$\frac{dM}{dt} = \beta - \alpha M, \qquad M(0) = 0$$

$$\frac{dH}{dt} = 0, \qquad H(0) = 0$$

**Analytical solutions:**

$$L(t) = L_0\, e^{-\alpha t}$$

$$M(t) = \frac{\beta}{\alpha}\!\left(1 - e^{-\alpha t}\right) = P_{ss}\!\left(1 - e^{-\alpha t}\right)$$

**Total protein at end of pulse 1** (substituting $L_0 = P_{ss}$):

$$\text{Total}(T_1) = L_0\,e^{-\alpha T_1} + P_{ss}\!\left(1 - e^{-\alpha T_1}\right) = P_{ss}\,e^{-\alpha T_1} + P_{ss} - P_{ss}\,e^{-\alpha T_1} = P_{ss} \checkmark$$

The total is already conserved at $T_1$ — a direct consequence of $L_0 = P_{ss}$.

---

## Phase 2: Heavy pulse ($T_1 \leq t \leq T_2$, let $\tau = t - T_1$)

$$\frac{dL}{dt} = -\alpha L, \qquad L(T_1) = P_{ss}\,e^{-\alpha T_1}$$

$$\frac{dM}{dt} = -\alpha M, \qquad M(T_1) = P_{ss}\!\left(1 - e^{-\alpha T_1}\right)$$

$$\frac{dH}{dt} = \beta - \alpha H, \qquad H(T_1) = 0$$

**Analytical solutions (in $\tau$):**

$$L(T_1 + \tau) = P_{ss}\,e^{-\alpha(T_1 + \tau)}$$

$$M(T_1 + \tau) = P_{ss}\!\left(1 - e^{-\alpha T_1}\right)e^{-\alpha\tau}$$

$$H(T_1 + \tau) = P_{ss}\!\left(1 - e^{-\alpha\tau}\right)$$

---

## Endpoint values at $t = T_2$ (i.e. $\tau = T_1 = 8$ h)

$$\boxed{L(T_2) = P_{ss}\,e^{-\alpha T_2}}$$

$$\boxed{M(T_2) = P_{ss}\!\left(1 - e^{-\alpha T_1}\right)e^{-\alpha T_1}}$$

$$\boxed{H(T_2) = P_{ss}\!\left(1 - e^{-\alpha T_1}\right)}$$

**Verify total is conserved:**

$$\text{Total}(T_2) = P_{ss}\!\left[\,e^{-\alpha T_2} + \left(1 - e^{-\alpha T_1}\right)e^{-\alpha T_1} + \left(1 - e^{-\alpha T_1}\right)\right]$$

$$= P_{ss}\!\left[\,e^{-2\alpha T_1} + e^{-\alpha T_1} - e^{-2\alpha T_1} + 1 - e^{-\alpha T_1}\right] = P_{ss} \checkmark$$

---

## What we can infer from the measurement

### 1. Estimating $\alpha$ from the H/M ratio

$$\frac{H(T_2)}{M(T_2)} = \frac{P_{ss}(1 - e^{-\alpha T_1})}{P_{ss}(1 - e^{-\alpha T_1})\,e^{-\alpha T_1}} = e^{\alpha T_1}$$

$$\boxed{\alpha = \frac{\ln\!\left(H / M\right)}{T_1}}$$

This is independent of $\beta$, $L_0$, and whether the protein is at steady state. The H/M ratio is a **pure readout of the degradation rate**.

### 2. Estimating $\beta$ from H

$$H(T_2) = \frac{\beta}{\alpha}\!\left(1 - e^{-\alpha T_1}\right) \implies \boxed{\beta = \frac{\alpha\, H(T_2)}{1 - e^{-\alpha T_1}}}$$

### 3. Estimating $L_0$ from L

$$L(T_2) = L_0\, e^{-\alpha T_2} \implies \boxed{L_0 = L(T_2)\,e^{\alpha T_2}}$$

### 4. Steady state test

Compare the recovered $L_0$ to $P_{ss} = \beta/\alpha$:

$$\Delta = L_0 - P_{ss} = L(T_2)\,e^{\alpha T_2} - \frac{H(T_2)}{1 - e^{-\alpha T_1}}$$

| $\Delta$ | Interpretation |
|----------|----------------|
| $\Delta = 0$ | Steady state — protein abundance unchanged |
| $\Delta < 0$ | Protein **increasing** — $L_0 < P_{ss}$, synthesis outpaces starting level |
| $\Delta > 0$ | Protein **decreasing** — $L_0 > P_{ss}$, starting level exceeds synthesis capacity |

---

## Isotope fractions at steady state

At steady state ($L_0 = P_{ss}$), each fraction depends only on $\alpha$ and $T_1$:

$$f_L = \frac{L(T_2)}{P_{ss}} = e^{-\alpha T_2}$$

$$f_M = \frac{M(T_2)}{P_{ss}} = \left(1 - e^{-\alpha T_1}\right)e^{-\alpha T_1}$$

$$f_H = \frac{H(T_2)}{P_{ss}} = 1 - e^{-\alpha T_1}$$

**Verification:** $f_L + f_M + f_H = e^{-\alpha T_2} + e^{-\alpha T_1} - e^{-\alpha T_2} + 1 - e^{-\alpha T_1} = 1\ \checkmark$

**Numerical example** ($\alpha = 0.1\ \text{h}^{-1}$, $T_1 = 8\ \text{h}$, $T_2 = 16\ \text{h}$):

$$f_H = 1 - e^{-0.8} \approx 55.1\%$$

$$f_M = (1 - e^{-0.8})\,e^{-0.8} \approx 24.8\%$$

$$f_L = e^{-1.6} \approx 20.2\%$$

$$\frac{H}{M} = e^{0.8} \approx 2.23 \implies \alpha = \frac{\ln 2.23}{8} = 0.1\ \text{h}^{-1}\ \checkmark$$
