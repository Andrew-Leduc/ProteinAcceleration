#!/usr/bin/env python3
"""
single_pulse_pseudotime_simulation.py

Single metabolic pulse + pseudotime recovery of protein dynamics.

Problem
-------
A single pulse (Light + Medium only) gives 2 measurements per protein per cell
but 3 unknowns (α, β, L₀) — under-determined at the single-cell level.

Solution
--------
Use a population of cells distributed along a biological trajectory.
Cells at earlier pseudotime represent the past state of later cells —
providing the missing L₀ estimate and making the system exactly determined.
This is the proteomics analogue of RNA velocity.

Algorithm
---------
1.  Measure (L_p, M_p) for N proteins across N_cells cells.
2.  Compute UMAP on total protein (L+M); assign diffusion pseudotime (DPT).
3.  For each cell i at pseudotime τ_i:
        L₀_p(i) ≈ weighted mean of Total_p(j) for cells j at τ_j ≈ τ_i − δτ
4.  Recover per-cell parameters:
        α_p   = −ln(L_p / L₀_p) / T
        β_p   = α_p · M_p / (1 − e^{−α_p T})
        P_ss  = β_p / α_p

Key assumption: δτ (pseudotime interval) corresponds to the pulse duration T.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
import scanpy as sc
import anndata as ad

np.random.seed(42)
rng = np.random.default_rng(42)
sc.settings.verbosity = 0

# ── Parameters ────────────────────────────────────────────────────────────────
N_CELLS   = 600
T_PULSE   = 8.0    # pulse duration (h)
ALPHA     = 0.1    # true degradation rate (h⁻¹), uniform across proteins
DELTA_TAU = 0.15   # pseudotime interval ≡ T_PULSE h  (calibration constant)
NOISE_CV  = 0.15   # lognormal measurement noise (15% CV)

N_STATIC  = 20     # proteins at steady state along trajectory
N_UP      = 5      # proteins increasing along trajectory
N_DOWN    = 5      # proteins decreasing along trajectory
N_PROTS   = N_STATIC + N_UP + N_DOWN

EPS = 1e-9

# ── Protein trajectory model ──────────────────────────────────────────────────
def sigmoid(x, k=10, x0=0.5):
    return 1.0 / (1.0 + np.exp(-k * (x - x0)))

def protein_ss(tau, kind, lo, hi):
    """Steady-state protein abundance as a function of pseudotime τ."""
    if kind == 'static': return np.full_like(tau, (lo + hi) / 2.0)
    if kind == 'up':     return lo + (hi - lo) * sigmoid(tau)
    if kind == 'down':   return hi - (hi - lo) * sigmoid(tau)

# Protein properties
kinds = ['static'] * N_STATIC + ['up'] * N_UP + ['down'] * N_DOWN
p_lo  = rng.uniform(50,  150, N_PROTS)
p_hi  = rng.uniform(400, 600, N_PROTS)

# True pseudotime (random draw across [0,1])
tau_true = rng.uniform(0, 1, N_CELLS)
tau_past = np.clip(tau_true - DELTA_TAU, 0, 1)

# Ground-truth protein levels: current and past steady state
P_cur  = np.column_stack([protein_ss(tau_true, kinds[p], p_lo[p], p_hi[p])
                           for p in range(N_PROTS)])
P_prev = np.column_stack([protein_ss(tau_past, kinds[p], p_lo[p], p_hi[p])
                           for p in range(N_PROTS)])

# ── Simulate L / M measurements ───────────────────────────────────────────────
L_true = P_prev * np.exp(-ALPHA * T_PULSE)
M_true = P_cur  * (1.0 - np.exp(-ALPHA * T_PULSE))

def add_noise(X):
    sigma = np.sqrt(np.log(1.0 + NOISE_CV ** 2))
    return X * np.exp(rng.normal(0, sigma, X.shape))

L     = add_noise(L_true)
M     = add_noise(M_true)
Total = L + M          # observable total protein per cell per protein

# ── UMAP + diffusion pseudotime ───────────────────────────────────────────────
adata = ad.AnnData(X=np.log1p(Total))
sc.pp.pca(adata, n_comps=15)
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca')
sc.tl.umap(adata)

# Root at cell with minimum total abundance across dynamic proteins
dyn_cols = list(range(N_STATIC, N_PROTS))
adata.uns['iroot'] = int(Total[:, dyn_cols].sum(axis=1).argmin())
sc.tl.dpt(adata)

tau_est = adata.obs['dpt_pseudotime'].values.copy()
tau_est = (tau_est - tau_est.min()) / (tau_est.max() - tau_est.min() + EPS)

# ── L₀ estimation via pseudotime neighbourhood ───────────────────────────────
def estimate_L0(tau, Total_obs, delta_tau, bw_frac=0.35):
    """
    Gaussian-weighted mean of Total_obs at pseudotime τ − δτ.
    Each row i: weights cells j by exp(−0.5·((τ_j − (τ_i−δτ))/bw)²).
    """
    target = tau[:, None] - delta_tau          # (N,1) target pseudotime
    dist   = tau[None, :] - target             # (N,N) signed distances
    bw     = delta_tau * bw_frac
    W      = np.exp(-0.5 * (dist / bw) ** 2)  # (N,N) Gaussian weights
    W     /= W.sum(axis=1, keepdims=True)
    return W @ Total_obs                       # (N,P) estimated L₀

L0_hat = estimate_L0(tau_est, Total, DELTA_TAU)

# ── Parameter recovery ────────────────────────────────────────────────────────
ratio      = np.clip(L / np.clip(L0_hat, EPS, None), EPS, 1.0)
alpha_hat  = -np.log(ratio) / T_PULSE
denom      = np.clip(1.0 - np.exp(-alpha_hat * T_PULSE), EPS, None)
beta_hat   = alpha_hat * M / denom
Pss_hat    = np.clip(beta_hat / np.clip(alpha_hat, EPS, None), 0, None)

# ── Smoothing helper ──────────────────────────────────────────────────────────
def smooth(y, win=40):
    """Simple boxcar smooth."""
    kernel = np.ones(win) / win
    return np.convolve(y, kernel, mode='valid')

def smooth_along(tau_sorted, values, win=40):
    ys = smooth(values, win)
    xs = tau_sorted[win // 2: win // 2 + len(ys)]
    return xs, ys

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 11))
gs  = GridSpec(2, 3, fig, hspace=0.44, wspace=0.38)

ax_u1   = fig.add_subplot(gs[0, 0])   # UMAP coloured by pseudotime
ax_u2   = fig.add_subplot(gs[0, 1])   # UMAP coloured by M-fraction
ax_mfrc = fig.add_subplot(gs[0, 2])   # M-fraction along pseudotime
ax_traj = fig.add_subplot(gs[1, 0])   # True vs recovered trajectory
ax_alph = fig.add_subplot(gs[1, 1])   # Recovered α along pseudotime
ax_corr = fig.add_subplot(gs[1, 2])   # True vs recovered Pss scatter

xy = adata.obsm['X_umap']
sort_est  = np.argsort(tau_est)
sort_true = np.argsort(tau_true)

# (A) UMAP — pseudotime -------------------------------------------------------
sc0 = ax_u1.scatter(*xy.T, c=tau_est, s=6, cmap='plasma', rasterized=True)
plt.colorbar(sc0, ax=ax_u1, label='DPT pseudotime', shrink=0.8)
ax_u1.set_title('(A)  UMAP — diffusion pseudotime', fontsize=10, fontweight='bold')
ax_u1.set_xlabel('UMAP 1'); ax_u1.set_ylabel('UMAP 2')

# (B) UMAP — M/(L+M) for an increasing protein --------------------------------
ex_up  = N_STATIC                   # first increasing protein
m_frac = M[:, ex_up] / np.clip(Total[:, ex_up], EPS, None)
sc1 = ax_u2.scatter(*xy.T, c=m_frac, s=6, cmap='RdYlGn', rasterized=True,
                     vmin=0, vmax=1)
plt.colorbar(sc1, ax=ax_u2, label='M / (L+M)', shrink=0.8)
ax_u2.set_title('(B)  UMAP — new-protein fraction\n(increasing protein)',
                fontsize=10, fontweight='bold')
ax_u2.set_xlabel('UMAP 1'); ax_u2.set_ylabel('UMAP 2')

# (C) M/(L+M) fraction along estimated pseudotime ----------------------------
steady_expected = 1.0 - np.exp(-ALPHA * T_PULSE)
palette = {'up': '#2ca02c', 'down': '#d62728', 'static': '#1f77b4'}
labels  = {'up': 'Increasing protein', 'down': 'Decreasing protein',
           'static': 'Steady-state protein'}

tau_s = tau_est[sort_est]
for kind, col_idx in [('up', N_STATIC), ('down', N_STATIC + N_UP), ('static', 0)]:
    frac = (M[:, col_idx] / np.clip(Total[:, col_idx], EPS, None))[sort_est]
    xs, ys = smooth_along(tau_s, frac)
    ax_mfrc.plot(xs, ys, color=palette[kind], lw=2.2, label=labels[kind])

ax_mfrc.axhline(steady_expected, color='k', ls='--', lw=1.5,
                label=f'Steady-state ≈ {steady_expected:.2f}\n(= 1−e^{{−αT}})')
ax_mfrc.set_xlabel('Pseudotime', fontsize=9)
ax_mfrc.set_ylabel('M / (L+M)', fontsize=9)
ax_mfrc.set_title('(C)  New-protein fraction along pseudotime', fontsize=10,
                  fontweight='bold')
ax_mfrc.legend(fontsize=8, framealpha=0.9)
ax_mfrc.set_xlim(0, 1); ax_mfrc.set_ylim(0, 1)

# (D) True vs recovered protein trajectory ------------------------------------
tau_ts = tau_true[sort_true]
for kind, col_idx, color in [('up',   N_STATIC,       '#2ca02c'),
                               ('down', N_STATIC+N_UP,  '#d62728')]:
    # True P_ss
    ax_traj.plot(tau_ts, P_cur[sort_true, col_idx],
                 color=color, lw=2.5, label=f'{labels[kind]} (true)')
    # Recovered P_ss (smoothed along estimated pseudotime)
    phat = Pss_hat[sort_est, col_idx]
    xs, ys = smooth_along(tau_s, phat)
    ax_traj.plot(xs, ys, color=color, lw=2, ls='--',
                 label=f'{labels[kind]} (recovered)')

ax_traj.set_xlabel('Pseudotime', fontsize=9)
ax_traj.set_ylabel('Protein abundance (AU)', fontsize=9)
ax_traj.set_title('(D)  True vs recovered protein trajectory', fontsize=10,
                  fontweight='bold')
ax_traj.legend(fontsize=8, framealpha=0.9)

# (E) Recovered α along estimated pseudotime ----------------------------------
for kind, col_idx, color in [('up',     N_STATIC,       '#2ca02c'),
                               ('down',   N_STATIC+N_UP,  '#d62728'),
                               ('static', 0,              '#1f77b4')]:
    ahat = alpha_hat[sort_est, col_idx]
    xs, ys = smooth_along(tau_s, np.clip(ahat, 0, 0.5))
    ax_alph.plot(xs, ys, color=color, lw=2, label=labels[kind])

ax_alph.axhline(ALPHA, color='k', ls='--', lw=1.5,
                label=f'True α = {ALPHA} h⁻¹')
ax_alph.set_xlabel('Pseudotime', fontsize=9)
ax_alph.set_ylabel('Recovered α (h⁻¹)', fontsize=9)
ax_alph.set_title('(E)  Recovered degradation rate α', fontsize=10,
                  fontweight='bold')
ax_alph.legend(fontsize=8, framealpha=0.9)

# (F) True vs recovered P_ss per protein (mid-trajectory cells) ---------------
mid = (tau_est > 0.4) & (tau_est < 0.6)
true_Pss  = P_cur[mid, :].mean(axis=0)
recov_Pss = Pss_hat[mid, :].mean(axis=0)

c_scatter = (['#1f77b4'] * N_STATIC + ['#2ca02c'] * N_UP + ['#d62728'] * N_DOWN)
ax_corr.scatter(true_Pss, recov_Pss, c=c_scatter, s=45, alpha=0.85, zorder=3)

lo = min(true_Pss.min(), recov_Pss.min()) * 0.85
hi = max(true_Pss.max(), recov_Pss.max()) * 1.10
ax_corr.plot([lo, hi], [lo, hi], 'k--', lw=1.5)
ax_corr.set_xlabel('True P_ss (AU)', fontsize=9)
ax_corr.set_ylabel('Recovered P_ss (AU)', fontsize=9)
ax_corr.set_title('(F)  Recovery accuracy\n(cells at τ = 0.4–0.6)',
                  fontsize=10, fontweight='bold')

r = np.corrcoef(true_Pss, recov_Pss)[0, 1]
ax_corr.text(0.05, 0.95, f'r = {r:.3f}', transform=ax_corr.transAxes,
             fontsize=11, va='top', fontweight='bold')
ax_corr.legend(handles=[
    Patch(fc='#1f77b4', label='Steady-state'),
    Patch(fc='#2ca02c', label='Increasing'),
    Patch(fc='#d62728', label='Decreasing'),
    plt.Line2D([0],[0], c='k', ls='--', label='Identity'),
], fontsize=8, framealpha=0.9)

fig.suptitle(
    'Single Metabolic Pulse + Pseudotime  |  Protein Dynamics Recovery\n'
    fr'α = {ALPHA} h⁻¹,  T_pulse = {T_PULSE} h,  δτ = {DELTA_TAU},  '
    fr'N_cells = {N_CELLS},  N_proteins = {N_PROTS}',
    fontsize=11, fontweight='bold', y=1.01
)

plt.savefig('single_pulse_pseudotime.png', dpi=150, bbox_inches='tight')
plt.savefig('single_pulse_pseudotime.pdf', bbox_inches='tight')
print('Saved single_pulse_pseudotime.png / .pdf')
plt.show()
