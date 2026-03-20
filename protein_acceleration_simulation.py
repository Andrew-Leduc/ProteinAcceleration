#!/usr/bin/env python3
"""
protein_acceleration_simulation.py

Dual metabolic pulse labelling simulation for single-cell proteomics.

Experimental design
-------------------
  t = 0–8 h   : Pulse 1 — medium amino acids added (new protein → M pool)
  t = 8–16 h  : Pulse 2 — heavy  amino acids added (new protein → H pool)
  Pre-existing protein at t = 0 → L pool (decays, never re-synthesised)

ODE system
----------
Phase 1  (0 ≤ t < 8 h):
    dL/dt = −α · L(t)
    dM/dt =  β − α · M(t)      ← synthesis + degradation
    dH/dt =  0

Phase 2  (8 ≤ t ≤ 16 h):
    dL/dt = −α · L(t)
    dM/dt = −α · M(t)          ← medium no longer synthesised
    dH/dt =  β − α · H(t)      ← synthesis + degradation

Analytical solutions
--------------------
    L(t)    = L₀ · e^{−αt}
    M(16)   = (β/α)(1 − e^{−8α}) · e^{−8α}      ← independent of L₀!
    H(16)   = (β/α)(1 − e^{−8α})                 ← independent of L₀!
    Total   = β/α + (L₀ − β/α) · e^{−αt}

Key insight: H and M at endpoint are identical across all three scenarios
(same α, β). Only L(16) reflects L₀, revealing whether the cell's protein
was increasing (L₀ < P_ss), at steady state (L₀ = P_ss), or
decreasing (L₀ > P_ss).
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

# ── Parameters ────────────────────────────────────────────────────────────────
ALPHA  = 0.1          # degradation rate constant (h⁻¹); half-life ≈ 6.9 h
P_SS   = 500          # steady-state protein abundance = β/α  (arbitrary units)
BETA   = ALPHA * P_SS # synthesis rate = 50 AU h⁻¹

T1, T2 = 8.0, 16.0   # end of pulse 1 / end of experiment (h)
DT     = 0.05         # Euler integration step (h)
t      = np.arange(0, T2 + DT, DT)

# Colour palette
CL = '#636363'   # Light  — grey
CM = '#E07B39'   # Medium — orange
CH = '#5178A8'   # Heavy  — blue
CLINE = 'black'  # Total protein line

SCENARIOS = [
    dict(name='Increasing',   L0=100,  col='#2ca02c'),
    dict(name='Steady State', L0=500,  col='#1f77b4'),
    dict(name='Decreasing',   L0=1000, col='#d62728'),
]


# ── ODE Integration ───────────────────────────────────────────────────────────
def simulate(L0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Forward Euler integration of the dual-pulse system."""
    L, M, H = np.zeros(len(t)), np.zeros(len(t)), np.zeros(len(t))
    L[0] = L0
    for i in range(len(t) - 1):
        if t[i] < T1:                             # Pulse 1: medium
            dL = -ALPHA * L[i]
            dM =  BETA  - ALPHA * M[i]
            dH =  0.0
        else:                                     # Pulse 2: heavy
            dL = -ALPHA * L[i]
            dM = -ALPHA * M[i]
            dH =  BETA  - ALPHA * H[i]
        L[i+1] = L[i] + dL * DT
        M[i+1] = M[i] + dM * DT
        H[i+1] = H[i] + dH * DT
    return L, M, H


# ── Run simulations ───────────────────────────────────────────────────────────
for s in SCENARIOS:
    L, M, H = simulate(s['L0'])
    s.update(L=L, M=M, H=H, total=L + M + H)


# ── Figure layout ─────────────────────────────────────────────────────────────
#
#   Row 0  │ [Increasing]  │ [Steady State] │ [Decreasing]  │   ← L/M/H traces
#   Row 1  │        [Total protein comparison — full width]  │
#   Row 2  │        [Endpoint bar chart     — full width]    │
#
fig = plt.figure(figsize=(15, 12))
gs  = GridSpec(3, 3, figure=fig,
               hspace=0.52, wspace=0.34,
               height_ratios=[1.5, 1.1, 0.85])

axes_top = [fig.add_subplot(gs[0, j]) for j in range(3)]
ax_total = fig.add_subplot(gs[1, :])
ax_bar   = fig.add_subplot(gs[2, :2])   # absolute amounts (left 2/3)
ax_frac  = fig.add_subplot(gs[2, 2])    # relative fractions (right 1/3)


# ── Helper: shade pulse windows ───────────────────────────────────────────────
def shade_pulses(ax, fontsize=8):
    ax.axvspan(0,  T1, alpha=0.07, color=CM, zorder=0)
    ax.axvspan(T1, T2, alpha=0.07, color=CH, zorder=0)
    ax.axvline(T1, color='k', lw=0.9, ls='--', zorder=1)
    yl = ax.get_ylim()
    ypos = yl[0] + 0.97 * (yl[1] - yl[0])
    ax.text(T1 / 2,       ypos, 'Pulse 1\n(Medium)', ha='center', va='top',
            fontsize=fontsize, color='#7B3F00')
    ax.text((T1 + T2) / 2, ypos, 'Pulse 2\n(Heavy)',  ha='center', va='top',
            fontsize=fontsize, color='#1E3F6F')


# ── Row 0: L / M / H dynamics for each scenario ──────────────────────────────
for ax, s in zip(axes_top, SCENARIOS):
    ax.plot(t, s['L'],     color=CL, lw=2.2, label='Light (L)')
    ax.plot(t, s['M'],     color=CM, lw=2.2, label='Medium (M)')
    ax.plot(t, s['H'],     color=CH, lw=2.2, label='Heavy (H)')
    ax.plot(t, s['total'], color=s['col'], lw=2.0, ls=':', label='Total')
    ax.axhline(P_SS, color='k', lw=0.8, ls='--', alpha=0.35, label=f'P_ss={P_SS}')

    ax.set_xlim(0, T2)
    ax.set_ylim(bottom=0)
    ax.set_title(f"{s['name']}\n(L₀ = {s['L0']} AU)",
                 fontsize=10.5, fontweight='bold', color=s['col'])
    ax.set_xlabel('Time (h)', fontsize=9)
    ax.set_ylabel('Protein abundance (AU)', fontsize=9)
    shade_pulses(ax, fontsize=7)

    # Endpoint readout box
    H16, M16, L16 = s['H'][-1], s['M'][-1], s['L'][-1],
    tot16 = s['total'][-1]
    txt = (f"t = 16 h\n"
           f"  H = {H16:.1f}\n"
           f"  M = {M16:.1f}\n"
           f"  L = {L16:.1f}\n"
           f"  Total = {tot16:.1f}")
    ax.text(0.97, 0.97, txt, transform=ax.transAxes, fontsize=7.5,
            ha='right', va='top', family='monospace',
            bbox=dict(boxstyle='round', fc='white', alpha=0.80, ec='#bbbbbb'))

# Shared legend (on middle panel)
handles_lmh = [
    mpatches.Patch(fc=CL, label='Light — pre-existing (L)'),
    mpatches.Patch(fc=CM, label='Medium — pulse 1 synthesis (M)'),
    mpatches.Patch(fc=CH, label='Heavy — pulse 2 synthesis (H)'),
    plt.Line2D([0], [0], color='k', ls=':', lw=2, label='Total protein'),
    plt.Line2D([0], [0], color='k', ls='--', lw=1, alpha=0.5, label=f'P_ss = {P_SS} AU'),
]
axes_top[1].legend(handles=handles_lmh, fontsize=8, loc='lower right', framealpha=0.9)


# ── Row 1: Total protein comparison ──────────────────────────────────────────
for s in SCENARIOS:
    ax_total.plot(t, s['total'], color=s['col'], lw=2.5,
                  label=f"{s['name']}  (L₀ = {s['L0']})")

ax_total.axhline(P_SS, color='k', lw=1.2, ls='--', alpha=0.5,
                 label=f'P_ss = β/α = {P_SS} AU')
ax_total.set_xlim(0, T2)
ax_total.set_xlabel('Time (h)', fontsize=10)
ax_total.set_ylabel('Total protein (AU)', fontsize=10)
ax_total.set_title('Total protein dynamics across scenarios', fontsize=11, fontweight='bold')
ax_total.legend(fontsize=9, framealpha=0.88)
shade_pulses(ax_total, fontsize=8.5)


# ── Row 2 left: Absolute amounts at t = 16 h ─────────────────────────────────
bar_width = 0.22
x = np.arange(len(SCENARIOS))
labels = [s['name'] for s in SCENARIOS]

H_vals = [s['H'][-1] for s in SCENARIOS]
M_vals = [s['M'][-1] for s in SCENARIOS]
L_vals = [s['L'][-1] for s in SCENARIOS]

bars_H = ax_bar.bar(x - bar_width, H_vals, bar_width, color=CH, label='Heavy (H)',  zorder=3)
bars_M = ax_bar.bar(x,             M_vals, bar_width, color=CM, label='Medium (M)', zorder=3)
bars_L = ax_bar.bar(x + bar_width, L_vals, bar_width, color=CL, label='Light (L)',  zorder=3)

ax_bar.set_xticks(x)
ax_bar.set_xticklabels(labels, fontsize=10)
ax_bar.set_ylabel('Protein abundance at t = 16 h (AU)', fontsize=9)
ax_bar.set_title('Absolute amounts at sampling\n(H and M identical — only L differs)',
                 fontsize=10, fontweight='bold')
ax_bar.legend(fontsize=8.5, framealpha=0.88)
ax_bar.set_ylim(bottom=0)
ax_bar.grid(axis='y', alpha=0.3, zorder=0)

for bars in (bars_H, bars_M, bars_L):
    for bar in bars:
        h = bar.get_height()
        ax_bar.text(bar.get_x() + bar.get_width() / 2, h + 2, f'{h:.0f}',
                    ha='center', va='bottom', fontsize=7)


# ── Row 2 right: Relative fractions at t = 16 h (stacked bars) ───────────────
totals = [H + M + L for H, M, L in zip(H_vals, M_vals, L_vals)]
H_frac = [H / tot * 100 for H, tot in zip(H_vals, totals)]
M_frac = [M / tot * 100 for M, tot in zip(M_vals, totals)]
L_frac = [L / tot * 100 for L, tot in zip(L_vals, totals)]

x2 = np.arange(len(SCENARIOS))
ax_frac.bar(x2, H_frac, color=CH, label='Heavy (H)',  zorder=3)
ax_frac.bar(x2, M_frac, color=CM, label='Medium (M)', bottom=H_frac, zorder=3)
ax_frac.bar(x2, L_frac, color=CL, label='Light (L)',
            bottom=[h + m for h, m in zip(H_frac, M_frac)], zorder=3)

# Fraction labels inside stacked bars
for i, (hf, mf, lf) in enumerate(zip(H_frac, M_frac, L_frac)):
    ax_frac.text(i, hf / 2,              f'{hf:.1f}%', ha='center', va='center',
                 fontsize=8, color='white', fontweight='bold')
    ax_frac.text(i, hf + mf / 2,         f'{mf:.1f}%', ha='center', va='center',
                 fontsize=8, color='white', fontweight='bold')
    ax_frac.text(i, hf + mf + lf / 2,    f'{lf:.1f}%', ha='center', va='center',
                 fontsize=8, color='white', fontweight='bold')

ax_frac.set_xticks(x2)
ax_frac.set_xticklabels(labels, fontsize=9)
ax_frac.set_ylabel('Relative fraction at t = 16 h (%)', fontsize=9)
ax_frac.set_ylim(0, 100)
ax_frac.set_title('Isotope composition\nat sampling', fontsize=10, fontweight='bold')
ax_frac.grid(axis='y', alpha=0.3, zorder=0)
ax_frac.legend(fontsize=8, framealpha=0.88, loc='upper right')

# ── Supertitle ───────────────────────────────────────────────────────────────
fig.suptitle(
    'Protein Acceleration — Dual Metabolic Pulse Simulation\n'
    r'$\alpha = 0.1\ h^{-1}$ (half-life $\approx 6.9\ h$),  '
    r'$\beta = 50\ AU\ h^{-1}$,  $P_{ss} = \beta/\alpha = 500\ AU$',
    fontsize=12, fontweight='bold', y=1.01
)

# ── Save ─────────────────────────────────────────────────────────────────────
fig.savefig('protein_acceleration_simulation.png', dpi=150, bbox_inches='tight')
fig.savefig('protein_acceleration_simulation.pdf', bbox_inches='tight')
print("Saved protein_acceleration_simulation.png / .pdf")
plt.show()
