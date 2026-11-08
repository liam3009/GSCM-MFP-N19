import warnings; warnings.filterwarnings("ignore")
import os, numpy as np, pandas as pd
import matplotlib; import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import font_manager
display = Image = lambda *a, **k: None

_fams = {f.name for f in font_manager.fontManager.ttflist}
FONT = next((f for f in ["Liberation Sans", "Arial", "DejaVu Sans"] if f in _fams), "DejaVu Sans")
plt.rcParams.update({"font.family": FONT, "font.size": 10, "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#444444",
                     "axes.linewidth": 0.8, "xtick.color": "#333333", "ytick.color": "#333333"})
INK, MUTED, EDGE, GRID = "#1a1a1a", "#555555", "#333333", "#e3e3e3"
C_MP, C_FP = "#1f6fb2", "#c8553d"
O = "output/"; F = "figures/"; os.makedirs(F, exist_ok=True)
N_FULL = int(pd.read_csv(O + "T1_sample_flow.csv").set_index("Step").loc["Consented (full sample)", "n"])
N_CLEAN = int(pd.read_csv(O + "T1_sample_flow.csv").set_index("Step").loc["Clean sample", "n"])

def save(fig, name):
    for ext in ("png", "pdf", "tiff"):
        kw = {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tiff" else {"dpi": 300}
        fig.savefig(F + name + "." + ext, bbox_inches="tight", facecolor="white", **kw)
    plt.close(fig); display(Image(F + name + ".png", width=760))

def fmt(x):
    s = f"{x:.2f}"; s = s.replace("0.", ".", 1) if abs(x) < 1 else s; return s.replace("-", "−")

lat = pd.read_csv(O + "X4_latent_paths.csv")
def b(s, l, r):
    row = lat[(lat.Sample == s) & (lat.lhs == l) & (lat.rhs == r)].iloc[0]; return row["est.std"], row["pvalue"]


def box(ax, x, y, w, h, t, sub=None, bold=True, ls="-", fs=10.5):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0,rounding_size=0.08", fc="white", ec=EDGE, lw=1.2, ls=ls))
    if sub:
        ax.text(x, y + .14, t, ha="center", va="center", fontsize=fs, fontweight="bold" if bold else "normal", color=INK)
        ax.text(x, y - .18, sub, ha="center", va="center", fontsize=8.5, color=MUTED)
    else: ax.text(x, y, t, ha="center", va="center", fontsize=fs, fontweight="bold" if bold else "normal", color=INK)
def arr(ax, p, q, lw, color=EDGE, ls="-", rad=0.0):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=12, lw=lw, color=color, ls=ls,
                                 connectionstyle=f"arc3,rad={rad}", shrinkA=0, shrinkB=0))

# ---------------- Figure 1: conceptual model ----------------
fig, ax = plt.subplots(figsize=(7.4, 4.5)); ax.set_xlim(0.1, 10.3); ax.set_ylim(-0.05, 6.0); ax.axis("off")
gx, ox = 2.35, 8.65
box(ax, gx, 5.0, 4.3, 0.95, "Green supply chain management", "green design · purchasing · internal mgmt", fs=10)
box(ax, gx, 1.2, 4.3, 0.95, "Lean management", "lean quality · just-in-time · TPM", fs=10)
box(ax, gx, 3.1, 3.1, 0.85, "GSCM × Lean", "H3a synergy / H3b trade-off", bold=False, ls=(0, (4, 3)), fs=10)
box(ax, ox, 5.0, 3.0, 0.8, "Market performance")
box(ax, ox, 1.2, 3.0, 0.8, "Financial performance")
box(ax, ox, 3.1, 3.0, 0.95, "Controls", "CSR · size · industry · position", bold=False, ls=(0, (1, 2)), fs=9.5)
L, R = gx + 2.15, ox - 1.5
arr(ax, (L, 5.0), (R, 5.0), 2.6, C_MP)                       # GSCM -> MP (H1, larger side)
arr(ax, (L, 1.2), (R, 1.2), 2.6, C_FP)                       # LEAN -> FP (H2, larger side)
arr(ax, (L, 4.7), (R, 1.45), 1.6, "#5a5a5a")                 # GSCM -> FP (reference path in H1 and H2)
arr(ax, (L, 1.5), (R, 4.75), 0.9, "#9a9a9a")                 # LEAN -> MP (not hypothesised)
arr(ax, (gx + 1.55, 3.25), (R, 4.85), 1.0, EDGE, (0, (4, 3)), -0.12)
arr(ax, (gx + 1.55, 2.95), (R, 1.35), 1.0, EDGE, (0, (4, 3)), 0.12)
ax.text((L + R) / 2, 5.36, "H1", ha="center", va="bottom", fontsize=10, color=C_MP, fontweight="bold")
ax.text((L + R) / 2, 5.12, "> GSCM → financial", ha="center", va="bottom", fontsize=8, color=C_MP)
ax.text((L + R) / 2, 1.08, "> GSCM → financial", ha="center", va="top", fontsize=8, color=C_FP)
ax.text((L + R) / 2, 0.80, "H2", ha="center", va="top", fontsize=10, color=C_FP, fontweight="bold")
ax.text(5.0, 0.18, "Thick arrows: the larger path in each contrast. Dark grey: GSCM → financial, the reference path shared by H1 and H2.",
        ha="center", fontsize=7.8, color=MUTED)
save(fig, "Figure1_conceptual_model")

# ---------------- Figure 2: sample screening flow ----------------
fl = pd.read_csv(O + "T1_sample_flow.csv").set_index("Step")["n"]
n_raw, n_nc, n_sd0 = int(fl["Raw responses"]), int(fl["Did not consent"]), int(fl["Straightlining (SD = 0)"])
n_car = int(fl[[k for k in fl.index if k.startswith("Longstring")][0]]); n_dup = int(fl["Exact duplicate copies (all)"])
fig, ax = plt.subplots(figsize=(7.0, 4.4)); ax.set_xlim(0, 11); ax.set_ylim(0.7, 7.0); ax.axis("off")
mx = 2.6
box(ax, mx, 6.3, 4.4, 0.8, "Responses in public dataset", f"N = {n_raw}", fs=10)
box(ax, mx, 4.3, 4.4, 0.8, "Full sample (consented)", f"n = {N_FULL}", fs=10)
box(ax, mx, 1.4, 4.4, 0.8, "Clean sample", f"n = {N_CLEAN}", fs=10)
box(ax, 8.3, 5.3, 5.2, 0.62, f"Did not consent: n = {n_nc}", bold=False, fs=9.5)
box(ax, 8.3, 2.85, 5.2, 1.5, "", bold=False)
x0 = 5.95
ax.text(x0, 3.3, f"Careless responding: n = {n_car}", fontsize=9.5, fontweight="bold", color=INK, va="center")
ax.text(x0, 2.92, f"• same answer to all 34 items: {n_sd0}", fontsize=8.5, color=MUTED, va="center")
ax.text(x0, 2.6, f"• same answer to ≥ 20 consecutive items: {n_car - n_sd0}", fontsize=8.5, color=MUTED, va="center")
ax.text(x0, 2.28, f"  (all {n_dup} exact duplicates fall in these groups)", fontsize=8.5, color=MUTED, va="center")
arr(ax, (mx, 5.9), (mx, 4.7), 1.3); arr(ax, (mx, 3.9), (mx, 1.8), 1.3)
arr(ax, (mx, 5.3), (5.7, 5.3), 1.3); arr(ax, (mx, 2.85), (5.7, 2.85), 1.3)
save(fig, "Figure2_sample_flow")

# ---------------- Figure 3: structural model results ----------------
fig, ax = plt.subplots(figsize=(7.4, 4.0)); ax.set_xlim(0.1, 10.3); ax.set_ylim(0.3, 5.9); ax.axis("off")
gx, ox = 2.35, 8.65
r2f = lat[lat.Sample == "full"].iloc[0]; r2c = lat[lat.Sample == "clean"].iloc[0]
box(ax, gx, 4.7, 4.3, 0.9, "Green supply chain mgmt", "second-order: GD, GP, IEM")
box(ax, gx, 1.3, 4.3, 0.9, "Lean management", "second-order: LQ, JIT, TPM")
box(ax, ox, 4.7, 3.1, 0.9, "Market performance", f"R² = {fmt(r2f.R2_MP)} [{fmt(r2c.R2_MP)}]")
box(ax, ox, 1.3, 3.1, 0.9, "Financial performance", f"R² = {fmt(r2f.R2_FP)} [{fmt(r2c.R2_FP)}]")
L, R = gx + 2.15, ox - 1.55
paths = [("MP", "GSCM", (L, 4.7), (R, 4.7), 2.6, (6.05, 4.88), "bottom"),
         ("FP", "LEAN", (L, 1.3), (R, 1.3), 2.6, (6.05, 1.12), "top"),
         ("FP", "GSCM", (L, 4.42), (R, 1.55), 1.0, (L + 0.22 * (R - L), 4.42 - 0.22 * 2.87), "center"),
         ("MP", "LEAN", (L, 1.58), (R, 4.45), 1.0, (L + 0.22 * (R - L), 1.58 + 0.22 * 2.87), "center")]
for out, pred, p, q, lw, (tx, ty), va in paths:
    ef, pf = b("full", out, pred); ec, pc = b("clean", out, pred)
    arr(ax, p, q, lw, EDGE if lw > 2 else "#8a8a8a")
    ax.text(tx, ty, f"β = {fmt(ef)} [{fmt(ec)}]", ha="center", va=va, fontsize=9.5 if lw > 2 else 9,
            fontweight="bold" if lw > 2 else "normal", color=INK,
            bbox=None if lw > 2 else dict(fc="white", ec="none", pad=1.5))
ax.text(5.0, 5.55, f"Standardized estimates: full sample (n = {N_FULL}) [clean sample, n = {N_CLEAN}]", ha="center", fontsize=8.5, color=MUTED)
save(fig, "Figure3_structural_model")

# ---------------- Figure 4: robustness of the H1 and H2 contrasts ----------------
t6 = pd.read_csv(O + "T6_RQ1_latent_SEM.csv"); rob = pd.read_csv(O + "X5_robustness.csv")
rows = []
for s, lab in [("full", "Latent SEM, full sample (main)"), ("clean", "Latent SEM, clean sample")]:
    g = t6[(t6.Sample == s) & (t6.Parameter == "H1_GSCM_MPminusFP")].iloc[0]; l = t6[(t6.Sample == s) & (t6.Parameter == "H2_FP_LEANminusGSCM")].iloc[0]
    n = N_FULL if s == "full" else N_CLEAN
    rows.append((lab, n, g.Estimate, g.CI_low, g.CI_high, l.Estimate, l.CI_low, l.CI_high))
bl = pd.read_csv(O + "B_latent_checks.csv")
for chk, s, lab in [("RQ1 SEM with method factor (ULMC)", "full", "Latent SEM + method factor, full sample"),
                    ("Without JIT2", "full", "Latent SEM without JIT2, full sample"), ("Without JIT2", "clean", "Latent SEM without JIT2, clean sample"),
                    ("Without multivariate outliers", "full", "Latent SEM without 48 outliers")]:
    q = bl[(bl.Check == chk) & (bl.Sample == s)].set_index("Parameter")
    rows.append((lab, int(q.n.iloc[0]), q.loc["H1", "Estimate"], q.loc["H1", "CI_low"], q.loc["H1", "CI_high"], q.loc["H2", "Estimate"], q.loc["H2", "CI_low"], q.loc["H2", "CI_high"]))
NLAT = len(rows)
for _, r in rob.iterrows():
    sp = r.Specification.replace("Main: full sample, CSR + controls", "full sample, CSR + controls")
    sp = sp if sp.startswith("SMEs") else sp[0].lower() + sp[1:]
    lab = "Composite: " + sp.replace("financial firms", "Financial firms")
    rows.append((lab, r.n, r.H1_diff, r.H1_lo, r.H1_hi, r.H2_diff, r.H2_lo, r.H2_hi))
df = pd.DataFrame(rows, columns=["lab", "n", "g", "gl", "gh", "l", "ll", "lh"])
fig, axs = plt.subplots(1, 2, figsize=(7.6, 5.6), sharey=True, gridspec_kw={"wspace": 0.10})
y = np.arange(len(df))[::-1]
for ax, (c, lo, hi, title, col) in zip(axs, [("g", "gl", "gh", "A. H1 — GSCM: market − financial", C_MP), ("l", "ll", "lh", "B. H2 — Financial: lean − GSCM", C_FP)]):
    ax.axvline(0, color="#888888", lw=0.9)
    for yi, (_, r) in zip(y, df.iterrows()):
        main = yi >= len(df) - NLAT
        ax.plot([r[lo], r[hi]], [yi, yi], color=col, lw=2.2 if main else 1.5, solid_capstyle="round")
        ax.plot(r[c], yi, "o" if main else "s", ms=7 if main else 5.5, color=col, mec="white", mew=1.2, zorder=3)
    ax.set_title(title, fontsize=10, fontweight="bold", color=INK, loc="left")
    ax.grid(axis="x", color=GRID, lw=0.7); ax.set_axisbelow(True)
    ax.set_xlabel("Difference in standardized β (95% bootstrap CI)", fontsize=9, color=INK)
    ax.tick_params(axis="y", length=0)
    ax.axhline(len(df) - NLAT - 0.5, color="#bbbbbb", lw=0.7, ls=(0, (3, 3)))
axs[0].set_yticks(y); axs[0].set_yticklabels([f"{a} (n = {int(n)})" for a, n in zip(df.lab, df.n)], fontsize=8.5)
for ax, (lo, hi) in zip(axs, [("gl", "gh"), ("ll", "lh")]): ax.set_xlim(min(df[lo].min(), 0) - 0.05, df[hi].max() + 0.05)
save(fig, "Figure5_robustness")


# ---------------- Figure 5: RQ2 interaction (equivalence test and simple slopes) ----------------
from matplotlib.lines import Line2D
ci5 = pd.read_csv(O + "X6_rq2_ci.csv"); SESOI = float(pd.read_csv(O + "settings.csv").SESOI.iloc[0])
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.3), gridspec_kw={"width_ratios": [1, 1.25], "wspace": 0.35})
ax = axs[0]
ax.axvspan(-SESOI, SESOI, color="#efefef", zorder=0); ax.axvline(0, color="#888888", lw=0.9)
for yi, (_, r) in zip([1, 0], ci5.iterrows()):
    ax.plot([r.lo95, r.hi95], [yi, yi], color="#9a9a9a", lw=1.2)
    ax.plot([r.lo90, r.hi90], [yi, yi], color=INK, lw=3, solid_capstyle="round")
    ax.plot(r.b, yi, "o", ms=6.5, color=INK, mec="white", mew=1.2, zorder=3)
ax.set_yticks([1, 0]); ax.set_yticklabels([f"Full sample\n(n = {N_FULL})", f"Clean sample\n(n = {N_CLEAN})"], fontsize=9)
ax.set_ylim(-0.6, 1.6); ax.tick_params(axis="y", length=0)
ax.set_xlabel("GSCM × Lean interaction (standardized β)", fontsize=9)
for xv in (-0.05, 0.05): ax.axvline(xv, color="#9a9a9a", lw=0.8, ls=(0, (2, 2)))
ax.text(0, 1.45, f"shaded: ±{fmt(SESOI)}; dotted: ±.05", ha="center", fontsize=8, color=MUTED)
ax.set_title("A. Interaction estimate", fontsize=10, fontweight="bold", loc="left", color=INK)
ax.legend(handles=[Line2D([0], [0], color=INK, lw=3, label="90% CI (TOST)"), Line2D([0], [0], color="#9a9a9a", lw=1.2, label="95% CI")],
          fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(0, -0.2), ncol=2)
sl = pd.read_csv(O + "X9_slopes.csv").set_index("term"); k = ["zGSCM", "zLEAN", "zGSCM:zLEAN"]
bb, V = sl["b"], sl[k].loc[k].values
ax = axs[1]; xs = np.linspace(-2, 2, 81)
for lv, ls, col, lab in [(-1, (0, (5, 3)), "#8a8a8a", "Low lean (−1 SD)"), (0, "-", "#5a5a5a", "Mean lean"), (1, "-", INK, "High lean (+1 SD)")]:
    yhat = bb["zGSCM"] * xs + bb["zLEAN"] * lv + bb["zGSCM:zLEAN"] * xs * lv
    se = np.array([np.sqrt(np.array([x, lv, x * lv]) @ V @ np.array([x, lv, x * lv])) for x in xs])
    if lv != 0: ax.fill_between(xs, yhat - 1.96 * se, yhat + 1.96 * se, color=col, alpha=0.10, lw=0)
    ax.plot(xs, yhat, ls=ls, color=col, lw=1.8 if lv == 1 else 1.4, label=lab)
    ax.text(2.05, yhat[-1], f"slope {fmt(bb['zGSCM'] + bb['zGSCM:zLEAN'] * lv)}", fontsize=8, color=col, va="center")
ax.set_xlim(-2, 2.75); ax.set_xticks([-2, -1, 0, 1, 2])
ax.set_xlabel("GSCM (standardized)", fontsize=9); ax.set_ylabel("Predicted performance (SD)", fontsize=9)
ax.grid(color=GRID, lw=0.7); ax.set_axisbelow(True); ax.legend(fontsize=8, frameon=False, loc="upper left")
ax.set_title("B. Simple slopes, full sample", fontsize=10, fontweight="bold", loc="left", color=INK)
save(fig, "Figure4_RQ2_interaction")


# ---------------- Figure S1: HTMT heatmap ----------------
order = ["GD", "GP", "IEM", "LQ", "JIT", "TPM", "CSR", "MP", "FP"]
names = {"GD": "Green design", "GP": "Green purchasing", "IEM": "Internal env. mgmt", "LQ": "Lean quality", "JIT": "Just-in-time",
         "TPM": "TPM", "CSR": "CSR", "MP": "Market perf.", "FP": "Financial perf."}
fig, axs = plt.subplots(1, 2, figsize=(9.2, 4.4), gridspec_kw={"wspace": 0.05})
cmap = LinearSegmentedColormap.from_list("seq", ["#f4f7fb", "#9cc0e0", "#1f5f9e"])
for ax, s, ttl in zip(axs, ["full", "clean"], [f"A. Full sample (n = {N_FULL})", f"B. Clean sample (n = {N_CLEAN})"]):
    h = pd.read_csv(O + f"T5_htmt_dimensions_{s}.csv", index_col=0).reindex(index=order, columns=order)
    M = h.values.astype(float); M = np.where(np.isnan(M), M.T, M)
    mask = np.triu(np.ones_like(M, bool))
    show = np.where(mask, np.nan, M)[1:, :-1]
    ax.imshow(show, cmap=cmap, vmin=0.5, vmax=1.0)
    for i in range(1, 9):
        for j in range(i):
            v = M[i, j]; ax.text(j, i - 1, fmt(v), ha="center", va="center", fontsize=7.5,
                                 color="white" if v > 0.85 else INK, fontweight="bold" if round(v, 2) >= 0.90 else "normal")
    ax.set_xticks(range(8)); ax.set_xticklabels([names[k] for k in order[:-1]], rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(8)); ax.set_yticklabels([names[k] for k in order[1:]] if s == "full" else [], fontsize=8)
    ax.tick_params(length=0); [sp.set_visible(False) for sp in ax.spines.values()]
    ax.set_title(ttl, fontsize=10, fontweight="bold", loc="left", color=INK)
save(fig, "FigureS1_HTMT")
print("Figures 1–5 and S1 saved to figures/ (PNG 300 dpi, TIFF 600 dpi, PDF)")
