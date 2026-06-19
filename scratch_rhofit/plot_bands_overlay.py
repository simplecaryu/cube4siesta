"""Band-overlay plot: converged baseline vs restart methods, Fermi-aligned.
  argv: title  out.png  ref.bands  [label::file.bands ...]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_bands(path):
    toks = open(path).read().split()
    ef = float(toks[0]); i = 5
    nb, ns, nk = int(toks[i]), int(toks[i+1]), int(toks[i+2]); i += 3
    per = 1 + nb*ns
    vals = np.array(toks[i:i+nk*per], dtype=float).reshape(nk, per)
    kx = vals[:, 0]; E = vals[:, 1:].reshape(nk, ns, nb)
    i += nk*per
    nlab = int(toks[i]); i += 1
    ticks, labels = [], []
    for j in range(nlab):
        ticks.append(float(toks[i])); lab = toks[i+1].strip("'")
        labels.append(r"$\Gamma$" if lab.lower() in ("gamma", "\\gamma") else lab)
        i += 2
    return ef, kx, E[:, 0, :], ticks, labels


title, out, ref = sys.argv[1], sys.argv[2], sys.argv[3]
methods = sys.argv[4:]
ef0, kx, E0, ticks, labels = read_bands(ref)
E0 -= ef0

fig, ax = plt.subplots(figsize=(6, 6))
# baseline: black, thick, behind
for b in range(E0.shape[1]):
    ax.plot(kx, E0[:, b], color="black", lw=2.4, alpha=0.85,
            label="converged (baseline)" if b == 0 else None, zorder=1)

colors = {"new": "#d62728", "raw": "#1f77b4", "proj": "#7f7f7f"}
styles = {"new": "-", "raw": "--", "proj": ":"}
for m in methods:
    label, f = m.split("::")
    key = "new" if "new" in label.lower() or "filt" in label.lower() else \
          ("proj" if "proj" in label.lower() else "raw")
    ef, kx2, E, _, _ = read_bands(f)
    E -= ef
    for b in range(E.shape[1]):
        ax.plot(kx2, E[:, b], color=colors[key], lw=1.3, ls=styles[key],
                alpha=0.9, label=label if b == 0 else None, zorder=2)

ax.axhline(0, color="0.5", lw=0.7, ls="-")
for t in ticks:
    ax.axvline(t, color="0.8", lw=0.7)
ax.set_xticks(ticks); ax.set_xticklabels(labels)
ax.set_xlim(kx.min(), kx.max()); ax.set_ylim(-7, 7)
ax.set_ylabel(r"$E - E_F$ (eV)"); ax.set_title(title)
ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
fig.tight_layout(); fig.savefig(out, dpi=150)
print(f"wrote {out}")
