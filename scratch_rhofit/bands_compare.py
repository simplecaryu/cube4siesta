"""Compare SIESTA .bands files to a reference, aligned at each run's Fermi level.
  argv: ref.bands  label1 file1.bands  label2 file2.bands ...
Reports RMS band deviation (eV) over all bands and within +-6 eV of Ef."""
import sys
import numpy as np


def read_bands(path):
    toks = open(path).read().split()
    ef = float(toks[0])
    # toks[1..4]=kmin,kmax,emin,emax ; toks[5,6,7]=nbands,nspin,nk
    i = 5
    nb, ns, nk = int(toks[i]), int(toks[i+1]), int(toks[i+2]); i += 3
    per = 1 + nb * ns
    vals = np.array(toks[i:i + nk*per], dtype=float).reshape(nk, per)
    kx = vals[:, 0]
    E = vals[:, 1:].reshape(nk, ns, nb)
    return ef, kx, E  # E in eV (absolute)


def main():
    ref = sys.argv[1]
    ef0, kx0, E0 = read_bands(ref)
    E0r = E0 - ef0
    win = np.abs(E0r) <= 6.0     # chemically relevant window (per element)
    print(f"reference {ref}: Ef={ef0:.3f}  bands {E0.shape}")
    print(f"{'method':<30s} {'RMS_all(eV)':>12s} {'RMS_<6eV(eV)':>13s} {'maxdev(eV)':>11s}")
    args = sys.argv[2:]
    for j in range(0, len(args), 2):
        label, f = args[j], args[j+1]
        ef, kx, E = read_bands(f)
        Er = E - ef
        if Er.shape != E0r.shape:
            print(f"{label:<30s}  SHAPE {Er.shape} vs {E0r.shape}"); continue
        d = Er - E0r
        rms_all = np.sqrt(np.mean(d**2))
        rms_win = np.sqrt(np.mean(d[win]**2))
        print(f"{label:<30s} {rms_all:12.4f} {rms_win:13.4f} {np.abs(d).max():11.4f}")


main()
