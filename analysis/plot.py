#!/usr/bin/env python3
"""
plot.py -- figures from a goy run directory.

    python3 analysis/plot.py runs/sharp_inverse_cascade
    python3 analysis/plot.py runs/sharp_inverse_cascade --outdir figs

Reads spectra.csv and scalars.csv and writes:

    spectrum_t<TIME>.png    one panel per checkpoint
    spectrum_overlay.png    all checkpoints together
    invariants.png          energy and enstrophy vs time
    kd_vs_time.png          dissipation wavenumber vs time
    drift.png               invariant drift -- inviscid unforced runs only

The simulation writes raw values and all log-taking and normalisation happens
here. That means a normalisation error can be fixed without re-running anything.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BLUE, RED, GREEN = "#1f4e9c", "#c0392b", "#27ae60"


def load(run_dir):
    run = Path(run_dir)
    spectra = pd.read_csv(run / "spectra.csv")
    scalars = pd.read_csv(run / "scalars.csv")
    params = {}
    manifest = run / "parameters.txt"
    if manifest.exists():
        for line in manifest.read_text().splitlines():
            parts = line.split()
            if len(parts) == 2:
                params[parts[0]] = parts[1]
    return spectra, scalars, params


def describe(params):
    """One-line summary of the run, for figure titles."""
    bits = []
    nu = params.get("nu")
    if nu is not None:
        bits.append(f"nu={float(nu):g}")
    ne = params.get("nu_ekman")
    if ne is not None:
        bits.append(f"nu_ek={float(ne):g}")
    amp = params.get("force_amplitude")
    forced = amp is not None and float(amp) != 0.0
    a, b = params.get("force_shell_a"), params.get("force_shell_b")
    if forced and a is not None and b is not None:
        bits.append(f"forced {a},{b}")
    elif amp is not None and not forced:
        bits.append("unforced")
    es = params.get("ekman_shells")
    if es is not None and int(es) > 0:
        bits.append(f"ekman<{es}")
    integ = params.get("integrator")
    if integ:
        bits.append(integ)
    return "  ".join(bits)


def is_inviscid_unforced(params):
    """True only when energy and enstrophy are exact invariants.

    Defaults are deliberately non-zero: a run with no parameters.txt is
    treated as viscous, so the drift figure is skipped.
    """
    try:
        return (float(params.get("nu", 1.0)) == 0.0
                and float(params.get("nu_ekman", 1.0)) == 0.0
                and float(params.get("force_amplitude", 1.0)) == 0.0)
    except ValueError:
        return False


def reference(ax, x, slope, x0, y0, label, colour):
    ax.plot(x, y0 + slope * (x - x0), "--", color=colour, lw=1.3, label=label)


def plot_spectrum(ax, frame, refs=True):
    valid = frame[(frame.E_norm > 0) & (frame.k_norm > 0)]
    if valid.empty:
        return
    x = np.log10(valid.k_norm.to_numpy())
    y = np.log10(valid.E_norm.to_numpy())
    ax.plot(x, y, "o-", color=BLUE, ms=3.5, lw=1.2, label="E (simulated)")

    if refs and len(x) > 2:
        peak = int(np.argmax(y))
        # Forward enstrophy cascade lives above the forcing scale (BKL, -3);
        # the inverse energy cascade lives below it (Kolmogorov, -5/3).
        reference(ax, x[peak:], -3.0, x[peak], y[peak], "BKL (-3)", RED)
        reference(ax, x[: peak + 1], -5 / 3, x[peak], y[peak],
                  "Kolmogorov (-5/3)", GREEN)

    ax.set_xlabel(r"$\log_{10}(k/k_d)$")
    ax.set_ylabel(r"$\log_{10}(E)$  (normalised)")
    ax.grid(True, alpha=0.3, lw=0.5)


def plot_drift(scalars, out, label, subtitle, params):
    """Relative drift of the inviscid invariants.

    With nu = 0 and no forcing, energy and enstrophy are exact invariants of
    the coupling coefficients, so everything plotted here is numerical error.
    Log axes because the interesting range is ~1e-11 to ~1e-9, invisible on
    the linear scale used by invariants.png.

    Skipped for viscous or forced runs: there the energy change is physical,
    and plotting it as "drift" would be misleading.
    """
    if not is_inviscid_unforced(params):
        return
    if len(scalars) < 2:
        return

    e0 = scalars.energy.iloc[0]
    w0 = scalars.enstrophy.iloc[0]
    if e0 == 0 or w0 == 0:
        return

    t = scalars.time.to_numpy()[1:]
    de = np.abs(scalars.energy.to_numpy()[1:] - e0) / abs(e0)
    dw = np.abs(scalars.enstrophy.to_numpy()[1:] - w0) / abs(w0)

    keep = (de > 0) & (dw > 0) & (t > 0)
    if keep.sum() < 2:
        return

    fig, ax = plt.subplots(figsize=(7, 4.8))
    ax.plot(t[keep], de[keep], "o-", color=BLUE, label="energy")
    ax.plot(t[keep], dw[keep], "o-", color=RED, label="enstrophy")

    # Linear accumulation is a slope-1 line on log-log axes. Each series gets
    # its own guide anchored at its own first point -- a single shared guide
    # would sit at the wrong height for the other series and invite a false
    # reading of its slope.
    t_ref = t[keep]
    ax.plot(t_ref, de[keep][0] * (t_ref / t_ref[0]), "--",
            color=BLUE, lw=1.0, alpha=0.45)
    ax.plot(t_ref, dw[keep][0] * (t_ref / t_ref[0]), "--",
            color=RED, lw=1.0, alpha=0.45,
            label="dashed: linear accumulation")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("simulated time (s)")
    ax.set_ylabel("relative drift from initial value")
    ax.set_title(f"{label} -- invariant drift"
                 + (f"\n{subtitle}" if subtitle else ""))
    ax.grid(True, alpha=0.3, lw=0.5, which="both")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "drift.png", dpi=150)
    plt.close(fig)
    print(f"  drift over the run: energy {de[-1]:.3e}, enstrophy {dw[-1]:.3e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    spectra, scalars, params = load(args.run_dir)
    out = Path(args.outdir or Path(args.run_dir) / "figures")
    out.mkdir(parents=True, exist_ok=True)
    label = params.get("run_name", Path(args.run_dir).name)
    subtitle = describe(params)

    times = sorted(spectra.time.unique())
    print(f"{label}: {len(times)} checkpoints")

    # create a separate panel for each checkpoint
    for t in times:
        frame = spectra[spectra.time == t]
        if frame[(frame.E_norm > 0) & (frame.k_norm > 0)].empty:
            continue
        fig, ax = plt.subplots(figsize=(7, 5))
        plot_spectrum(ax, frame)
        ax.set_title(f"{label} -- t = {t:g} s"
                     + (f"\n{subtitle}" if subtitle else ""))
        ax.legend(frameon=False, fontsize=9)
        fig.tight_layout()
        fig.savefig(out / f"spectrum_t{t:g}.png", dpi=150)
        plt.close(fig)

    # overlay all the energy spectra
    drawn = 0
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    cmap = plt.get_cmap("viridis")
    for i, t in enumerate(times):
        frame = spectra[(spectra.time == t) & (spectra.E_norm > 0)]
        if frame.empty:
            continue
        ax.plot(np.log10(frame.k_norm), np.log10(frame.E_norm), "o-",
                ms=3, lw=1.1, color=cmap(i / max(len(times) - 1, 1)),
                label=f"t = {t:g} s")
        drawn += 1
    final = spectra[(spectra.time == times[-1]) & (spectra.E_norm > 0)]
    if not final.empty:
        x = np.log10(final.k_norm.to_numpy())
        y = np.log10(final.E_norm.to_numpy())
        peak = int(np.argmax(y))
        reference(ax, x[peak:], -3.0, x[peak], y[peak], "BKL (-3)", RED)
        reference(ax, x[: peak + 1], -5 / 3, x[peak], y[peak],
                  "Kolmogorov (-5/3)", GREEN)
    if drawn:
        ax.set_xlabel(r"$\log_{10}(k/k_d)$")
        ax.set_ylabel(r"$\log_{10}(E)$  (normalised)")
        ax.set_title(f"{label} -- spectrum evolution"
                     + (f"\n{subtitle}" if subtitle else ""))
        ax.grid(True, alpha=0.3, lw=0.5)
        ax.legend(frameon=False, fontsize=8, ncol=2)
        fig.tight_layout()
        fig.savefig(out / "spectrum_overlay.png", dpi=150)
    plt.close(fig)

    # inviscid invariants (total energy and enstrophy) plotted against time
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 6.5), sharex=True)
    ax1.plot(scalars.time, scalars.energy, "o-", color=BLUE)
    ax1.set_ylabel("total energy")
    ax1.set_yscale("log")
    ax1.grid(True, alpha=0.3, lw=0.5)
    ax1.set_title(f"{label} -- invariants"
                  + (f"\n{subtitle}" if subtitle else ""))
    ax2.plot(scalars.time, scalars.enstrophy, "o-", color=RED)
    ax2.set_ylabel("total enstrophy")
    ax2.set_yscale("log")
    ax2.set_xlabel("time (s)")
    ax2.grid(True, alpha=0.3, lw=0.5)
    fig.tight_layout()
    fig.savefig(out / "invariants.png", dpi=150)
    plt.close(fig)

    # dissipation wavenumber
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(scalars.time, scalars.k_d, "o-", color=BLUE)
    ax.set_xlabel("time (s)")
    ax.set_ylabel(r"$k_d$")
    ax.set_title(f"{label} -- dissipation wavenumber"
                 + (f"\n{subtitle}" if subtitle else ""))
    ax.grid(True, alpha=0.3, lw=0.5)
    fig.tight_layout()
    fig.savefig(out / "kd_vs_time.png", dpi=150)
    plt.close(fig)

    plot_drift(scalars, out, label, subtitle, params)

    print(f"wrote figures to {out}")


if __name__ == "__main__":
    main()