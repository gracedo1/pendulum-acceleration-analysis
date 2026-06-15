"""Analysis script for the smartphone pendulum experiment."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, savgol_filter

# Experimental constants
L = 1.20  # m
G = 9.80  # m/s^2
TOP_SEPARATION_CM = 18.0
ENVELOPE_HOLE_SEPARATION_CM = 18.0
MASS_BASIC_G = 239.0
MASS_ADDED_TOTAL_G = 445.0
CARDBOARD_HEIGHT_CM = 8.5

THEORY_PERIOD = 2 * np.pi * np.sqrt(L / G)

# Stable time ranges selected from each measured signal.
WINDOWS = {
    "angle10_trial1": (18, 45),
    "angle10_trial2": (8, 33),
    "angle10_trial3": (11, 25),
    "angle20_trial1": (10, 55),
    "angle20_trial2": (15, 82),
    "angle20_trial3": (13, 70),
    "angle30_trial1": (15, 130),
    "angle30_trial2": (10, 75),
    "angle30_trial3": (13, 135),
    "cardboard20_trial1": (10, 75),
    "cardboard20_trial2": (8, 73),
    "cardboard20_trial3": (10, 60),
    "mass20_trial1": (10, 75),
    "mass20_trial2": (20, 130),
    "mass20_trial3": (13, 125),
}

CONDITION_LABELS = {
    "angle10": "10 deg",
    "angle20": "20 deg",
    "angle30": "30 deg",
    "cardboard20": "20 deg + cardboard",
    "mass20": "20 deg + added mass",
}

PLOT_ORDER = ["angle10", "angle20", "angle30", "cardboard20", "mass20"]


def load_trial(csv_path: Path) -> pd.DataFrame:
    return pd.read_csv(csv_path)


def estimate_period_from_y(df: pd.DataFrame, start: float, end: float):
    t = df["Time (s)"].to_numpy(float)
    y = df["Linear Acceleration y (m/s^2)"].to_numpy(float)

    mask = (t >= start) & (t <= end)
    tt = t[mask]
    yy = y[mask] - np.mean(y[mask])

    dt = np.median(np.diff(tt))
    win = max(5, int(0.12 / dt))
    win += 1 - (win % 2)
    ys = savgol_filter(yy, win, 3) if win < len(yy) else yy

    n = len(ys)
    freqs = np.fft.rfftfreq(n, dt)
    spectrum = np.abs(np.fft.rfft(ys * np.hanning(n)))
    band = (freqs > 0.2) & (freqs < 1.0)
    T_fft = 1 / freqs[band][np.argmax(spectrum[band])]

    distance = int(0.65 * T_fft / dt)
    prominence = max(np.std(ys) * 0.25, np.max(np.abs(ys)) * 0.06)
    peaks, _ = find_peaks(ys, distance=distance, prominence=prominence)
    troughs, _ = find_peaks(-ys, distance=distance, prominence=prominence)

    periods = []
    if len(peaks) > 2:
        periods.extend(p for p in np.diff(tt[peaks]) if 0.75 * T_fft < p < 1.25 * T_fft)
    if len(troughs) > 2:
        periods.extend(p for p in np.diff(tt[troughs]) if 0.75 * T_fft < p < 1.25 * T_fft)

    zero_crossings = []
    signs = np.signbit(ys)
    for i in np.where(signs[:-1] != signs[1:])[0]:
        if ys[i + 1] != ys[i]:
            z = tt[i] - ys[i] * (tt[i + 1] - tt[i]) / (ys[i + 1] - ys[i])
        else:
            z = tt[i]
        zero_crossings.append(z)

    if len(zero_crossings) > 3:
        half_periods = np.diff(np.array(zero_crossings))
        periods.extend(2 * h for h in half_periods if 0.35 * T_fft < h < 0.65 * T_fft)

    periods = np.array(periods if periods else [T_fft])
    median = np.median(periods)
    periods = periods[(periods > 0.85 * median) & (periods < 1.15 * median)]

    return {
        "T_exp_s": float(np.mean(periods)),
        "period_std_s": float(np.std(periods, ddof=1)) if len(periods) > 1 else np.nan,
        "n_intervals_used": int(len(periods)),
        "T_fft_s": float(T_fft),
    }


def condition_from_trial(trial_name: str) -> str:
    return trial_name.split("_trial")[0]


def analyze_trials(data_dir="data/csv") -> pd.DataFrame:
    rows = []
    for csv_path in sorted(Path(data_dir).glob("*.csv")):
        trial = csv_path.stem
        if trial not in WINDOWS:
            continue
        start, end = WINDOWS[trial]
        result = estimate_period_from_y(load_trial(csv_path), start, end)
        condition = condition_from_trial(trial)
        T_exp = result["T_exp_s"]
        rows.append({
            "trial": trial,
            "condition": condition,
            "condition_label": CONDITION_LABELS.get(condition, condition),
            "analysis_start_s": start,
            "analysis_end_s": end,
            **result,
            "T_theory_s": THEORY_PERIOD,
            "error_rate_percent": abs(T_exp - THEORY_PERIOD) / THEORY_PERIOD * 100,
            "effective_length_m": G * (T_exp / (2 * np.pi)) ** 2,
        })
    return pd.DataFrame(rows)


def summarize(trial_df: pd.DataFrame) -> pd.DataFrame:
    summary = trial_df.groupby(["condition", "condition_label"], as_index=False).agg(
        mean_T_exp_s=("T_exp_s", "mean"),
        sd_between_trials_s=("T_exp_s", "std"),
        mean_error_percent=("error_rate_percent", "mean"),
        mean_effective_length_m=("effective_length_m", "mean"),
        trials=("trial", "count"),
    )
    summary["plot_order"] = summary["condition"].map({c: i for i, c in enumerate(PLOT_ORDER)})
    return summary.sort_values("plot_order").drop(columns="plot_order")


def save_outputs(trial_df: pd.DataFrame, summary: pd.DataFrame, results_dir="results"):
    results_dir = Path(results_dir)
    results_dir.mkdir(exist_ok=True)
    trial_df.to_csv(results_dir / "trial_results.csv", index=False)
    summary.to_csv(results_dir / "group_summary.csv", index=False)

    md = [
        "# Analysis Summary",
        "",
        "## Constants",
        "",
        f"- Pendulum length: {L:.2f} m",
        f"- Gravitational acceleration: {G:.2f} m/s^2",
        f"- Theoretical period: {THEORY_PERIOD:.3f} s",
        f"- Upper support separation: {TOP_SEPARATION_CM:.0f} cm",
        f"- Envelope hole separation: {ENVELOPE_HOLE_SEPARATION_CM:.0f} cm",
        f"- Basic mass: {MASS_BASIC_G:.0f} g",
        f"- Mass after adding weight: {MASS_ADDED_TOTAL_G:.0f} g",
        f"- Smartphone height change with cardboard: {CARDBOARD_HEIGHT_CM:.1f} cm",
        "",
        "## Group summary",
        "",
        summary[["condition_label", "mean_T_exp_s", "sd_between_trials_s", "mean_error_percent", "mean_effective_length_m", "trials"]].round(3).to_markdown(index=False),
        "",
    ]
    (results_dir / "analysis_summary.md").write_text("\n".join(md), encoding="utf-8")


def make_figures(data_dir="data/csv", figures_dir="figures"):
    data_dir = Path(data_dir)
    figures_dir = Path(figures_dir)
    figures_dir.mkdir(exist_ok=True)

    files = sorted(data_dir.glob("*.csv"))
    cols = 3
    rows = int(np.ceil(len(files) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(15, 3 * rows))
    axes = np.array(axes).reshape(-1)
    for ax, csv_path in zip(axes, files):
        df = load_trial(csv_path)
        ax.plot(df["Time (s)"], df["Linear Acceleration y (m/s^2)"], lw=0.5)
        ax.set_title(csv_path.stem, fontsize=9)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("$a_y$ (m/s$^2$)")
        ax.grid(alpha=0.25)
    for ax in axes[len(files):]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(figures_dir / "overview_y_axis.png", dpi=200)
    plt.close(fig)

    rep = data_dir / "angle20_trial2.csv"
    df = load_trial(rep)
    start, end = WINDOWS[rep.stem]
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(df["Time (s)"], df["Linear Acceleration y (m/s^2)"], lw=0.7)
    ax.axvspan(start, end, alpha=0.2)
    ax.set_title("Representative acceleration data: angle20_trial2")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("$a_y$ (m/s$^2$)")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(figures_dir / "representative_angle20_trial2.png", dpi=200)
    plt.close(fig)

    summary = pd.read_csv("results/group_summary.csv")
    fig, ax = plt.subplots(figsize=(9, 4.8))
    x = np.arange(len(summary))
    ax.bar(x, summary["mean_T_exp_s"], yerr=summary["sd_between_trials_s"].fillna(0), capsize=5)
    ax.axhline(THEORY_PERIOD, linestyle="--", label=f"Theory ({THEORY_PERIOD:.3f} s)")
    ax.set_xticks(x)
    ax.set_xticklabels(summary["condition_label"], rotation=20, ha="right")
    ax.set_ylabel("Period (s)")
    ax.set_title("Measured period by condition")
    ax.grid(axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures_dir / "period_by_condition.png", dpi=200)
    plt.close(fig)

    angle_summary = summary[summary["condition"].isin(["angle10", "angle20", "angle30"])]
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.errorbar([10, 20, 30], angle_summary["mean_T_exp_s"], yerr=angle_summary["sd_between_trials_s"], marker="o", capsize=5)
    ax.axhline(THEORY_PERIOD, linestyle="--", label="Small-angle model")
    ax.set_xlabel("Initial angle (degrees)")
    ax.set_ylabel("Period (s)")
    ax.set_title("Period measured at different initial angles")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures_dir / "period_vs_angle.png", dpi=200)
    plt.close(fig)


def main():
    trial_df = analyze_trials()
    summary = summarize(trial_df)
    save_outputs(trial_df, summary)
    make_figures()
    print(f"Theoretical period: {THEORY_PERIOD:.3f} s")
    print(summary[["condition_label", "mean_T_exp_s", "sd_between_trials_s", "mean_error_percent", "mean_effective_length_m", "trials"]].round(3))


if __name__ == "__main__":
    main()
