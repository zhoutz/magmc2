"""Draw FD11 Figure 7 from data/output.h5; no smoothing or intensity rescaling."""

import argparse
from pathlib import Path
import warnings

import h5py
import numpy as np


ENERGY_EDGES = np.geomspace(0.04, 40.0, 121)
MU_EDGES = np.linspace(-1.0, 1.0, 65)
GROUPS = ("Total", "Unscattered", "Scattered")
FORMAT_COLUMNS = {"fd11_fig7_photons_v1": 7, "fd11_photons_v1": 9}
PHOTON_FIELDS = ("id", "energy_inf", "muk", "Q", "U", "V", "n_scattered")
PROJECT_ROOT = Path(__file__).resolve().parents[1]


def observables(stokes):
    """Input is (..., 4, n_mu, n_energy), with photon-number Stokes sums."""
    count, q, u, _ = np.moveaxis(stokes, -3, 0)
    linear = np.hypot(q, u)
    fraction = np.full_like(count, np.nan, dtype=float)
    np.divide(100 * linear, count, out=fraction, where=count > 0)
    angle = np.mod(np.degrees(0.5 * np.arctan2(u, q)), 180.0)
    angle = np.where((count > 0) & (linear > 1e-12 * count), angle, np.nan)
    return fraction, angle


def accumulate(stokes, rows):
    """Bin the seven photon fields, optionally followed by two mode labels."""
    if rows.ndim != 2 or rows.shape[1] not in (7, 9):
        raise ValueError("Expected seven or nine columns per photon")
    if not np.isfinite(rows).all():
        raise ValueError("Non-finite photon data")
    if np.any(rows[:, 1] <= 0) or np.any(np.abs(rows[:, 2]) > 1):
        raise ValueError("Invalid energy or direction cosine")
    norm = np.linalg.norm(rows[:, 3:6], axis=1)
    if np.any(np.abs(norm - 1) > 1e-7):
        raise ValueError("Invalid single-photon Stokes norm")
    for column in (0, 6):
        if np.any(rows[:, column] < 0) or np.any(rows[:, column] != np.floor(rows[:, column])):
            raise ValueError("Photon IDs and scattering counts must be nonnegative integers")
    if rows.shape[1] == 9:
        if not np.isin(rows[:, 7], (0, 1)).all():
            raise ValueError("Invalid initial photon mode")
        valid_last_mode = np.where(rows[:, 6] == 0, rows[:, 8] == -1,
                                   np.isin(rows[:, 8], (0, 1)))
        if not valid_last_mode.all():
            raise ValueError("Last scattering mode disagrees with scattering count")
    for group, select in enumerate((np.ones(len(rows), dtype=bool), rows[:, 6] == 0, rows[:, 6] > 0)):
        part = rows[select]
        for component in range(4):
            weights = None if component == 0 else part[:, component + 2]
            stokes[group, component] += np.histogram2d(
                part[:, 2], part[:, 1], bins=(MU_EDGES, ENERGY_EDGES), weights=weights
            )[0]


def read_hdf5_photons(paths, chunk_size=50000):
    """Read the seed_* compound datasets written by OutputHDF5, in chunks."""
    stokes = np.zeros((3, 4, len(MU_EDGES) - 1, len(ENERGY_EDGES) - 1))
    metadata = []
    seen_paths, seen_seeds = set(), set()
    for path in paths:
        path = Path(path).resolve()
        if path in seen_paths:
            raise ValueError(f"Duplicate input: {path}")
        seen_paths.add(path)
        seeds = []
        totals = dict(emitted=0, escaped=0, absorbed=0, failed=0)
        with h5py.File(path, "r") as stream:
            names = sorted(name for name in stream if name.startswith("seed_"))
            if not names:
                raise ValueError(f"No seed_* photon datasets in {path}")
            for name in names:
                dataset = stream[name]
                location = f"{path}:{name}"
                if (not isinstance(dataset, h5py.Dataset) or dataset.ndim != 1
                        or not set(PHOTON_FIELDS).issubset(dataset.dtype.names or ())):
                    raise ValueError(f"Invalid photon dataset in {location}")
                required = ("seed", "n_escaped", "n_absorbed", "n_failed")
                if any(key not in dataset.attrs for key in required):
                    raise ValueError(f"Missing photon accounting attributes in {location}")
                seed, escaped, absorbed, failed = (int(dataset.attrs[key]) for key in required)
                if name != f"seed_{seed}":
                    raise ValueError(f"Dataset name disagrees with seed in {location}")
                if seed in seen_seeds:
                    raise ValueError(f"Input files reuse random seed {seed}")
                seen_seeds.add(seed)
                seeds.append(seed)
                if min(escaped, absorbed, failed) < 0 or escaped != len(dataset):
                    raise ValueError(f"Photon accounting mismatch in {location}")
                emitted = escaped + absorbed + failed
                last_id = -1
                for start in range(0, len(dataset), chunk_size):
                    photons = dataset[start:start + chunk_size]
                    # energy_inf is already in keV; muk is cos(theta) at infinity.
                    rows = np.column_stack([photons[field] for field in PHOTON_FIELDS])
                    if (rows[0, 0] <= last_id or np.any(np.diff(rows[:, 0]) <= 0)
                            or rows[-1, 0] >= emitted):
                        raise ValueError(f"Invalid photon IDs in {location}")
                    accumulate(stokes, rows)
                    last_id = rows[-1, 0]
                for key, value in zip(totals, (emitted, escaped, absorbed, failed)):
                    totals[key] += value
        # The file has no run-wide completion marker or model settings. These
        # totals describe only the stored seed groups, including failed photons.
        metadata.append(dict(
            path=str(path), seed=f"{min(seeds)}..{max(seeds)} ({len(seeds)} groups)",
            complete="true", rows_read=totals["escaped"],
            stored_groups=len(seeds), **totals,
        ))
        if totals["failed"]:
            warnings.warn(f"{path}: {totals['failed']:,} failed photons have no records",
                          stacklevel=2)
    if not np.allclose(stokes[0], stokes[1] + stokes[2], rtol=1e-11, atol=1e-9):
        raise ValueError("Total Stokes sums disagree with the two components")
    return stokes, metadata


def read_photons(paths, chunk_size=50000):
    """Read HDF5 output, retaining support for the original text records."""
    is_hdf5 = [h5py.is_hdf5(path) for path in paths]
    if all(is_hdf5):
        return read_hdf5_photons(paths, chunk_size)
    if any(is_hdf5):
        raise ValueError("HDF5 and text inputs must be plotted separately")
    return read_text_photons(paths, chunk_size)


def read_text_photons(paths, chunk_size=50000):
    stokes = np.zeros((3, 4, 64, 120))
    metadata = []
    seen_paths, seen_seeds = set(), set()
    for path in paths:
        path = Path(path).resolve()
        if path in seen_paths:
            raise ValueError(f"Duplicate input: {path}")
        seen_paths.add(path)
        info, buffer = {}, []
        rows_read, last_id = 0, -1

        def consume():
            nonlocal rows_read, last_id
            rows = np.loadtxt(buffer, ndmin=2)
            expected = FORMAT_COLUMNS.get(info.get("format"))
            if expected is None:
                raise ValueError(f"Missing or unsupported photon format in {path}")
            if rows.shape[1] != expected:
                raise ValueError(f"Expected {expected} columns per photon in {path}")
            if rows[0, 0] <= last_id or np.any(np.diff(rows[:, 0]) <= 0):
                raise ValueError(f"Photon IDs must increase strictly in {path}")
            accumulate(stokes, rows)
            last_id = rows[-1, 0]
            rows_read += len(rows)
            buffer.clear()

        with path.open() as stream:
            for line in stream:
                if line.startswith("#"):
                    if " = " in line:
                        key, value = line[1:].strip().split(" = ", 1)
                        if key in info:
                            raise ValueError(f"Duplicate metadata key {key} in {path}")
                        info[key] = value
                elif line.strip():
                    buffer.append(line)
                    if len(buffer) == chunk_size:
                        consume()
            if buffer:
                consume()
        if info.get("format") not in FORMAT_COLUMNS:
            raise ValueError(f"Missing or unsupported photon format in {path}")
        requested = int(info["requested"])
        if requested < 0 or last_id >= requested:
            raise ValueError(f"Photon accounting mismatch in {path}")
        info["rows_read"] = rows_read
        if info.get("complete") == "true":
            emitted, escaped, absorbed = (int(info[k]) for k in ("emitted", "escaped", "absorbed"))
            if not (requested == emitted == escaped + absorbed and escaped == rows_read
                    and last_id < emitted and absorbed >= 0):
                raise ValueError(f"Photon accounting mismatch in {path}")
        else:
            info["complete"] = "false"
            warnings.warn(
                f"{path}: missing completion marker; plotting {rows_read:,} recorded escaped "
                f"photons from an incomplete sample ({requested:,} requested). "
                "Actual emitted/absorbed totals are unknown.",
                stacklevel=2,
            )
        if info["seed"] in seen_seeds:
            raise ValueError("Input files reuse a random seed; these are not independent samples")
        seen_seeds.add(info["seed"])
        if metadata:
            excluded = {"seed", "requested", "emitted", "escaped", "absorbed",
                        "scatterings", "path", "complete", "rows_read"}
            previous = {k: v for k, v in metadata[0].items() if k not in excluded}
            current = {k: v for k, v in info.items() if k not in excluded}
            if previous != current:
                raise ValueError("Input files have different model settings or tolerances")
        info["path"] = str(path)
        metadata.append(info)
    if not np.allclose(stokes[0], stokes[1] + stokes[2], rtol=1e-11, atol=1e-9):
        raise ValueError("Total Stokes sums disagree with the two components")
    return stokes, metadata


def save_bins(prefix, stokes, metadata):
    fraction, angle = observables(stokes)
    rows = []
    for group in range(3):
        for j in range(64):
            for i in range(120):
                rows.append((group, ENERGY_EDGES[i], ENERGY_EDGES[i + 1],
                             MU_EDGES[j], MU_EDGES[j + 1], *stokes[group, :, j, i],
                             fraction[group, j, i], angle[group, j, i]))
    header = "\n".join([
        "FD11 Figure 7: raw photon counts and Stokes sums; no rescaling or smoothing",
        "groups: 0=Total 1=Unscattered 2=Scattered; NaN means undefined",
        *(f"input = {item['path']} (seed {item['seed']}, complete {item['complete']}, "
          f"recorded escaped {item['rows_read']}, "
          f"emitted {item['emitted'] if item['complete'] == 'true' else 'unknown'}, "
          f"failed {item.get('failed', 0)})"
          for item in metadata),
        "group energy_lo_keV energy_hi_keV mu_lo mu_hi N sum_Q sum_U sum_V linear_percent angle_deg",
    ])
    np.savetxt(str(prefix) + "_bins.txt", rows, fmt=["%d"] + ["%.12e"] * 4 + ["%d"] + ["%.12e"] * 5, header=header)


def draw(prefix, stokes, metadata):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.ticker import MultipleLocator

    fraction, angle = observables(stokes)
    counts = stokes[:, 0]
    log_counts = np.full_like(counts, np.nan)
    np.log10(counts, out=log_counts, where=counts > 0)
    max_log = max(1.0, np.log10(max(1.0, counts.max())))
    norms = (Normalize(0, max_log), Normalize(0, 100), Normalize(90, 125))
    labels = (r"$\log_{10} N$", "Pol. Fraction (%)", "Pol. Angle (deg)")
    # Match the paper's purple/blue/cyan/green/yellow/red scale, including
    # bright red at 100% rather than jet's dark-maroon endpoint.
    cmap = LinearSegmentedColormap.from_list("fd11_rainbow", [
        (0.0, "#8000ff"), (0.125, "#0000ff"), (0.25, "#00ffff"),
        (0.5, "#00ff00"), (0.75, "#ffff00"), (1.0, "#ff0000"),
    ])
    cmap.set_bad("black")
    plt.rcParams.update({"font.family": "serif", "font.size": 12, "mathtext.fontset": "dejavuserif"})
    fig, axes = plt.subplots(3, 3, figsize=(10.5, 8.5), sharex=True, sharey=True,
                             gridspec_kw={"hspace": 0.06, "wspace": 0.06})
    fig.subplots_adjust(left=0.09, right=0.86, bottom=0.09, top=0.88)
    for row, values in enumerate((log_counts, fraction, angle)):
        for col in range(3):
            ax = axes[row, col]
            mesh = ax.pcolormesh(np.log10(ENERGY_EDGES), MU_EDGES,
                                 np.ma.masked_invalid(values[col]), cmap=cmap,
                                 norm=norms[row], shading="flat", rasterized=True)
            ax.text(0.94, 0.94, f"({chr(97 + row * 3 + col)})", ha="right", va="top",
                    transform=ax.transAxes, color="black",
                    bbox={"facecolor": "white", "alpha": 0.5, "edgecolor": "none", "pad": 1})
            ax.set_ylim(-1, 1)
            ax.set_xlim(np.log10(0.04), np.log10(40))
            ax.xaxis.set_major_locator(MultipleLocator(0.5))
            ax.xaxis.set_minor_locator(MultipleLocator(0.1))
            ax.yaxis.set_major_locator(MultipleLocator(0.5))
            ax.tick_params(which="both", direction="in", top=True, right=True)
            if row == 0:
                ax.set_title(GROUPS[col], fontsize=15)
            if row == 2:
                ax.set_xlabel(r"$\log_{10}(E/\mathrm{keV})$")
            if col == 0:
                ax.set_ylabel(r"$\cos\theta$")
        box = axes[row, 2].get_position()
        cax = fig.add_axes([0.875, box.y0, 0.023, box.height])
        bar = fig.colorbar(mesh, cax=cax, extend="both" if row == 2 else "neither")
        bar.set_label(labels[row])
        if row == 1:
            bar.set_ticks([0, 25, 50, 75, 100])
        elif row == 2:
            bar.set_ticks([90, 98, 107, 116, 125])
    if all(item["complete"] == "true" for item in metadata):
        emitted = sum(int(item["emitted"]) for item in metadata)
        sample = f"{emitted:,} emitted photons"
        if all("stored_groups" in item for item in metadata):
            groups = sum(item["stored_groups"] for item in metadata)
            sample += f" in {groups:,} stored seed groups"
            failed = sum(item["failed"] for item in metadata)
            if failed:
                sample += f" ({failed:,} failed)"
    else:
        recorded = sum(item["rows_read"] for item in metadata)
        sample = f"INCOMPLETE SAMPLE · {recorded:,} recorded escaped photons"
    model = metadata[0].get("table1_model")
    title = "FD11 Figure 7" + (f" — {model}" if model else "")
    fig.suptitle(f"{title}\n{sample}", fontsize=13)
    fig.savefig(str(prefix) + ".png", dpi=240)
    fig.savefig(str(prefix) + ".pdf", dpi=240)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="*", type=Path,
                        default=[PROJECT_ROOT / "data/output.h5"],
                        help="Input HDF5 files (default: project data/output.h5)")
    parser.add_argument("--output", type=Path, default=PROJECT_ROOT / "data/fd11_fig7",
                        help="Output prefix (default: project data/fd11_fig7)")
    args = parser.parse_args()
    stokes, metadata = read_photons(args.inputs)
    if stokes[0, 0].sum() == 0:
        raise ValueError("No photons in the plotted energy range")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    save_bins(args.output, stokes, metadata)
    draw(args.output, stokes, metadata)
    if all(m["complete"] == "true" for m in metadata):
        print(f"Emitted: {sum(int(m['emitted']) for m in metadata):,}")
    else:
        print("Emitted: unknown (incomplete sample)")
    print(f"Recorded escaped: {sum(m['rows_read'] for m in metadata):,}")
    if all("failed" in m for m in metadata):
        print(f"Absorbed: {sum(m['absorbed'] for m in metadata):,}")
        print(f"Failed: {sum(m['failed'] for m in metadata):,}")
    for group, name in enumerate(GROUPS):
        print(f"{name} in 0.04–40 keV: {int(stokes[group, 0].sum()):,}")
    print(f"Saved {args.output}.png, {args.output}.pdf, {args.output}_bins.txt")


if __name__ == "__main__":
    main()
