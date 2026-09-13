#!/usr/bin/env python3
"""Neutral BayPass POD, an explicitly identified NumPy port of the v3.1 model.

This is NOT an execution of R simulate.baypass, and is NOT bitwise equivalent
to its R RNG. It uses NumPy PCG64. The statistical model is transcribed from
the hash-verified official utility: sample complete observed N rows; draw
pi~Beta(a,b); alpha=pi+sqrt(pi*(1-pi))*L*z, L L'=Omega; clip alpha to [0,1];
Y~Binomial(N,alpha). Then apply global MAC and match observed global MAF bins,
using exactly the same ascertainment design as this project's 06_make_pod.R.

No environmental effect is simulated. Outputs are null calibration data.
The supplied Omega and Beta posterior means are fixed plug-in estimates;
posterior uncertainty in those parameters is not integrated out.
"""

import argparse
import collections
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import sys

import numpy as np

OFFICIAL_UTILS_SHA256 = "a52262f648aeb30ea97330e197cade08b0b835878172823ee50c7823c12102ce"
BREAKS = np.array([0, .01, .02, .03, .05, .075, .10, .15, .20, .30, .40, .5000001])


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def file_info(path):
    p = Path(path).resolve(strict=True)
    return {"path": str(p), "sha256": sha256(p), "bytes": p.stat().st_size}


def read_beta(path):
    lines = Path(path).read_text().splitlines()
    header = lines[0].split()
    canonical = [re.sub(r"[^a-z0-9]", "", h.lower()) for h in header]
    if canonical.count("mean") != 1 or sum(h in ("param", "parameter") for h in canonical) != 1:
        raise ValueError("Beta parameter table must have PARAM and Mean columns")
    ip = next(i for i, h in enumerate(canonical) if h in ("param", "parameter"))
    im = canonical.index("mean")
    values = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        row = line.split()
        key = re.sub(r"[^a-z0-9]", "", row[ip].lower())
        if key in values:
            raise ValueError("Duplicated Beta parameter")
        values[key] = float(row[im])
    if set(values) != {"abetapi", "bbetapi"}:
        raise ValueError("Expected exactly a_beta_pi and b_beta_pi")
    a, b = values["abetapi"], values["bbetapi"]
    if not np.isfinite([a, b]).all() or min(a, b) <= 0:
        raise ValueError("Beta parameters must be finite and positive")
    return a, b


def write_tsv(path, rows):
    rows = list(rows)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def largest_remainder(target, counts):
    raw = target * np.asarray(counts, dtype=float) / np.sum(counts)
    out = np.floor(raw).astype(int)
    # R order(..., decreasing=TRUE) preserves original order for tied remainders.
    order = np.argsort(-(raw - out), kind="stable")
    out[order[:target - out.sum()]] += 1
    return out


def simulate(rng, n, observed_n, cholesky, a, b):
    """Port of individual-genotyping neutral model, pi.maf=0, no fixed filter."""
    source_rows = rng.integers(0, len(observed_n), size=n)
    nn = observed_n[source_rows]
    # Official R utility special-cases a=b=1 as runif().
    pi = rng.uniform(size=n) if a == 1 and b == 1 else rng.beta(a, b, size=n)
    z = rng.standard_normal((n, observed_n.shape[1]))
    # NumPy cholesky is LOWER triangular. Row-vector form is z @ L.T.
    alpha = pi[:, None] + np.sqrt(pi * (1 - pi))[:, None] * (z @ cholesky.T)
    clipped = np.clip(alpha, 0., 1.)
    yy = rng.binomial(nn, clipped)
    return yy, nn, source_rows, pi, alpha


def self_test():
    """Moment checks independently target covariance orientation and sampling."""
    rng = np.random.Generator(np.random.PCG64(928132))
    omega = np.array([[.04, .015, -.007], [.015, .09, .006], [-.007, .006, .03]])
    lower = np.linalg.cholesky(omega)
    z = rng.standard_normal((300000, 3)) @ lower.T
    cov = np.cov(z, rowvar=False, ddof=0)
    cov_error = float(np.max(np.abs(cov - omega)))
    assert cov_error < .0015, "Cholesky covariance mismatch"
    pi0 = .3
    alpha = pi0 + np.sqrt(pi0 * (1 - pi0)) * z
    assert np.max(np.abs(alpha.mean(axis=0) - pi0)) < .001
    assert np.max(np.abs(np.cov(alpha, rowvar=False, ddof=0) - pi0 * (1 - pi0) * omega)) < .0004
    a, b = 2.4, 1.3
    beta = rng.beta(a, b, 300000)
    beta_mean_error = float(abs(beta.mean() - a / (a + b)))
    assert beta_mean_error < .002
    assert abs(beta.var() - a * b / ((a + b)**2 * (a + b + 1))) < .0005
    draws = rng.binomial(18, .27, 300000)
    assert abs(draws.mean() - 18 * .27) < .015
    assert abs(draws.var() - 18 * .27 * .73) < .035
    nn0 = np.array([[0, 2, 20], [4, 0, 8]])
    yy, nn, idx, pi, raw = simulate(rng, 10000, nn0, lower, .4, .8)
    assert np.array_equal(nn, nn0[idx])
    assert (yy[nn == 0] == 0).all()
    assert ((yy >= 0) & (yy <= nn)).all()
    assert np.array_equal(np.clip(np.array([-1., 0., .4, 1., 2.]), 0, 1), [0, 0, .4, 1, 1])
    assert np.array_equal(largest_remainder(7, [1, 1, 1]), [3, 2, 2])
    # In-memory deterministic replay proves exact PCG64 reproducibility.
    rep1 = simulate(np.random.Generator(np.random.PCG64(61)), 100, nn0, lower, a, b)
    rep2 = simulate(np.random.Generator(np.random.PCG64(61)), 100, nn0, lower, a, b)
    assert all(np.array_equal(x, y) for x, y in zip(rep1, rep2))
    return {"status": "PASS", "moment_seed": 928132, "moment_draws": 300000,
            "max_abs_gaussian_covariance_error": cov_error,
            "beta_mean_error": beta_mean_error,
            "checks": ["Gaussian covariance orientation", "conditional alpha mean/covariance",
                       "Beta mean/variance", "binomial mean/variance", "whole-row N resampling",
                       "zero N", "clipping endpoints", "largest remainder ties", "RNG replay"]}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test", action="store_true")
    for flag in ("counts", "omega", "beta-params", "population-order", "outdir", "prefix"):
        p.add_argument("--" + flag)
    p.add_argument("--official-utils", default="/tmp/pod_official_v31/baypass_public-v3.1/utils/baypass_utils.R")
    p.add_argument("--core-command-manifest")
    p.add_argument("--target", type=int, default=100000)
    p.add_argument("--batch-size", type=int, default=10000)
    p.add_argument("--min-mac", type=int, default=3)
    p.add_argument("--seed", type=int, default=61001)
    opt = p.parse_args()
    checks = self_test()
    if opt.self_test:
        print(json.dumps(checks, indent=2))
        return
    for field in ("counts", "omega", "beta_params", "population_order", "outdir", "prefix"):
        if not getattr(opt, field):
            p.error("--" + field.replace("_", "-") + " is required")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", opt.prefix):
        raise ValueError("Unsafe prefix")
    if min(opt.target, opt.batch_size, opt.min_mac) < 1 or opt.seed < 0:
        raise ValueError("target, batch-size, min-mac must be positive; seed nonnegative")
    sources = {k: file_info(getattr(opt, k)) for k in ("counts", "omega", "beta_params", "population_order", "official_utils")}
    sources["generator_script"] = file_info(__file__)
    if sources["official_utils"]["sha256"] != OFFICIAL_UTILS_SHA256:
        raise ValueError("Official v3.1 utility SHA256 mismatch")
    observed = np.loadtxt(opt.counts, ndmin=2)
    if not np.isfinite(observed).all() or (observed < 0).any() or (observed != np.floor(observed)).any() or observed.shape[1] % 2:
        raise ValueError("Counts must be nonnegative integers with 2 columns per population")
    observed = observed.astype(np.int64)
    yy0, nn0 = observed[:, ::2], observed[:, ::2] + observed[:, 1::2]
    npop = nn0.shape[1]
    with open(opt.population_order) as handle:
        order_rows = list(csv.DictReader(handle, delimiter="\t"))
    population_ids = [row["population_id"] for row in order_rows]
    if len(population_ids) != npop or len(set(population_ids)) != npop or not all(population_ids):
        raise ValueError("Population order incompatible with counts")
    if opt.core_command_manifest:
        sources["core_command_manifest"] = file_info(opt.core_command_manifest)
        with open(opt.core_command_manifest) as handle:
            core_rows = list(csv.DictReader(handle, delimiter="\t"))
        beta_name = Path(opt.beta_params).name
        matches = [row for row in core_rows if Path(row["outprefix"]).name + "_summary_beta_params.out" == beta_name]
        if len(matches) != 1:
            raise ValueError("Beta file not uniquely linked to Core command manifest")
        core = matches[0]
        if core["status"] != "COMPLETED" or core["phase"] != "core" or core["return_code"] != "0":
            raise ValueError("Matched Core run was not completed")
        if core["count_sha256"] != sources["counts"]["sha256"] or core["population_order_sha256"] != sources["population_order"]["sha256"]:
            raise ValueError("Core inputs differ from POD observed inputs")
        if int(core["n_loci"]) != len(observed) or core["population_order"] != ",".join(population_ids):
            raise ValueError("Core population order / count differs")
    else:
        core = None
    omega = np.loadtxt(opt.omega, ndmin=2)
    if omega.shape != (npop, npop) or not np.isfinite(omega).all() or not np.allclose(omega, omega.T, rtol=0, atol=1e-8):
        raise ValueError("Omega must be finite and symmetric with matching populations")
    lower = np.linalg.cholesky(omega)
    a, b = read_beta(opt.beta_params)
    totals0 = nn0.sum(axis=1)
    mac0 = np.minimum(yy0.sum(axis=1), (nn0 - yy0).sum(axis=1))
    if (totals0 <= 0).any() or (mac0 < opt.min_mac).any():
        raise ValueError("Observed rows must pass specified minimum global MAC")
    maf0 = mac0 / totals0
    nb = len(BREAKS) - 1
    obs_counts = np.histogram(maf0, bins=BREAKS)[0]
    targets = largest_remainder(opt.target, obs_counts)
    rng = np.random.Generator(np.random.PCG64(opt.seed))
    ybins, nbins = [[] for _ in range(nb)], [[] for _ in range(nb)]
    filled = np.zeros(nb, dtype=int)
    attempted = 0
    monomorphic = 0
    below_mac = 0
    attempts_by_bin = np.zeros(nb, dtype=int)
    outdir = Path(opt.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = outdir / (opt.prefix + ".numpy_POD_manifest.json")
    if manifest_path.exists() or list(outdir.glob("G." + opt.prefix + "_b*")):
        raise ValueError("Refusing to overwrite existing files for this prefix")
    lock = outdir / (opt.prefix + ".lock")
    with open(lock, "x") as handle:
        handle.write(str(os.getpid()))
    try:
        for iteration in range(1, 501):
            need_total = int((targets - filled).sum())
            if need_total == 0:
                break
            ntry = min(1000000, max(10000, int(np.ceil(2.5 * need_total))))
            yy, nn, _, _, _ = simulate(rng, ntry, nn0, lower, a, b)
            attempted += ntry
            mac = np.minimum(yy.sum(axis=1), (nn - yy).sum(axis=1))
            maf = mac / nn.sum(axis=1)
            eligible = mac >= opt.min_mac
            monomorphic += int(np.sum(mac == 0))
            below_mac += int(np.sum(~eligible))
            bins = np.searchsorted(BREAKS, maf, side="right") - 1
            attempts_by_bin += np.bincount(bins[eligible], minlength=nb)
            for k in range(nb):
                need = int(targets[k] - filled[k])
                if need > 0:
                    take = np.flatnonzero(eligible & (bins == k))[:need]
                    if len(take):
                        ybins[k].append(yy[take])
                        nbins[k].append(nn[take])
                        filled[k] += len(take)
        if not np.array_equal(filled, targets):
            raise RuntimeError("Could not fill MAF bins: " + str(targets - filled))
        yy = np.vstack([v for chunks in ybins for v in chunks])
        nn = np.vstack([v for chunks in nbins for v in chunks])
        shuffle = rng.permutation(opt.target)
        yy, nn = yy[shuffle], nn[shuffle]
        mac = np.minimum(yy.sum(axis=1), (nn - yy).sum(axis=1))
        maf = mac / nn.sum(axis=1)
        actual = np.histogram(maf, bins=BREAKS)[0]
        assert np.array_equal(actual, targets) and (mac >= opt.min_mac).all()
        old_patterns = collections.Counter(map(tuple, nn0.tolist()))
        new_patterns = collections.Counter(map(tuple, nn.tolist()))
        assert set(new_patterns) <= set(old_patterns)
        pattern_rows = [{"N_pattern": ",".join(map(str, k)), "observed_count": old_patterns[k],
                         "POD_count": new_patterns[k], "observed_fraction": old_patterns[k] / len(nn0),
                         "POD_fraction": new_patterns[k] / opt.target}
                        for k in sorted(old_patterns)]
        tv = .5 * sum(abs(row["observed_fraction"] - row["POD_fraction"]) for row in pattern_rows)
        batch_rows = []
        for i, start in enumerate(range(0, opt.target, opt.batch_size), 1):
            end = min(start + opt.batch_size, opt.target)
            geno = np.empty((end - start, 2 * npop), dtype=np.int64)
            geno[:, ::2] = yy[start:end]
            geno[:, 1::2] = nn[start:end] - yy[start:end]
            path = outdir / f"G.{opt.prefix}_b{i:03d}"
            np.savetxt(path, geno, fmt="%d")
            batch_rows.append({"batch": f"b{i:03d}", "count_file": str(path), "n_loci": end-start,
                               "first_null_index": start+1, "last_null_index": end, "sha256": sha256(path)})
        maf_rows = [{"bin": k+1, "lower_inclusive": float(BREAKS[k]), "upper_exclusive": float(BREAKS[k+1]),
                     "observed_count": int(obs_counts[k]), "target_count": int(targets[k]),
                     "POD_count": int(actual[k]), "eligible_attempt_count": int(attempts_by_bin[k])} for k in range(nb)]
        size_rows = []
        for dataset, ns in [("observed", nn0), ("POD", nn)]:
            for j, pop in enumerate(population_ids):
                size_rows.append({"dataset": dataset, "population_id": pop, "min_N": int(ns[:,j].min()),
                                  "mean_N": float(ns[:,j].mean()), "median_N": float(np.median(ns[:,j])),
                                  "max_N": int(ns[:,j].max()), "zero_N_fraction": float(np.mean(ns[:,j] == 0))})
        diagnostics = {"target": opt.target, "attempted": attempted, "monomorphic_attempts": monomorphic,
                       "below_min_mac_attempts": below_mac, "acceptance_fraction": opt.target / attempted,
                       "N_pattern_total_variation_distance": tv, "all_N_patterns_observed": True,
                       "MAF_histogram_exact_target": True, "minimum_retained_MAC": int(mac.min()),
                       "N_pattern_warning": "MAF-bin retention may change N-pattern proportions. Inspect TV and population means; no universal pass threshold.",
                       "moment_and_edge_tests": checks}
        diagnostic_files = []
        for suffix, rows in [("numpy_POD_batches.tsv", batch_rows), ("POD_MAF_bin_matching.tsv", maf_rows),
                             ("POD_N_pattern_summary.tsv", pattern_rows), ("POD_sample_size_summary.tsv", size_rows)]:
            path = outdir / (opt.prefix + "." + suffix)
            write_tsv(path, rows)
            diagnostic_files.append(file_info(path))
        diagnostic_path = outdir / (opt.prefix + ".numpy_POD_diagnostics.json")
        diagnostic_path.write_text(json.dumps(diagnostics, indent=2) + "\n")
        diagnostic_files.append(file_info(diagnostic_path))
        manifest = {"schema": "coral_gea_numpy_neutral_POD_v1", "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                    "status": "COMPLETE", "implementation": "Independent NumPy port of verified official BayPass v3.1 individual-count neutral simulation model",
                    "R_executed": False, "bitwise_equivalent_to_R": False, "rng": "numpy.random.Generator(PCG64)",
                    "seed": opt.seed, "python_version": sys.version, "numpy_version": np.__version__, "platform": platform.platform(),
                    "argv": sys.argv, "sources": sources, "core_link": core, "population_order": population_ids,
                    "beta_prior": {"a": a, "b": b, "mode": "fixed Core posterior means"}, "pi_maf": 0,
                    "omega_eigenvalues": np.linalg.eigvalsh(omega).tolist(), "minimum_global_MAC": opt.min_mac,
                    "MAF_breaks": BREAKS.tolist(), "MAF_bin_targets": targets.tolist(), "target": opt.target,
                    "batches": batch_rows, "diagnostics": diagnostics, "diagnostic_files": diagnostic_files,
                    "limitations": ["Plugin Omega and Beta estimates; parameter uncertainty not integrated", "No environmental effect simulated", "MAF ascertainment can alter sample-size proportions", "R source verified and transcribed, not executed"]}
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        print(json.dumps({"manifest": str(manifest_path), "batches": len(batch_rows), **diagnostics}, indent=2))
    finally:
        lock.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
