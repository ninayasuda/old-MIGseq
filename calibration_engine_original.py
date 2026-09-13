#!/usr/bin/env python3
"""Run and calibrate independent NumPy BayPass PODs against observed IS results.

This supplementary workflow leaves scripts 06/07 and their historical outputs
unchanged. All probabilities are conditional on the fitted neutral model.
Clopper-Pearson intervals describe POD Monte Carlo uncertainty, not uncertainty
in Omega, demographic history, genotyping or the environmental measurements.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.stats import beta as beta_distribution

SEEDS = (73001, 73002, 73003)
MCMC = dict(npilot=20, pilotlength=500, burnin=5000, nval=1000, thin=20)
IS = dict(minbeta=-0.3, maxbeta=0.3, nbetagrid=201)
BINARY_SHA256 = "a2b70c12fb4a032a18aefa12773492b2eed96bebbe8535616ba3900a78c95835"


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stamp():
    return datetime.now(timezone.utc).isoformat()


def file_path(value, base=None):
    path = Path(value).expanduser()
    if not path.is_absolute() and base is not None:
        path = base / path
    path = path.resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty file: {path}")
    return path


def checked(path, expected, base=None):
    path = file_path(path, base)
    if sha(path) != expected:
        raise ValueError(f"SHA-256 mismatch: {path}")
    return path


def atomic_json(path, obj):
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")
    os.replace(tmp, path)


def atomic_tsv(path, rows, headers=None):
    if headers is None:
        headers = list(rows[0])
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    with tmp.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, headers, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def tsv(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def population_order(path):
    rows = tsv(path)
    values = [r["population_id"] for r in rows]
    if not values or len(set(values)) != len(values):
        raise ValueError("Population order has missing or duplicate population IDs")
    return values


def numeric_matrix(path):
    values = np.loadtxt(path, ndmin=2)
    if not np.isfinite(values).all():
        raise ValueError(f"Non-finite numeric value in {path}")
    return values


def beta_params(path):
    with path.open() as stream:
        lines = [line.split() for line in stream if line.strip()]
    headers = [re.sub("[^a-z0-9]", "", x.lower()) for x in lines[0]]
    p = headers.index("param") if "param" in headers else headers.index("parameter")
    mean = headers.index("mean")
    params = {re.sub("[^a-z0-9]", "", row[p].lower()): float(row[mean]) for row in lines[1:]}
    a, b = params["abetapi"], params["bbetapi"]
    if not math.isfinite(a + b) or min(a, b) <= 0:
        raise ValueError("Pi-Beta prior requires two positive finite parameters")
    return a, b


def read_is(path, expected_n=None):
    with path.open() as stream:
        lines = [line.split() for line in stream if line.strip()]
    headers = lines[0]
    required = ["COVARIABLE", "MRK", "BF(dB)", "Beta_is"]
    indices = [headers.index(x) for x in required]
    output = {}
    for line in lines[1:]:
        cov, mrk = int(line[indices[0]]), int(line[indices[1]])
        bf, coefficient = float(line[indices[2]]), float(line[indices[3]])
        if not math.isfinite(bf + coefficient) or mrk in output.setdefault(cov, {}):
            raise ValueError(f"Invalid or duplicate IS result: {path}")
        output[cov][mrk] = (bf, coefficient)
    if not output:
        raise ValueError(f"No IS results: {path}")
    for cov, rows in output.items():
        n = expected_n if expected_n is not None else len(rows)
        if sorted(rows) != list(range(1, n + 1)):
            raise ValueError(f"Nonconsecutive or incomplete MRK rows, covariate {cov}: {path}")
    return {cov: np.array([rows[i] for i in sorted(rows)]) for cov, rows in output.items()}


def bh(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p, kind="stable")
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty(len(p))
    result[order] = np.minimum(adjusted, 1)
    return result


def cp(k, n, alpha):
    k = np.asarray(k)
    lower = np.where(k == 0, 0., beta_distribution.ppf(alpha / 2, np.maximum(k, 1), n - k + 1))
    upper = np.where(k == n, 1., beta_distribution.ppf(1 - alpha / 2, k + 1, np.maximum(n - k, 1)))
    return lower, upper


def option(command, key):
    if command.count(key) != 1:
        raise ValueError(f"Observed command requires exactly one {key}")
    return command[command.index(key) + 1]


def infer_settings(args, manifest, source_paths):
    order = manifest["population_order"]
    a, b = beta_params(source_paths["beta_params"])
    prior = manifest["beta_prior"]
    if not np.allclose([a, b], [prior["a"], prior["b"]], rtol=0, atol=1e-10):
        raise ValueError("POD prior does not match its source Beta file")
    source_counts = source_paths["counts"]
    observed_runs = []
    provenance = {}
    if args.observed_manifest:
        path = file_path(args.observed_manifest)
        rows = tsv(path)
        if len(rows) != 3 or len({r["seed"] for r in rows}) != 3:
            raise ValueError("Observed manifest requires three distinct IS seeds")
        rows.sort(key=lambda x: int(x["seed"]))
        envs, scaling_modes = [], []
        for row in rows:
            if row["status"] != "COMPLETED" or int(row["return_code"]) != 0:
                raise ValueError("Observed run was not completed")
            command = json.loads(row["command_json"])
            if any(x in command for x in ("-d0yij", "-auxmodel", "-covmcmc", "-poolsizefile")):
                raise ValueError("Observed model is incompatible with standard allele-count IS")
            if row["population_order"].split(",") != order:
                raise ValueError("Observed and POD population orders differ")
            for key, hashkey, source in [("count_file", "count_sha256", "counts"),
                                         ("omega_file", "omega_sha256", "omega"),
                                         ("beta_prior_source_file", "beta_prior_source_sha256", "beta_params"),
                                         ("population_order_file", "population_order_sha256", "population_order")]:
                current = checked(row[key], row[hashkey], path.parent)
                if sha(current) != sha(source_paths[source]):
                    raise ValueError(f"Observed / POD {source} differ")
            environment = checked(row["environment_file"], row["environment_sha256"], path.parent)
            env_order = checked(row["environment_population_order_file"], row["environment_population_order_sha256"], path.parent)
            if population_order(env_order) != order:
                raise ValueError("Environment / genotype population order differs")
            for key, value in {**MCMC, **IS}.items():
                if float(option(command, "-" + key)) != value:
                    raise ValueError(f"Observed {key} differs from required calibration settings")
            if command.count("-setpibetapar") != 1:
                raise ValueError("Observed prior must be fixed")
            i = command.index("-betapiprior")
            if not np.allclose([float(command[i + 1]), float(command[i + 2])], [a, b], rtol=0, atol=1e-10):
                raise ValueError("Observed prior values differ from POD")
            if sha(file_path(command[0])) != BINARY_SHA256:
                raise ValueError("Observed BayPass binary differs from pinned version")
            for flag, intended in [("-efile", environment), ("-omegafile", source_paths["omega"]), ("-countdatafile", source_counts)]:
                if sha(file_path(option(command, flag))) != sha(intended):
                    raise ValueError(f"Observed command / manifest {flag} disagree")
            envs.append(environment)
            scaling_modes.append("-nocovscaling" in command)
            observed_runs.append(file_path(row["outprefix"] + "_summary_betai_reg.out"))
        if len({sha(x) for x in envs}) != 1 or len(set(scaling_modes)) != 1:
            raise ValueError("Observed IS runs use different environments / scaling")
        env, nocovscaling = envs[0], scaling_modes[0]
        provenance["observed_manifest"] = {"path": str(path), "sha256": sha(path)}
        if args.efile and sha(file_path(args.efile)) != sha(env):
            raise ValueError("Explicit environment differs from observed manifest")
    else:
        if not args.efile or not args.population_order or not args.observed_counts:
            raise ValueError("Explicit mode requires --efile, --population-order and --observed-counts")
        env = file_path(args.efile)
        if population_order(file_path(args.population_order)) != order:
            raise ValueError("Explicit population order differs from POD")
        if sha(file_path(args.observed_counts)) != sha(source_counts):
            raise ValueError("Explicit observed counts differ from POD source")
        nocovscaling = args.nocovscaling
        if args.observed_is:
            observed_runs = [file_path(p) for p in args.observed_is]
        provenance["explicit_settings"] = True
    if args.efile_override:
        env = file_path(args.efile_override)
        nocovscaling = args.nocovscaling or nocovscaling
        if not args.observed_is:
            raise ValueError("--efile-override requires --observed-is with the three matching observed outputs")
        observed_runs = [file_path(p) for p in args.observed_is]
        provenance["environment_override"] = {"path": str(env), "sha256": sha(env),
                                               "nocovscaling": nocovscaling}
    for supplied, source in [(args.omega, "omega"), (args.beta_params, "beta_params")]:
        if supplied and sha(file_path(supplied)) != sha(source_paths[source]):
            raise ValueError(f"Explicit {source} differs from POD source")
    envmat = numeric_matrix(env)
    if envmat.shape[1] != len(order) or np.any(np.std(envmat, axis=1) == 0):
        raise ValueError("Environment must be covariates x populations with nonzero variation")
    if args.covariate_index < 1 or args.covariate_index > envmat.shape[0]:
        raise ValueError("Covariate index is outside environment rows")
    return env, nocovscaling, a, b, observed_runs, provenance


def run_one(command, outprefix, inputs, expected_n):
    completed_file = outprefix.with_name(outprefix.name + ".completed.json")
    record = {"command": command, "inputs": inputs}
    if completed_file.exists():
        completed = json.loads(completed_file.read_text())
        if any(completed.get(k) != v for k, v in record.items()) or completed.get("return_code") != 0:
            raise ValueError(f"Resume command/input mismatch: {completed_file}")
        for output in completed["outputs"].values():
            checked(output["path"], output["sha256"])
        return completed
    existing = list(outprefix.parent.glob(outprefix.name + "*"))
    if existing:
        raise ValueError(f"Partial outputs exist; refusing to overwrite {outprefix}: {existing}")
    lock = outprefix.with_name(outprefix.name + ".running.json")
    with lock.open("x") as stream:
        json.dump({**record, "started_utc": stamp()}, stream, indent=2)
    log = outprefix.with_name(outprefix.name + ".stdout.log")
    print(f"START {outprefix.name}", flush=True)
    with log.open("x") as stream:
        stream.write("COMMAND_JSON " + json.dumps(command) + "\n")
        stream.flush()
        result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                                env={**os.environ, "OMP_NUM_THREADS": command[command.index("-nthreads") + 1], "LC_ALL": "C"})
    if result.returncode:
        raise RuntimeError(f"BayPass failed ({result.returncode}); partial outputs retained: {log}")
    is_file = file_path(str(outprefix) + "_summary_betai_reg.out")
    xtx_file = file_path(str(outprefix) + "_summary_pi_xtx.out")
    read_is(is_file, expected_n)
    outputs = {name: {"path": str(path), "sha256": sha(path), "bytes": path.stat().st_size}
               for name, path in [("is", is_file), ("xtx", xtx_file), ("stdout", log)]}
    completed = {**record, "outputs": outputs, "return_code": 0,
                 "completed_utc": stamp(), "runner_sha256": sha(Path(__file__))}
    atomic_json(completed_file, completed)
    lock.unlink()
    print(f"COMPLETE {outprefix.name}", flush=True)
    return completed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pod-manifest", required=True)
    parser.add_argument("--observed-manifest")
    parser.add_argument("--observed-summary", required=True)
    parser.add_argument("--observed-is", nargs=3)
    parser.add_argument("--observed-counts")
    parser.add_argument("--efile")
    parser.add_argument("--efile-override", help="Use this matching supplementary environment with --observed-is")
    parser.add_argument("--omega")
    parser.add_argument("--beta-params")
    parser.add_argument("--population-order")
    parser.add_argument("--nocovscaling", action="store_true")
    parser.add_argument("--baypass-bin", default="/tmp/baypass_unpacked/bin/g_baypass")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--batches", type=int, default=1)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--nthreads", type=int, default=1)
    parser.add_argument("--covariate-index", type=int, default=1)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.workers <= 8 or args.nthreads < 1 or args.batches < 1:
        parser.error("--workers must be 1..8; nthreads and batches must be positive")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", args.prefix):
        parser.error("Unsafe prefix")
    manifest_path = file_path(args.pod_manifest)
    manifest = json.loads(manifest_path.read_text())
    schema = manifest.get("schema", manifest.get("schema_version", manifest.get("workflow")))
    if schema != "coral_gea_numpy_neutral_POD_v1":
        raise ValueError(f"Unsupported NumPy POD schema: {schema}")
    sources = {key: checked(value["path"], value["sha256"], manifest_path.parent)
               for key, value in manifest["sources"].items()}
    pops = manifest["population_order"]
    if population_order(sources["population_order"]) != pops:
        raise ValueError("POD population order differs from source file")
    observed_counts = numeric_matrix(sources["counts"])
    if observed_counts.shape[1] != 2 * len(pops) or np.any(observed_counts < 0) or np.any(observed_counts != np.floor(observed_counts)):
        raise ValueError("Observed allele-count matrix has invalid shape or values")
    omega = numeric_matrix(sources["omega"])
    if omega.shape != (len(pops), len(pops)) or not np.allclose(omega, omega.T, atol=1e-7):
        raise ValueError("Invalid Omega dimensions / symmetry")
    np.linalg.cholesky(omega)
    env, nocovscaling, a, b, observed_is, obs_provenance = infer_settings(args, manifest, sources)
    binary = checked(args.baypass_bin, BINARY_SHA256)
    observed_summary = file_path(args.observed_summary)
    observed = tsv(observed_summary)
    if [int(r["MRK"]) for r in observed] != list(range(1, len(observed) + 1)) or len(observed) != observed_counts.shape[0]:
        raise ValueError("Observed summary MRK rows differ from POD source counts")
    bf_observed = np.array([float(r["BFdB_median"]) for r in observed])
    if not np.isfinite(bf_observed).all():
        raise ValueError("Observed BF contains nonfinite numbers")
    if observed_is:
        actual = [read_is(path, len(observed))[args.covariate_index][:, 0] for path in observed_is]
        if not np.allclose(np.median(actual, axis=0), bf_observed, atol=1e-7, rtol=0):
            raise ValueError("Observed BF summary differs from actual three-seed IS output")
        obs_provenance["observed_is"] = [{"path": str(p), "sha256": sha(p)} for p in observed_is]
    else:
        obs_provenance["observed_is_verified"] = False
    batches = manifest["batches"]
    if args.batches > len(batches):
        raise ValueError("Requested more batches than the generator manifest contains")
    chosen = batches[:args.batches]
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    jobs = []
    seen_hashes = set()
    seeds_by_batch = {}
    for batch_index, batch in enumerate(chosen):
        batch_seeds = tuple(seed + batch_index * 1000 for seed in SEEDS)
        seeds_by_batch[batch["batch"]] = batch_seeds
        counts = checked(batch["count_file"], batch["sha256"], manifest_path.parent)
        if batch["sha256"] in seen_hashes:
            raise ValueError("Duplicate POD batch content")
        seen_hashes.add(batch["sha256"])
        arr = numeric_matrix(counts)
        if arr.shape != (batch["n_loci"], 2 * len(pops)) or np.any(arr < 0) or np.any(arr != np.floor(arr)):
            raise ValueError("Invalid POD count dimensions / values")
        totals = arr[:, ::2] + arr[:, 1::2]
        obs_patterns = {tuple(row) for row in observed_counts[:, ::2] + observed_counts[:, 1::2]}
        if any(tuple(row) not in obs_patterns for row in totals):
            raise ValueError("POD sample-size pattern absent from observed matrix")
        for seed in batch_seeds:
            prefix = outdir / f"{args.prefix}_{batch['batch']}_s{seed}"
            command = [str(binary), "-countdatafile", str(counts), "-efile", str(env),
                       "-omegafile", str(sources["omega"]), "-outprefix", str(prefix), "-seed", str(seed)]
            for key, value in {**MCMC, "nthreads": args.nthreads}.items():
                command.extend(["-" + key, str(value)])
            command.extend(["-setpibetapar", "-betapiprior", str(a), str(b)])
            for key, value in IS.items():
                command.extend(["-" + key, str(value)])
            if nocovscaling:
                command.append("-nocovscaling")
            inputs = {"binary_sha256": sha(binary), "counts_sha256": sha(counts),
                      "environment_sha256": sha(env), "omega_sha256": sha(sources["omega"]),
                      "beta_params_sha256": sha(sources["beta_params"]),
                      "population_order": pops, "generator_manifest_sha256": sha(manifest_path)}
            jobs.append((batch, seed, command, prefix, inputs))
    if args.validate_only:
        print(json.dumps({"status": "VALIDATED", "n_jobs": len(jobs), "n_null": sum(b["n_loci"] for b in chosen),
                          "n_observed": len(observed), "population_order": pops,
                          "commands": [j[2] for j in jobs]}, indent=2))
        return
    completed = {}
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(run_one, command, prefix, inputs, batch["n_loci"]): (batch["batch"], seed)
                   for batch, seed, command, prefix, inputs in jobs}
        for future in as_completed(futures):
            completed[futures[future]] = future.result()
    null_by_cov = {}
    for batch in chosen:
        batch_seeds = seeds_by_batch[batch["batch"]]
        values = [read_is(Path(completed[(batch["batch"], seed)]["outputs"]["is"]["path"]), batch["n_loci"]) for seed in batch_seeds]
        if any(set(run) != set(values[0]) for run in values[1:]):
            raise ValueError("POD covariate IDs differ between seeds")
        for cov in sorted(values[0]):
            for i in range(batch["n_loci"]):
                bfs = [run[cov][i, 0] for run in values]
                betas = [run[cov][i, 1] for run in values]
                null_by_cov.setdefault(cov, []).append({
                    "batch": batch["batch"], "MRK_within_batch": i + 1,
                    "seed1": batch_seeds[0], "seed2": batch_seeds[1], "seed3": batch_seeds[2],
                    "BFdB_run1": bfs[0], "BFdB_run2": bfs[1], "BFdB_run3": bfs[2],
                    "BFdB_median": float(np.median(bfs)), "Beta_is_median": float(np.median(betas))})
    for cov, null_rows in null_by_cov.items():
        atomic_tsv(outdir / f"{args.prefix}.cov{cov}.null_median3.tsv", null_rows)
    null = np.sort(np.array([r["BFdB_median"] for r in null_by_cov[args.covariate_index]]))
    n, m = len(null), len(observed)
    k = n - np.searchsorted(null, bf_observed, side="left")
    p = (k + 1) / (n + 1)
    adjusted = bh(p)
    lower, upper = cp(k, n, .05)
    sim_lower, sim_upper = cp(k, n, .05 / m)
    lower_bh, upper_bh = bh(sim_lower), bh(sim_upper)
    final_rows = []
    for i, row in enumerate(observed):
        final_rows.append({**row, "POD_calibrated": True, "POD_n": n, "POD_upper_tail_count": int(k[i]),
                          "POD_empirical_p": float(p[i]), "POD_BH_adjusted_p": float(adjusted[i]),
                          "POD_tail_p_CP95_lower": float(lower[i]), "POD_tail_p_CP95_upper": float(upper[i]),
                          "POD_tail_p_simultaneous95_lower": float(sim_lower[i]), "POD_tail_p_simultaneous95_upper": float(sim_upper[i]),
                          "POD_BH_on_simultaneous_lower": float(lower_bh[i]), "POD_BH_on_simultaneous_upper": float(upper_bh[i]),
                          "POD_BH05": bool(adjusted[i] <= .05),
                          "POD_BH05_MC_supported": bool(upper_bh[i] <= .05),
                          "POD_BH05_MC_excluded": bool(lower_bh[i] > .05),
                          "POD_BH05_MC_unresolved": bool(lower_bh[i] <= .05 < upper_bh[i])})
    calibrated = outdir / f"{args.prefix}.cov{args.covariate_index}.calibrated.tsv"
    atomic_tsv(calibrated, final_rows)
    summary = {"schema": "coral_gea_final_numpy_POD_calibration_v1", "completed_utc": stamp(),
               "n_observed": m, "n_null": n, "n_batches": len(chosen), "covariate_index": args.covariate_index,
               "population_order": pops, "Pi_beta": [a, b], "nocovscaling": nocovscaling,
               "null_statistic": "median BF(dB) across three independent seeds per batch",
               "batch_seeds": seeds_by_batch,
               "n_BH05": int(np.sum(adjusted <= .05)), "n_BH05_MC_supported": int(np.sum(upper_bh <= .05)),
               "n_BH05_MC_excluded": int(np.sum(lower_bh > .05)),
               "n_BH05_MC_unresolved": int(np.sum((lower_bh <= .05) & (upper_bh > .05))),
               "minimum_empirical_p": float(np.min(p)), "minimum_BH_adjusted_p": float(np.min(adjusted)),
               "BF_null_quantiles": {str(q): float(np.quantile(null, q, method="inverted_cdf")) for q in (.95,.99,.999,.9999)},
               "calibrated_output": {"path": str(calibrated), "sha256": sha(calibrated)},
               "pod_manifest": {"path": str(manifest_path), "sha256": sha(manifest_path)},
               "observed_summary": {"path": str(observed_summary), "sha256": sha(observed_summary)},
               "runner": {"path": str(Path(__file__).resolve()), "sha256": sha(Path(__file__))},
               **obs_provenance,
               "limitations": ["Plug-in neutral model; Omega/Pi estimation uncertainty is not propagated.",
                               "Monte Carlo intervals treat generated retained POD loci as independent draws; fixed exact MAF quotas produce stratified sampling, so Clopper-Pearson intervals are a conservative diagnostic approximation for the pooled mixture, not an exact design-based guarantee.",
                               "Simultaneous tail intervals use Bonferroni alpha=0.05/m; BH applied to interval endpoints brackets BH-adjusted limiting tail probabilities conditional on interval coverage.",
                               "MC_supported/excluded describe numerical resolution only, not biological truth or guaranteed FDR.",
                               "Extending simulation after inspection is adaptive precision checking; fixed-N CI coverage is not anytime-valid.",
                               "Distinct seed triplets are used across batches; full model uncertainty is larger than the displayed intervals."]}
    atomic_json(outdir / f"{args.prefix}.cov{args.covariate_index}.calibration_summary.json", summary)
    print(json.dumps({k: summary[k] for k in ["n_observed","n_null","n_BH05","n_BH05_MC_supported","n_BH05_MC_excluded","n_BH05_MC_unresolved"]}), flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
