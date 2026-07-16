#!/usr/bin/env python3
import re
import subprocess
from pathlib import Path
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

# === Config ===
ROOT_METHODS = Path("../single-models") # dir with subdirs for each method, each containing PDB files
REF_DIR      = Path("../pdb/single") # dir with reference PDB files
OUT_DIR      = Path("../csv") # output CSV files
RNAALIGN     = "RNAalign"  # or full path if needed
MAX_WORKERS  = 8

OUT_DIR.mkdir(parents=True, exist_ok=True)

# Standard RNAalign output: prefer TM-score normalized by Chain_2
RE_TM_CHAIN2 = re.compile(
    r"^TM-score=\s*([0-9.]+)\s*\(if normalized by length of Chain_2",
    re.MULTILINE
)

# Any TM-score line
RE_TM_ANY = re.compile(r"^TM-score=\s*([0-9.]+)", re.MULTILINE)


def parse_standard_tmscore(stdout: str) -> float | None:
    """
    Standard -d (default settings):
    Prefer TM-score normalized by Chain_2.
    """
    m = RE_TM_CHAIN2.search(stdout)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            return None

    vals = RE_TM_ANY.findall(stdout)
    if len(vals) >= 2:
        try:
            return float(vals[1])
        except ValueError:
            return None

    return None


def parse_d_tmscore(stdout: str) -> float | None:
    """
    For runs with -d flag:
    take the LAST TM-score= line.
    """
    vals = RE_TM_ANY.findall(stdout)
    if not vals:
        return None
    try:
        return float(vals[-1])
    except ValueError:
        return None


def run_one_alignment(model_pdb: Path, ref_pdb: Path, d_value: int | None):
    """
    Run RNAalign once.
    Returns: (score, error)
    """
    cmd = [RNAALIGN, str(model_pdb), str(ref_pdb)]
    if d_value is not None:
        cmd += ["-d", str(d_value)]

    try:
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    except Exception as e:
        return None, f"exec error: {e}"

    if proc.returncode != 0:
        return None, f"retcode {proc.returncode}: {proc.stderr.strip()}"

    if d_value is None:
        score = parse_standard_tmscore(proc.stdout)
        if score is None:
            return None, "could not parse standard TM-score"
    else:
        score = parse_d_tmscore(proc.stdout)
        if score is None:
            return None, f"could not parse TM-score for -d {d_value}"

    return score, None


def run_rnaalign_all(method: str, pid: str):
    """
    Run standard + -d 1..5 for one PDB ID.
    Returns:
      pid,
      {
        "TM-scoreRNA": ...,
        "TM-scoreRNA-d1": ...,
        ...
      },
      [error strings]
    """
    model_pdb = ROOT_METHODS / method / f"{pid}.pdb"
    ref_pdb   = REF_DIR / f"{pid}.pdb"

    if not model_pdb.is_file():
        return pid, None, [f"missing model {model_pdb}"]
    if not ref_pdb.is_file():
        return pid, None, [f"missing ref {ref_pdb}"]

    row = {"PDB_ID": pid}
    errs = []

    # Standard run
    score, err = run_one_alignment(model_pdb, ref_pdb, d_value=None)
    row["TM-scoreRNA"] = score
    if err:
        errs.append(f"standard: {err}")

    # -d 1..5 runs
    for d in range(1, 6):
        col = f"TM-scoreRNA-d{d}"
        score, err = run_one_alignment(model_pdb, ref_pdb, d_value=d)
        row[col] = score
        if err:
            errs.append(f"-d {d}: {err}")

    # Keep row even if some columns failed, unless everything failed
    score_cols = [k for k in row if k != "PDB_ID"]
    if all(row[c] is None for c in score_cols):
        return pid, None, errs

    return pid, row, errs


def collect_ids(dirpath: Path) -> set[str]:
    return {p.stem for p in dirpath.glob("*.pdb")}


def process_method(method: str):
    method_dir = ROOT_METHODS / method
    if not method_dir.is_dir():
        return

    model_ids = collect_ids(method_dir)
    ref_ids   = collect_ids(REF_DIR)
    common    = sorted(model_ids & ref_ids)

    if not common:
        print(f"[{method}] no common IDs with reference")
        return

    results = []
    errors  = []

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        fut2id = {ex.submit(run_rnaalign_all, method, pid): pid for pid in common}
        for fut in as_completed(fut2id):
            pid, row, errs = fut.result()
            if row is not None:
                results.append(row)
            if errs:
                for err in errs:
                    errors.append((pid, err))

    out_csv = OUT_DIR / f"rnalign-{method}.csv"
    if results:
        df = pd.DataFrame(results)
        cols = [
            "PDB_ID",
            "TM-scoreRNA",
            "TM-scoreRNA-d1",
            "TM-scoreRNA-d2",
            "TM-scoreRNA-d3",
            "TM-scoreRNA-d4",
            "TM-scoreRNA-d5",
        ]
        df = df[cols]
        df.sort_values("PDB_ID", inplace=True)
        df.to_csv(out_csv, index=False)
        print(f"[{method}] wrote {out_csv} with {len(df)} rows.")
    else:
        print(f"[{method}] no successful TM-scores; nothing written.")

    if errors:
        log_path = OUT_DIR / f"rnalign-{method}.errors.txt"
        with open(log_path, "w") as fh:
            for pid, err in errors:
                fh.write(f"{pid}\t{err}\n")
        print(f"[{method}] {len(errors)} failures/issues → {log_path}")


def main():
    methods = [d.name for d in ROOT_METHODS.iterdir() if d.is_dir()]
    if not methods:
        print("No method directories found.")
        return

    for method in sorted(methods):
        process_method(method)


if __name__ == "__main__":
    main()