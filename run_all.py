"""Reproduce every claim the article makes, from scratch.

    python run_all.py            everything
    python run_all.py --quick    skip the two slow fits, keep the proof

Each check below names the sentence in the article it is testing. If a check
fails, the article is wrong and I want to know.

The proof half runs in exact rational and symbolic arithmetic - no floating
point anywhere in it - which is why the article is allowed to say "cannot"
rather than "did not".
"""

import argparse
import json
import statistics
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = []


def run(script):
    print(f"\n$ python {script}")
    proc = subprocess.run([sys.executable, str(HERE / script)],
                          cwd=HERE, capture_output=True, text=True,
                          encoding="utf-8", errors="replace")
    if proc.returncode != 0:
        print(proc.stdout[-2000:])
        print(proc.stderr[-2000:])
        raise SystemExit(f"{script} failed")
    print(proc.stdout.strip()[-700:])
    return proc.stdout


def metrics(name):
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def check(claim, ok, detail):
    RESULTS.append((claim, ok, detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {claim}: {detail}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="skip the two slow fits, keep the proof")
    args = ap.parse_args()

    # ---- one snapshot: the test that everything passes -------------------
    if args.quick:
        print("(skipping the two fits: --quick)")
    else:
        run("engine.py")
        run("single_moment_fit.py")
        s = metrics("metrics_single_moment_fit.json")
        world = s["stageB"]["max_born_error"]
        controls = {k: v["max_born_error"] for k, v in s["stageC"].items()
                    if isinstance(v, dict) and "max_born_error" in v}
        check("the quantum probability rule fits one snapshot of the clockwork "
              "(article: beautifully, far better than it needed to be)",
              world < 0.01, f"error {world:.6f} against a pass mark of 0.01")
        passed_too = [k for k, v in controls.items() if v < 0.01]
        check("and it fits the controls too, including pure random numbers "
              "(article: a test a random number generator passes is not a test)",
              len(passed_too) == 2 and len(controls) == 2,
              ", ".join(f"{k} {v:.6f} passed too" for k, v in controls.items()))

        # ---- two snapshots: the wheel wins -------------------------------
        run("two_moment_quantum.py")
        run("the_rotor.py")
        r = metrics("metrics_the_rotor.json")
        pairs = [(p["D_cls"], p["D_B_quantum"])
                 for f in r["stageB"]["frames"] for p in f["per_tau"]]
        wins = sum(1 for c, q in pairs if c < q)
        check("twenty two-moment comparisons, and the plain wheel wins every one "
              "(article: twenty out of twenty)",
              len(pairs) == 20 and wins == 20, f"{wins} of {len(pairs)}")
        med_c = statistics.median(c for c, _ in pairs)
        med_q = statistics.median(q for _, q in pairs)
        check("the wheel gets them right to about a percent, the quantum formula misses badly",
              med_c < 0.02 and med_q > 0.2,
              f"wheel {100 * med_c:.1f}% off, quantum formula {100 * med_q:.0f}% off")

    # ---- the proof -------------------------------------------------------
    run("the_proof.py")
    p = metrics("metrics_the_proof.json")
    lemmas = {k: v for k, v in p["stageA"].items()
              if k.startswith("L") and isinstance(v, dict)}
    check("every lemma of the proof checks out in exact arithmetic "
          "(article: no floating point anywhere in the proof half)",
          len(lemmas) == 5 and all(v["pass"] for v in lemmas.values()),
          f"{sum(1 for v in lemmas.values() if v['pass'])} of {len(lemmas)} lemmas")

    cases = lemmas["L2"]["cases"]
    values = sorted({c["value"] for c in cases})
    check("every possible kind of day is enumerated, and none beats the ceiling "
          "(article: it tallies every possible kind of day)",
          len(cases) == 8 and lemmas["L2"]["pointwise_bound_1"],
          f"{len(cases)} cases, scoring only {values} - nothing above 1")

    l4 = lemmas["L4"]
    check("and a quantum system does beat it, reaching exactly 3/2 at 60 degrees "
          "(article: flip to TOUCH, set the gap to 60 degrees, watch it go to 1.5)",
          l4["K_at_pi_3_equals_3_2"] and l4["is_maximum"]
          and l4["violates_classical_bound_1"],
          f"K3 at 60 degrees = {l4['K_at_pi_3']}, and it is the maximum")

    l5 = lemmas["L5"]
    check("no watcher, of any setting, reproduces what the toucher does "
          "(article: every possible passive setting, exhaustively - none match it)",
          l5["separation_holds"] and not l5["any_passive_matches_reset"],
          f"all {len(l5['passive_a_tested'])} passive settings tested, none match the touch")

    k3 = p["stageC"]
    score = k3["K3_at_delta"]
    check("my own watched clockwork sits exactly on the ceiling, never through it "
          "(article: pressed flat against it, and never once through)",
          0.99 < score <= 1.0 and not k3["K3_at_delta_exceeds_classical_bound"],
          f"score {score:.6f} against the ceiling of 1")

    print("\n" + "=" * 68)
    for claim, ok, _ in RESULTS:
        print(f"[{'PASS' if ok else 'FAIL'}] {claim}")
    print("=" * 68)
    bad = [c for c, ok, _ in RESULTS if not ok]
    if bad:
        raise SystemExit(f"{len(bad)} check(s) FAILED")
    print(f"all {len(RESULTS)} checks reproduced")


if __name__ == "__main__":
    main()
