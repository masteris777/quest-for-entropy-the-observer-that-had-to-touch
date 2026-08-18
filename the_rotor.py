# The rotor: a wheel turning at a steady rate, fitted to the same two-moment data.
# Run:  python the_rotor.py

import sys
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

HERE = Path(__file__).parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
import two_moment_quantum as quantum  # noqa: E402

TAU_GRID = quantum.TAU_GRID
GRID_N = 600
STAGE_C_SEEDS = list(range(200, 220))
BORN_MARGIN = 0.02


# ---------------------------------------------------------------------------
# Circular statistics (lab-goal.md recipe - never np.mean on raw angles)
# ---------------------------------------------------------------------------
def circular_mean_std(angles):
    m = np.angle(np.mean(np.exp(1j * angles)))
    R = np.abs(np.mean(np.exp(1j * angles)))
    R = min(R, 1.0)
    circ_std = float(np.sqrt(-2 * np.log(R))) if R > 0 else float("inf")
    return float(m), circ_std


# ---------------------------------------------------------------------------
# Classical model
# ---------------------------------------------------------------------------
def z_model(u, v, amp, p, q, delta):
    z0 = amp[0] * np.exp(1j * u)
    z1 = amp[1] * np.exp(1j * v)
    z2 = amp[2] * np.exp(1j * (p * u + q * v + delta))
    Z = np.stack([z0, z1, z2], axis=-1)
    norm = np.linalg.norm(Z, axis=-1, keepdims=True)
    return Z / norm


def label_fn(Z, Qf):
    overlaps = np.abs(Z @ Qf.conj()) ** 2
    return np.argmax(overlaps, axis=-1)


def build_C_cls(Qf, amp, p, q, delta, theta0, theta1, tau_grid, grid_n=GRID_N):
    u = np.linspace(0, 2 * np.pi, grid_n, endpoint=False)
    v = np.linspace(0, 2 * np.pi, grid_n, endpoint=False)
    U, V = np.meshgrid(u, v, indexing="ij")
    U = U.ravel()
    V = V.ravel()

    Z_grid = z_model(U, V, amp, p, q, delta)
    label0 = label_fn(Z_grid, Qf)
    counts0 = np.bincount(label0, minlength=3).astype(float)

    C_cls_by_tau = {}
    for tau in tau_grid:
        Uf = np.mod(U + theta0 * tau, 2 * np.pi)
        Vf = np.mod(V + theta1 * tau, 2 * np.pi)
        Z_flow = z_model(Uf, Vf, amp, p, q, delta)
        label1 = label_fn(Z_flow, Qf)
        idx = label0 * 3 + label1
        C = np.bincount(idx, minlength=9).reshape(3, 3).astype(float)
        row_sums = C.sum(axis=1, keepdims=True)
        row_sums_safe = np.where(row_sums == 0, 1.0, row_sums)
        C_cls_by_tau[tau] = C / row_sums_safe
    return C_cls_by_tau


# ---------------------------------------------------------------------------
# Efficient loss(U): diagonalize once, reuse via matrix power for all tau
# ---------------------------------------------------------------------------
def B_U_all_tau(U, Qf, tau_grid):
    lam, V = np.linalg.eig(U)
    Vinv = np.linalg.inv(V)
    out = {}
    for tau in tau_grid:
        A_tau = V @ np.diag(lam ** tau) @ Vinv
        M = Qf.conj().T @ A_tau @ Qf  # M[k, j] = phi_k^dagger A_tau phi_j
        B = (np.abs(M) ** 2).T  # B[j, k]
        out[tau] = B
    return out


def loss_of_U(U, precomputed_frames, tau_grid):
    total = 0.0
    n = 0
    for Qf, C_cls_by_tau in precomputed_frames:
        B_by_tau = B_U_all_tau(U, Qf, tau_grid)
        for tau in tau_grid:
            d = quantum.row_max_tv(C_cls_by_tau[tau], B_by_tau[tau])
            total += d
            n += 1
    return total / n


def params_to_U(params):
    d0, d1, d2, m01r, m01i, m02r, m02i, m12r, m12i = params
    D = np.diag([d0, d1, d2]).astype(complex)
    M = np.zeros((3, 3), dtype=complex)
    M[0, 1] = m01r + 1j * m01i
    M[0, 2] = m02r + 1j * m02i
    M[1, 2] = m12r + 1j * m12i
    H = D + M + M.conj().T
    return expm(1j * H)


def eigenframe_align(U):
    lam, V = np.linalg.eig(U)
    align_vals = []
    for k in range(3):
        v = V[:, k]
        v = v / np.linalg.norm(v)
        align_vals.append(float(np.max(np.abs(v) ** 2)))
    phases = np.sort(np.abs(np.angle(lam)))
    return phases, float(np.mean(align_vals))


# ---------------------------------------------------------------------------
# Stage A
# ---------------------------------------------------------------------------
def run_stage_a():
    y = quantum.build_stream("torus", quantum.T, np.random.default_rng(88))
    pipe = quantum.build_pipeline(y)

    A = pipe["A"]
    procrustes_orthogonality_residual = float(np.max(np.abs(A @ A.T - np.eye(A.shape[0]))))
    ph_lst = pipe["ph_lst"]
    lstsq_dev = (float(np.max(np.abs(ph_lst - np.sort(quantum.UNROUNDED_ANCHOR))))
                 if ph_lst is not None else None)
    reproduction_ok = bool(procrustes_orthogonality_residual < 1e-12
                            and lstsq_dev is not None and lstsq_dev < 1e-5)
    print(f"A1 reproduction: procrustes_orth_resid={procrustes_orthogonality_residual:.3e}, "
          f"lstsq_dev={lstsq_dev}, ok={reproduction_ok}")

    A_C_used = pipe["A_C_used"]
    Z_C = pipe["Z_C"]
    theta = np.angle(np.diag(A_C_used))  # signed, positional
    print(f"theta (signed) = {theta}")

    amp = np.array([float(np.mean(np.abs(Z_C[:, k]))) for k in range(3)])
    amp_cv = np.array([float(np.std(np.abs(Z_C[:, k])) / np.mean(np.abs(Z_C[:, k]))) for k in range(3)])
    print(f"amplitudes={amp}, amp_cv={amp_cv}")

    best_pq = None
    best_resid = np.inf
    for p in range(-2, 3):
        for q in range(-2, 3):
            if p == 0 and q == 0:
                continue
            resid = abs(p * theta[0] + q * theta[1] - theta[2])
            if resid < best_resid:
                best_resid = resid
                best_pq = (p, q)
    p, q = best_pq
    integer_relation_residual = float(best_resid)
    print(f"integer relation: (p,q)={best_pq}, residual={integer_relation_residual:.3e}")

    x_t = np.angle(Z_C[:, 2]) - p * np.angle(Z_C[:, 0]) - q * np.angle(Z_C[:, 1])
    delta, delta_circ_std = circular_mean_std(x_t)
    print(f"delta={delta:.6f}, delta_circ_std={delta_circ_std:.6f}")

    A6 = integer_relation_residual < 1e-5
    A7 = delta_circ_std < 0.2

    pass_a = bool(reproduction_ok and A6 and A7)

    stage_a_metrics = {
        "reproduction_ok": reproduction_ok,
        "integer_relation": {"p": int(p), "q": int(q), "residual": integer_relation_residual},
        "delta": delta, "delta_circ_std": delta_circ_std,
        "amplitudes": [float(a) for a in amp], "amp_cv": [float(a) for a in amp_cv],
        "pass": pass_a,
    }
    model_params = {"amp": amp, "p": p, "q": q, "delta": delta,
                     "theta0": float(theta[0]), "theta1": float(theta[1]), "theta_all": theta}
    return stage_a_metrics, pipe, model_params


# ---------------------------------------------------------------------------
# Stage B
# ---------------------------------------------------------------------------
def run_stage_b(pipe, model_params):
    Z_C = pipe["Z_C"]
    A_C_used = pipe["A_C_used"]
    amp, p, q, delta = model_params["amp"], model_params["p"], model_params["q"], model_params["delta"]
    theta0, theta1 = model_params["theta0"], model_params["theta1"]

    frames_out = []
    all_D_cls = []
    precomputed_frames = []

    for seed in [91, 93]:
        Qf = quantum.haar_frame(seed)
        s = quantum.symbol_stream(Z_C, Qf)

        u_t = np.angle(Z_C[:, 0])
        v_t = np.angle(Z_C[:, 1])
        Z_model_data = z_model(u_t, v_t, amp, p, q, delta)
        label_model = label_fn(Z_model_data, Qf)
        agreement = float(np.mean(label_model == s))
        print(f"frame {seed}: model_symbol_agreement={agreement:.4f}")

        C_cls_by_tau = build_C_cls(Qf, amp, p, q, delta, theta0, theta1, TAU_GRID)
        precomputed_frames.append((Qf, C_cls_by_tau))

        per_tau = []
        for tau in TAU_GRID:
            F, _ = quantum.counted_F(s, tau)
            B_tau = quantum.born_B(A_C_used, Qf, tau)
            D_B_quantum = quantum.row_max_tv(F, B_tau)
            D_cls = quantum.row_max_tv(F, C_cls_by_tau[tau])
            all_D_cls.append(D_cls)
            per_tau.append({"tau": tau, "D_cls": float(D_cls), "D_B_quantum": float(D_B_quantum)})
            print(f"  tau={tau}: D_cls={D_cls:.4f}, D_B_quantum={D_B_quantum:.4f}")

        frames_out.append({"seed": seed, "model_symbol_agreement": agreement, "per_tau": per_tau})

    all_D_cls = np.array(all_D_cls)
    frac_le_005 = float(np.mean(all_D_cls <= 0.05))
    median_D_cls = float(np.median(all_D_cls))
    agreement_ok = all(fr["model_symbol_agreement"] >= 0.95 for fr in frames_out)
    completion_ok = (frac_le_005 >= 0.90) and (median_D_cls <= 0.03)
    pass_b = bool(agreement_ok and completion_ok)

    print(f"Stage B: frac_D_cls_le_0.05={frac_le_005:.3f}, median_D_cls={median_D_cls:.4f}, "
          f"agreement_ok={agreement_ok}, pass={pass_b}")

    stage_b_metrics = {
        "frames": frames_out,
        "frac_cells_D_cls_le_005": frac_le_005,
        "median_D_cls": median_D_cls,
        "pass": pass_b,
    }
    return stage_b_metrics, precomputed_frames


# ---------------------------------------------------------------------------
# Stage C
# ---------------------------------------------------------------------------
def run_stage_c(pipe, model_params):
    Z_C = pipe["Z_C"]
    A_C_used = pipe["A_C_used"]
    theta_all = model_params["theta_all"]

    frames_out = []
    born_fractions = []
    all_frac_nulls_beaten = []
    n_degenerate = 0

    for seed in STAGE_C_SEEDS:
        Qf = quantum.haar_frame(seed)
        s = quantum.symbol_stream(Z_C, Qf)
        occ = quantum.occupancies(s)
        degenerate = bool(np.min(occ) < quantum.DEGENERACY_FLOOR)

        frame_entry = {"seed": seed, "degenerate": degenerate,
                        "occupancies": [float(o) for o in occ],
                        "born_fraction": None, "per_tau_frac_nulls_beaten": []}

        if degenerate:
            n_degenerate += 1
            frames_out.append(frame_entry)
            continue

        S = np.tile(occ, (3, 1))
        null_seed = 5000 + (seed - 200)
        null_rng = np.random.default_rng(null_seed)
        # draw 30 Haar-random W's from the seeded rng (manual, since haar_frame takes a seed not an rng)
        null_Us = []
        for _ in range(30):
            G = null_rng.standard_normal((3, 3)) + 1j * null_rng.standard_normal((3, 3))
            W, _ = np.linalg.qr(G)
            U_null = W @ np.diag(np.exp(1j * theta_all)) @ W.conj().T
            null_Us.append(U_null)

        n_preferred = 0
        for tau in TAU_GRID:
            F, _ = quantum.counted_F(s, tau)
            B_tau = quantum.born_B(A_C_used, Qf, tau)
            D_B = quantum.row_max_tv(F, B_tau)
            D_I = quantum.row_max_tv(F, np.eye(3))
            D_S = quantum.row_max_tv(F, S)
            born_preferred = bool(D_B <= min(D_I, D_S) - BORN_MARGIN)
            if born_preferred:
                n_preferred += 1

            n_beaten = 0
            for U_null in null_Us:
                B_null = quantum.born_B(U_null, Qf, tau)
                D_null = quantum.row_max_tv(F, B_null)
                if D_B < D_null:
                    n_beaten += 1
            frac_beaten = n_beaten / len(null_Us)
            all_frac_nulls_beaten.append(frac_beaten)
            frame_entry["per_tau_frac_nulls_beaten"].append({"tau": tau, "frac_beaten": frac_beaten})

        born_fraction = n_preferred / len(TAU_GRID)
        frame_entry["born_fraction"] = born_fraction
        born_fractions.append(born_fraction)
        frames_out.append(frame_entry)
        print(f"  Stage C frame seed={seed}: occ={occ}, born_fraction={born_fraction:.2f}")

    median_born_fraction = float(np.median(born_fractions)) if born_fractions else None
    frac_frames_ge_030 = float(np.mean([b >= 0.30 for b in born_fractions])) if born_fractions else None
    median_frac_nulls_beaten = float(np.median(all_frac_nulls_beaten)) if all_frac_nulls_beaten else None

    print(f"Stage C summary: n_degenerate={n_degenerate}, median_born_fraction={median_born_fraction}, "
          f"frac_frames_ge_030={frac_frames_ge_030}, median_frac_nulls_beaten={median_frac_nulls_beaten}")

    return {
        "n_degenerate": n_degenerate,
        "median_born_fraction": median_born_fraction,
        "frac_frames_ge_030": frac_frames_ge_030,
        "median_frac_nulls_beaten": median_frac_nulls_beaten,
        "frames": frames_out,
    }


# ---------------------------------------------------------------------------
# Stage D
# ---------------------------------------------------------------------------
def run_stage_d(pipe, model_params, precomputed_frames):
    A_C_used = pipe["A_C_used"]
    theta_all = model_params["theta_all"]

    loss_A_C = loss_of_U(A_C_used, precomputed_frames, TAU_GRID)
    print(f"loss_A_C = {loss_A_C:.5f}")

    # Haar ensemble
    rng_haar = np.random.default_rng(921)
    haar_losses = np.empty(2000)
    haar_Us = []
    for i in range(2000):
        G = rng_haar.standard_normal((3, 3)) + 1j * rng_haar.standard_normal((3, 3))
        Qh, _ = np.linalg.qr(G)
        haar_losses[i] = loss_of_U(Qh, precomputed_frames, TAU_GRID)
        haar_Us.append(Qh)
    print(f"Haar ensemble done. min={haar_losses.min():.5f}, mean={haar_losses.mean():.5f}")

    # Spectrum-matched ensemble
    rng_spec = np.random.default_rng(922)
    spec_losses = np.empty(2000)
    spec_Us = []
    D_true = np.diag(np.exp(1j * theta_all))
    for i in range(2000):
        G = rng_spec.standard_normal((3, 3)) + 1j * rng_spec.standard_normal((3, 3))
        W, _ = np.linalg.qr(G)
        U_spec = W @ D_true @ W.conj().T
        spec_losses[i] = loss_of_U(U_spec, precomputed_frames, TAU_GRID)
        spec_Us.append(U_spec)
    print(f"Spectrum-matched ensemble done. min={spec_losses.min():.5f}, mean={spec_losses.mean():.5f}")

    # Optimizer: 20 random starts, seed 923
    rng_opt = np.random.default_rng(923)
    starts = [rng_opt.standard_normal(9) for _ in range(20)]
    opt_results = []
    def objective(p_):
        return loss_of_U(params_to_U(p_), precomputed_frames, TAU_GRID)

    for i, x0 in enumerate(starts):
        res = minimize(objective, x0, method="Nelder-Mead",
                        options={"maxiter": 4000, "fatol": 1e-5})
        U_final = params_to_U(res.x)
        opt_results.append((float(res.fun), U_final))
        print(f"  optimizer start {i}: loss={res.fun:.5f}, nit={res.nit}")

    opt_losses = np.array([r[0] for r in opt_results])
    opt_Us = [r[1] for r in opt_results]

    # Best found overall (Haar + spec-matched + optimizer)
    candidates = [("haar", haar_losses[i], haar_Us[i]) for i in range(2000)] + \
                 [("spec", spec_losses[i], spec_Us[i]) for i in range(2000)] + \
                 [("opt", opt_losses[i], opt_Us[i]) for i in range(len(opt_Us))]
    best_source, loss_best_found, U_best = min(candidates, key=lambda c: c[1])
    print(f"Best found: source={best_source}, loss={loss_best_found:.5f}")

    haar_percentile_of_A_C = float(np.mean(haar_losses <= loss_A_C) * 100)
    spectrum_matched_percentile_of_A_C = float(np.mean(spec_losses <= loss_A_C) * 100)

    opt_phases, align = eigenframe_align(U_best)
    true_sorted = np.sort(np.abs(theta_all))
    opt_phase_max_dev = float(np.max(np.abs(opt_phases - true_sorted)))

    D1 = bool((loss_A_C - loss_best_found) <= 0.02)
    D2 = bool(haar_percentile_of_A_C <= 5.0)
    D3 = bool(opt_phase_max_dev < 0.05 and align >= 0.9)

    print(f"D1 near_optimal (gap={loss_A_C - loss_best_found:.5f}): {D1}")
    print(f"D2 haar_distinguished (percentile={haar_percentile_of_A_C:.2f}): {D2}")
    print(f"D3 recovery (phase_dev={opt_phase_max_dev:.4f}, align={align:.4f}): {D3}")

    return {
        "loss_A_C": float(loss_A_C), "loss_best_found": float(loss_best_found),
        "haar_percentile_of_A_C": haar_percentile_of_A_C,
        "spectrum_matched_percentile_of_A_C": spectrum_matched_percentile_of_A_C,
        "opt_phase_max_dev": opt_phase_max_dev, "opt_eigenframe_align": align,
        "D1_near_optimal": D1, "D2_haar_distinguished": D2, "D3_recovery": D3,
        "best_source": best_source,
    }


if __name__ == "__main__":
    print("=" * 60)
    print("STAGE A - Pipeline + classical-model gate")
    print("=" * 60)
    stageA, pipe, model_params = run_stage_a()
    print(f"\nStage A: {'PASS' if stageA['pass'] else 'FAIL'}")

    if not stageA["pass"]:
        metrics = {
            "experiment": "the_rotor", "stream": "quest-for-entropy", "kind": "hardened-simulation",
            "serves_node": "born_rule_emergence",
            "stageA": stageA,
            "stageB": {"frames": [], "frac_cells_D_cls_le_005": None, "median_D_cls": None, "pass": False},
            "stageC": {"n_degenerate": None, "median_born_fraction": None,
                       "frac_frames_ge_030": None, "median_frac_nulls_beaten": None, "frames": []},
            "stageD": {"loss_A_C": None, "loss_best_found": None,
                       "haar_percentile_of_A_C": None, "spectrum_matched_percentile_of_A_C": None,
                       "opt_phase_max_dev": None, "opt_eigenframe_align": None,
                       "D1_near_optimal": None, "D2_haar_distinguished": None, "D3_recovery": None},
            "verdict": "VOID",
            "verdict_reason": "Stage A gate failed (reproduction or classical-model freeze). Stopping.",
            "leakage_audit": {
                "wrote_outside_own_folder": False, "reimplemented_engine_functions": False,
                "posthoc_threshold_change": False, "c_cls_computed_from_data_counting": False,
                "model_params_reestimated_after_stageA": False, "optimization_target_switched_to_counts": False,
            },
        }
        with open(HERE / "metrics_the_rotor.json", "w") as f:
            json.dump(metrics, f, indent=2)
        print("\nVerdict: VOID")
        sys.exit(0)

    print("\n" + "=" * 60)
    print("STAGE B - Classical completion")
    print("=" * 60)
    stageB, precomputed_frames = run_stage_b(pipe, model_params)
    print(f"\nStage B: {'PASS' if stageB['pass'] else 'FAIL'}")

    print("\n" + "=" * 60)
    print("STAGE C - Frame ensemble (no gate)")
    print("=" * 60)
    stageC = run_stage_c(pipe, model_params)

    if not stageB["pass"]:
        metrics = {
            "experiment": "the_rotor", "stream": "quest-for-entropy", "kind": "hardened-simulation",
            "serves_node": "born_rule_emergence",
            "stageA": stageA, "stageB": stageB, "stageC": stageC,
            "stageD": {"loss_A_C": None, "loss_best_found": None,
                       "haar_percentile_of_A_C": None, "spectrum_matched_percentile_of_A_C": None,
                       "opt_phase_max_dev": None, "opt_eigenframe_align": None,
                       "D1_near_optimal": None, "D2_haar_distinguished": None, "D3_recovery": None},
            "verdict": "ANOMALY",
            "verdict_reason": (
                f"Stage B failed: frac_cells_D_cls_le_005={stageB['frac_cells_D_cls_le_005']:.3f} "
                f"(need >=0.90), median_D_cls={stageB['median_D_cls']:.4f} (need <=0.03), or a frame's "
                f"model_symbol_agreement was below 0.95. The exact classical rotation model does not fully "
                f"explain the counted conditionals - reported loudly per protocol. Stage D not run."
            ),
            "leakage_audit": {
                "wrote_outside_own_folder": False, "reimplemented_engine_functions": False,
                "posthoc_threshold_change": False, "c_cls_computed_from_data_counting": False,
                "model_params_reestimated_after_stageA": False, "optimization_target_switched_to_counts": False,
            },
        }
        with open(HERE / "metrics_the_rotor.json", "w") as f:
            json.dump(metrics, f, indent=2)
        print("\nVerdict: ANOMALY")
        print(metrics["verdict_reason"])
        sys.exit(0)

    print("\n" + "=" * 60)
    print("STAGE D - Compression test")
    print("=" * 60)
    stageD = run_stage_d(pipe, model_params, precomputed_frames)

    D1, D2, D3 = stageD["D1_near_optimal"], stageD["D2_haar_distinguished"], stageD["D3_recovery"]

    if D1 and D2:
        verdict = "PASS"
        verdict_reason = (
            f"Stages A, B passed. D1 (gap={stageD['loss_A_C']-stageD['loss_best_found']:.4f} <= 0.02) "
            f"and D2 (Haar percentile={stageD['haar_percentile_of_A_C']:.2f} <= 5) both hold. "
            f"{'D3 (strong form) also holds: optimizer rediscovered the dashboard eigenphases/eigenframe.' if D3 else 'D3 (strong form) did not hold.'}"
        )
    elif D1 or D2:
        verdict = "PARTIAL"
        verdict_reason = (
            f"Stages A, B passed. Exactly one of D1 ({D1}), D2 ({D2}) holds - a partial distinction. "
            f"gap={stageD['loss_A_C']-stageD['loss_best_found']:.4f}, "
            f"haar_percentile={stageD['haar_percentile_of_A_C']:.2f}."
        )
    else:
        verdict = "FAIL"
        verdict_reason = (
            f"Stages A, B passed. Neither D1 (gap={stageD['loss_A_C']-stageD['loss_best_found']:.4f}) "
            f"nor D2 (haar_percentile={stageD['haar_percentile_of_A_C']:.2f}) holds - the collapse matrix "
            f"is not a distinguished compression of the classical law; the finite-tau preference dissolves "
            f"into generic smoothness."
        )

    metrics = {
        "experiment": "the_rotor", "stream": "quest-for-entropy", "kind": "hardened-simulation",
        "serves_node": "born_rule_emergence",
        "stageA": stageA, "stageB": stageB, "stageC": stageC, "stageD": stageD,
        "verdict": verdict, "verdict_reason": verdict_reason,
        "leakage_audit": {
            "wrote_outside_own_folder": False, "reimplemented_engine_functions": False,
            "posthoc_threshold_change": False, "c_cls_computed_from_data_counting": False,
            "model_params_reestimated_after_stageA": False, "optimization_target_switched_to_counts": False,
        },
    }
    with open(HERE / "metrics_the_rotor.json", "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"\nVerdict: {verdict}")
    print(f"Reason: {verdict_reason}")
