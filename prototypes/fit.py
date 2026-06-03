import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution, minimize


# -----------------------------
# Geometry: 20 triangular faces
# -----------------------------

def make_icosahedron_graph():
    faces = [
        (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
        (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
        (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
        (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
    ]

    edge_to_faces = {}
    for f_id, face in enumerate(faces):
        for i in range(3):
            edge = tuple(sorted((face[i], face[(i + 1) % 3])))
            edge_to_faces.setdefault(edge, []).append(f_id)

    neighbors = {i: set() for i in range(20)}
    for attached in edge_to_faces.values():
        if len(attached) == 2:
            i, j = attached
            neighbors[i].add(j)
            neighbors[j].add(i)

    return {i: sorted(v) for i, v in neighbors.items()}


def count_intact_bonds(remaining, neighbors):
    b = 0
    for i in remaining:
        for j in neighbors[i]:
            if j in remaining and j > i:
                b += 1
    return b


# -----------------------------
# Single Gillespie trajectory
# -----------------------------

def gillespie_disassembly_fast(beta_eps, rng=None, neighbors=None):
    if rng is None:
        rng = np.random.default_rng()
    if neighbors is None:
        neighbors = make_icosahedron_graph()

    N = 20
    alive = np.ones(N, dtype=bool)

    # current number of intact neighbors for each triangle
    m = np.array([len(neighbors[i]) for i in range(N)], dtype=np.int64)

    t = 0.0
    times = []
    n_values = []
    b_values = []
    removed_triangles = []
    broken_bonds = []

    n_alive = N
    b_alive = sum(len(neighbors[i]) for i in range(N)) // 2

    while n_alive > 0:
        idx = np.flatnonzero(alive)
        rates = np.exp(-beta_eps * m[idx])
        K = rates.sum()

        if K <= 0:
            break

        t += -np.log(rng.random()) / K

        # faster than rng.choice(..., p=rates/K)
        r = rng.random() * K
        choice_pos = np.searchsorted(np.cumsum(rates), r)
        removed = idx[choice_pos]

        m_removed = m[removed]

        alive[removed] = False
        n_alive -= 1
        b_alive -= m_removed

        # Removing this triangle reduces m by 1 for each still-alive neighbor
        for j in neighbors[removed]:
            if alive[j]:
                m[j] -= 1

        times.append(t)
        n_values.append(n_alive)
        b_values.append(b_alive)
        removed_triangles.append(removed)
        broken_bonds.append(m_removed)

    t_wait = times[0]
    times = np.asarray(times) - t_wait

    return {
        "time": times,
        "waiting_time": t_wait,
        "n": np.asarray(n_values),
        "b": np.asarray(b_values),
        "removed_triangle": np.asarray(removed_triangles),
        "broken_bonds": np.asarray(broken_bonds),
    }

def gillespie_disassembly(beta_eps, rng=None, neighbors=None):
    if rng is None:
        rng = np.random.default_rng()

    if neighbors is None:
        neighbors = make_icosahedron_graph()

    remaining = set(range(20))

    t = 0.0
    times = [0.0]
    n_values = []
    b_values = []
    removed_triangles = [None]
    broken_bonds = [0]

    while len(remaining) > 0:
        candidates = []
        rates = []
        m_values = []

        for i in remaining:
            m_i = sum(j in remaining for j in neighbors[i])
            rate_i = np.exp(-beta_eps * m_i)

            candidates.append(i)
            rates.append(rate_i)
            m_values.append(m_i)

        rates = np.asarray(rates, dtype=float)
        K = rates.sum()

        if K <= 0:
            break

        dt = -np.log(rng.random()) / K
        t += dt

        choice = rng.choice(len(candidates), p=rates / K)

        removed = candidates[choice]
        m_removed = m_values[choice]

        remaining.remove(removed)

        times.append(t)
        n_values.append(len(remaining))
        b_values.append(count_intact_bonds(remaining, neighbors))
        removed_triangles.append(removed)
        broken_bonds.append(m_removed)
        
    # we retain the information about the first opening, but then shift the time
    # so that t = 0 coincides with the time at which the first triangle detaches
    t_wait = times[1]
    times = np.array(times[1:]) - t_wait

    return {
        "time": np.asarray(times),
        "waiting_time": t_wait,
        "n": np.asarray(n_values),
        "b": np.asarray(b_values),
        "removed_triangle": np.asarray(removed_triangles, dtype=object),
        "broken_bonds": np.asarray(broken_bonds),
    }


def interpolate_step(times, values, grid):
    idx = np.searchsorted(times, grid, side="right") - 1
    idx = np.clip(idx, 0, len(values) - 1)
    return values[idx]

# -----------------------------
# Average theoretical curve
# -----------------------------

def simulate_mean_curve(
    beta_eps,
    t_grid,
    nu=1.0,
    n_traj=1000,
    seed=123,
    dt_sim=None
):
    rng = np.random.default_rng(seed)
    neighbors = make_icosahedron_graph()

    all_n = np.zeros((n_traj, len(t_grid)))
    all_b = np.zeros((n_traj, len(t_grid)))

    for r in range(n_traj):
        traj = gillespie_disassembly_fast(
            beta_eps=beta_eps,
            rng=rng,
            neighbors=neighbors,
        )

        all_n[r] = interpolate_step(traj["time"], traj["n"], t_grid)
        all_b[r] = interpolate_step(traj["time"], traj["b"], t_grid)

    return {
        "time": np.asarray(t_grid) / nu,
        "mean_n": all_n.mean(axis=0),
        "std_n": all_n.std(axis=0),
        "mean_b": all_b.mean(axis=0),
        "std_b": all_b.std(axis=0),
    }


def simulate_mean_curve_unit_nu(beta_eps, t_grid_unit, n_traj=1000, seed=123):
    return simulate_mean_curve(
        beta_eps=beta_eps,
        t_grid=t_grid_unit,
        n_traj=n_traj,
        seed=seed,
        dt_sim=None,
    )["mean_b"]

# -----------------------------
# Fitting global nu
# -----------------------------

def global_objective_beta_log_nu(
    params,
    datasets,
    n_traj=1000,
    seed=123,
):
    nu = np.exp(params[0])

    beta_values = params[1:]

    loss = 0.0

    for i, data in enumerate(datasets):
        t_sim = data["time"]
        b_sim = data["b"]

        t_theory = nu * t_sim

        b_th = simulate_mean_curve_unit_nu(
            beta_eps=beta_values[i],
            t_grid_unit=t_theory,
            n_traj=n_traj,
            seed=seed + 1000 * i,
        )

        loss += np.mean((b_sim - b_th) ** 2)

    return loss / len(datasets)


def print_progress(x, conv):
    print(x, conv)


def fit_global_nu_and_betas(
    datasets,
    beta_bounds,
    log_nu_bounds,
    n_traj=1000,
    seed=123,
):
    n_data = len(datasets)

    bounds = [log_nu_bounds]
    bounds += [beta_bounds] * n_data

    res = differential_evolution(
        global_objective_beta_log_nu,
        bounds=bounds,
        args=(datasets,),
        seed=1,
        polish=False,
        callback=print_progress
    )

    res2 = minimize(
        global_objective_beta_log_nu,
        res.x,
        args=(datasets,),
        method="Nelder-Mead",
    )

    nu_fit = np.exp(res2.x[0])
    beta_fits = res2.x[1:]
    
    loss = global_objective_beta_log_nu(res.x, datasets)
    
    fits = {}
    for i, data in enumerate(datasets):
        t_grid = data["time"]

        res = simulate_mean_curve(
            beta_eps=beta_fits[i],
            t_grid=t_grid,
            nu=nu_fit,
            n_traj=n_traj,
            seed=seed + 1000 * i,
        )
        
        fits[data["T"]] = np.c_[res["time"], res["mean_b"]]

    return {
        "nu": nu_fit,
        "beta_eps": beta_fits,
        "loss": loss,
        "fits": fits
    }


# -----------------------------
# Plotting best fit
# -----------------------------

def plot_fit(
    datasets,
    fit_result,
    n_traj=1000,
    seed=999,
):
    nu = fit_result["nu"]
    beta_values = fit_result["beta_eps"]

    plt.figure()

    for i, data in enumerate(datasets):
        t_sim = np.asarray(data["time"])
        n_sim = np.asarray(data["b"])
        dt_sim = np.median(np.diff(t_sim)) if len(t_sim) > 1 else None

        theory = simulate_mean_curve(
            beta_eps=beta_values[i],
            nu=nu,
            t_grid=t_sim,
            n_traj=n_traj,
            seed=seed + 1000 * i,
            dt_sim=dt_sim,
        )

        plt.plot(t_sim, n_sim, "o", alpha=0.6)
        plt.plot(
            theory["time"],
            theory["mean_n"],
            "-",
            label=rf"$\beta\epsilon={beta_values[i]:.3f}$, $\nu={nu:.3g}$",
        )

    plt.xlabel("time")
    plt.ylabel(r"$\langle n(t)\rangle$")
    plt.legend()
    plt.tight_layout()
    plt.show()


# -----------------------------
# Example usage
# -----------------------------

if __name__ == "__main__":
    datasets = []
    for T in ["0.085", "0.090", "0.095", "0.100", "0.105"]:
        filename = f"risultati/disassembly_T_{T}.txt"
        times, bs = np.loadtxt(filename, unpack=True)
        datasets.append({"T" : T, "time" : times, "b" : bs})

    fit = fit_global_nu_and_betas(
        datasets,
        beta_bounds=(4.0, 7.0),
        log_nu_bounds=(np.log(1e-4), np.log(1e4)),
        n_traj=1000,
        seed=1,
    )

    print("Fit with beta_eps free:")
    print("nu =", fit["nu"])
    print("beta_eps =", fit["beta_eps"])
    print("loss =", fit["loss"])
    
    for dataset in datasets:
        T = dataset["T"]
        filename = f"fit_disassembly_T_{T}.dat"
        np.savetxt(filename, fit["fits"][T])

    exit(0)
    plot_fit(
        datasets,
        fit,
        n_traj=5000,
        seed=2,
    )
