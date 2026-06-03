import numpy as np
import matplotlib.pyplot as plt


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


def gillespie_disassembly(beta_eps, nu=1.0, rng=None):
    if rng is None:
        rng = np.random.default_rng()

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
            rate_i = nu * np.exp(-beta_eps * m_i)

            candidates.append(i)
            rates.append(rate_i)
            m_values.append(m_i)

        rates = np.array(rates, dtype=float)
        K = rates.sum()

        dt = -np.log(rng.random()) / K
        t += dt

        probs = rates / K
        choice = rng.choice(len(candidates), p=probs)

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
        "time": times,
        "waiting_time": t_wait,
        "n": np.array(n_values),
        "b": np.array(b_values),
        "removed_triangle": np.array(removed_triangles, dtype=object),
        "broken_bonds": np.array(broken_bonds),
    }


def interpolate_step(times, values, grid):
    idx = np.searchsorted(times, grid, side="right") - 1
    idx = np.clip(idx, 0, len(values) - 1)
    return values[idx]


def run_many(beta_eps, n_traj=1000, nu=1.0, n_grid=1000, seed=1):
    rng = np.random.default_rng(seed)

    trajectories = []
    final_times = []

    for _ in range(n_traj):
        traj = gillespie_disassembly(beta_eps, nu=nu, rng=rng)
        trajectories.append(traj)
        final_times.append(traj["time"][-1])

    t_max = np.percentile(final_times, 99)
    grid = np.linspace(0, t_max, n_grid)

    all_n = np.zeros((n_traj, n_grid))
    all_b = np.zeros((n_traj, n_grid))

    for r, traj in enumerate(trajectories):
        all_n[r] = interpolate_step(traj["time"], traj["n"], grid)
        all_b[r] = interpolate_step(traj["time"], traj["b"], grid)

    return {
        "time_grid": grid,
        "mean_n": all_n.mean(axis=0),
        "std_n": all_n.std(axis=0),
        "mean_b": all_b.mean(axis=0),
        "std_b": all_b.std(axis=0),
        "trajectories": trajectories,
    }


if __name__ == "__main__":
    nu = np.exp(0.11269746)
    #for beta_eps in [6.45126437, 6.07169255, 5.58674352, 5.15753681, 4.91393197]:
    for beta_eps in [6.45126437, ]:
        print(f"beta_eps {beta_eps}")
        result = run_many(beta_eps=beta_eps, n_traj=10000, nu=1.0)

        t = result["time_grid"] / nu
        
        np.savetxt(f"n_triangles_{beta_eps:.2f}.dat", np.c_[(t, result["mean_n"])])
        np.savetxt(f"n_bonds_{beta_eps:.2f}.dat", np.c_[(t, result["mean_b"])])

