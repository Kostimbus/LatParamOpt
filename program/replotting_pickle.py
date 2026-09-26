import os, pickle, numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import plotly.graph_objects as go

# === CONFIG ===
OUTPUT_DIR = "outputs"
voxel_path = os.path.join(OUTPUT_DIR, "voxel_data_primal_hybrid.pkl")  # current
reference_pickle_paths = [
    os.path.join(OUTPUT_DIR, "voxel_data.pkl"),  # <- your old pickle (multi-attack or single)
    # add more paths if needed
]
attack_key = "primal_hybrid"  # only used if a dict is passed and you want a specific attack

n_values = list(range(50, 1051, 100))
q_values = list(range(300, 10001, 1000))
t_values = list(range(2, 9))

def arrays_from_pickle(p):
    with open(p, "rb") as f:
        obj = pickle.load(f)
    if isinstance(obj, np.ndarray):
        return [obj]
    if isinstance(obj, dict):
        # take all ndarray values (covers multi-attack files)
        return [v for v in obj.values() if isinstance(v, np.ndarray)]
    return []

def compute_global_scale(arrays, ignore_zeros=True):
    parts = []
    for a in arrays:
        a = np.asarray(a)
        mask = np.isfinite(a)
        if ignore_zeros:
            mask &= a > 0
        if np.any(mask):
            parts.append(a[mask].ravel())
    if parts:
        allvals = np.concatenate(parts)
        return float(allvals.min()), float(allvals.max())
    return 0.0, 1.0

# --- load current voxels ---
with open(voxel_path, "rb") as f:
    cur = pickle.load(f)
voxels = cur[attack_key] if isinstance(cur, dict) and attack_key in cur else (cur if isinstance(cur, np.ndarray) else None)
if voxels is None:
    raise RuntimeError("Could not find current voxels array in current pickle.")

# --- load reference arrays (old pickle[s]) ---
ref_arrays = []
for p in reference_pickle_paths:
    if os.path.exists(p):
        try:
            ref_arrays.extend(arrays_from_pickle(p))
        except Exception as e:
            print(f"⚠️ Failed to read reference pickle {p}: {e}")
    else:
        print(f"⚠️ Reference pickle not found: {p}")

# --- compute shared vmin/vmax from current + reference data ---
vmin, vmax = compute_global_scale([voxels] + ref_arrays, ignore_zeros=True)
print(f"Using global color scale vmin={vmin}, vmax={vmax}")

# === Visualization ===
X, Y, Z = np.meshgrid(n_values, q_values, t_values, indexing='ij')

# Matplotlib
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.set_title('primal_hybrid Attack')
ax.set_xlabel('n (Dimension)')
ax.set_ylabel('q (Modulus)')
ax.set_zlabel('t (CB Param)')
sc = ax.scatter(
    X.flatten(), Y.flatten(), Z.flatten(),
    c=voxels.flatten(),
    cmap='viridis',
    vmin=vmin, vmax=vmax,
    marker='o', s=60
)
cb = fig.colorbar(sc, ax=ax, pad=0.02, shrink=0.8)
cb.set_label('Bit Security')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "primal_hybrid_matplotlib_3d_summary.png"), dpi=150)

# Plotly
fig = go.Figure(data=[go.Scatter3d(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    mode='markers',
    marker=dict(
        size=6,
        color=voxels.flatten(),
        colorscale='Viridis',
        cmin=vmin, cmax=vmax,
        colorbar=dict(title='Bit Security'),
        opacity=0.85
    )
)])
fig.update_layout(
    title='primal_hybrid Attack - Interactive',
    scene=dict(xaxis_title='n (Dimension)', yaxis_title='q (Modulus)', zaxis_title='t (CB Param)')
)
fig.write_html(os.path.join(OUTPUT_DIR, "primal_hybrid_interactive.html"))
