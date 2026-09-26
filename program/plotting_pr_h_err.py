import numpy as np
import re
import pickle
import os
from itertools import product
from estimator import LWE, ND
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import plotly.graph_objects as go

# === Parameter Grid (same ranges as before) ===
n_values = list(range(50, 1051, 100))
q_values = list(range(300, 10001, 1000))
#t_values = list(range(2, 9))
t_values = list(range(2, 9))

checkpoint_every = 5

# === Output Directory ===
OUTPUT_DIR = "outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Use a dedicated filename to avoid clashing with your multi-attack file
voxel_data_file = os.path.join(OUTPUT_DIR, "voxel_data_primal_hybrid.pkl")

# === Grid & Resume Logic ===
shape = (len(n_values), len(q_values), len(t_values))
param_combinations = list(product(enumerate(n_values), enumerate(q_values), enumerate(t_values)))

def extract_log2_rop_from_str(s: str) -> float:
    """
    Parse 'rop: ≈2^<num>' from estimator's string output.
    Returns float(<num>) or 0.0 if not found.
    """
    m = re.search(r'rop:\s*≈?2\^([0-9.]+)', s)
    return float(m.group(1)) if m else 0.0

# Load or initialize voxel data
if os.path.exists(voxel_data_file):
    with open(voxel_data_file, "rb") as f:
        voxels = pickle.load(f)
    if voxels.shape != shape:
        # In case you tweak ranges, start fresh with the new shape
        voxels = np.zeros(shape, dtype=float)
    print("🔄 Resuming from saved voxel_data_primal_hybrid.pkl")
else:
    voxels = np.zeros(shape, dtype=float)
    print("🆕 Starting new data grid for primal_hybrid")

# Filter remaining work
remaining_jobs = []
for (i, n), (j, q), (k, t) in param_combinations:
    done = voxels[i, j, k] != 0.0
    if not done:
        remaining_jobs.append(((i, n), (j, q), (k, t)))

print(f"💡 Remaining parameter sets to evaluate: {len(remaining_jobs)}")

# === Main Execution Loop ===
processed = 0
for (i, n), (j, q), (k, t) in remaining_jobs:
    try:
        # Keep Xs=SparseTernary(16) and Xe=CenteredBinomial(t) like your first script.
        # If you prefer the example settings, change to SparseTernary(8) / CB(4).
        lwe_params = LWE.Parameters(
            n=n,
            q=q,
            Xs=ND.SparseTernary(16),    # <- change to 8 if you want the example's Xs
            Xe=ND.CenteredBinomial(t)   # <- change to 4 if you want the example's Xe
        )
        out = LWE.primal_hybrid(lwe_params)
        voxels[i, j, k] = extract_log2_rop_from_str(out.str())
    except Exception as e:
        voxels[i, j, k] = 0.0
        print(f"⚠️ Error on n={n}, q={q}, t={t}: {e}")

    processed += 1

    # Periodic Checkpoint
    if processed % checkpoint_every == 0:
        with open(voxel_data_file, "wb") as f:
            pickle.dump(voxels, f)
        with open(os.path.join(OUTPUT_DIR, f"primal_hybrid_progress_{processed}.txt"), "w") as log:
            log.write(f"{processed} evaluations completed.\n")
        print(f"📍 Checkpoint after {processed} evaluations")

# Final Save
with open(voxel_data_file, "wb") as f:
    pickle.dump(voxels, f)
print("✅ All evaluations complete and saved for primal_hybrid.")

# === Visualization ===
X, Y, Z = np.meshgrid(n_values, q_values, t_values, indexing='ij')

# --- Matplotlib 3D Plot ---
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.set_title('Primal Hybrid Attack')
ax.set_xlabel('n (Dimension)')
ax.set_ylabel('q (Modulus)')
ax.set_zlabel('t (CB Param)')

max_val = np.max(voxels)
if max_val > 0:
    norm = voxels / max_val
else:
    norm = np.zeros_like(voxels)

ax.scatter(X.flatten(), Y.flatten(), Z.flatten(), c=norm.flatten(), cmap='coolwarm', marker='o', s=60)
plt.tight_layout()
matplot_path = os.path.join(OUTPUT_DIR, "primal_hybrid_matplotlib_3d_summary.png")
plt.savefig(matplot_path, dpi=150)
print(f"📊 Saved: {matplot_path}")

# --- Plotly HTML ---
print("💾 Saving Plotly plot for primal_hybrid...")
fig = go.Figure(data=[go.Scatter3d(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    mode='markers',
    marker=dict(
        size=6,
        color=voxels.flatten(),
        colorscale='Viridis',
        colorbar=dict(title='log2(rop)'),
        opacity=0.85
    )
)])
fig.update_layout(
    title='primal_hybrid Attack - Interactive',
    scene=dict(
        xaxis_title='n (Dimension)',
        yaxis_title='q (Modulus)',
        zaxis_title='t (CB Param)'
    )
)
html_path = os.path.join(OUTPUT_DIR, "primal_hybrid_interactive.html")
fig.write_html(html_path)
print(f"✅ {html_path} saved.")
