import numpy as np
import re
import pickle
import os
from itertools import product
from estimator import LWE, ND
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go

# === Parameter Grid ===
n_values = list(range(50, 1051, 100))
q_values = list(range(300, 10001, 1000))
t_values = list(range(2, 9))
attack_types = ["arora-gb", "bkw", "usvp", "bdd", "dual", "dual_hybrid"]
checkpoint_every = 5

# === Output Directory ===
OUTPUT_DIR = "outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)
voxel_data_file = os.path.join(OUTPUT_DIR, "voxel_data.pkl")

# === Grid & Resume Logic ===
shape = (len(n_values), len(q_values), len(t_values))
param_combinations = list(product(enumerate(n_values), enumerate(q_values), enumerate(t_values)))

# Load or initialize voxel data
if os.path.exists(voxel_data_file):
    with open(voxel_data_file, "rb") as f:
        voxel_data = pickle.load(f)
    print("🔄 Resuming from saved voxel_data.pkl")
else:
    voxel_data = {atk: np.zeros(shape) for atk in attack_types}
    print("🆕 Starting new data grid")

# Filter remaining work
remaining_jobs = []
for (i, n), (j, q), (k, t) in param_combinations:
    done = all(voxel_data[atk][i, j, k] != 0 for atk in attack_types)
    if not done:
        remaining_jobs.append(((i, n), (j, q), (k, t)))

print(f"💡 Remaining parameter sets to evaluate: {len(remaining_jobs)}")

# === Main Execution Loop ===
processed = 0
for (i, n), (j, q), (k, t) in remaining_jobs:
    try:
        output_all = LWE.estimate(LWE.Parameters(n=n, q=q, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(t)))
        for attack in attack_types:
            if attack in output_all:
                match = re.search(r'rop: ≈2\^([0-9.]+)', output_all[attack].str())
                voxel_data[attack][i, j, k] = float(match.group(1)) if match else 0.0
            else:
                voxel_data[attack][i, j, k] = 0.0
    except Exception as e:
        for attack in attack_types:
            voxel_data[attack][i, j, k] = 0.0
        print(f"⚠️ Error on n={n}, q={q}, t={t}: {e}")

    processed += 1

    # Periodic Checkpoint
    if processed % checkpoint_every == 0:
        with open(voxel_data_file, "wb") as f:
            pickle.dump(voxel_data, f)
        with open(os.path.join(OUTPUT_DIR, f"progress_{processed}.txt"), "w") as log:
            log.write(f"{processed} evaluations completed.\n")
        print(f"📍 Checkpoint after {processed} evaluations")

# Final Save
with open(voxel_data_file, "wb") as f:
    pickle.dump(voxel_data, f)
print("✅ All evaluations complete and saved.")

# === Visualization ===
X, Y, Z = np.meshgrid(n_values, q_values, t_values, indexing='ij')

# --- Matplotlib 3D Plot ---
fig = plt.figure(figsize=(18, 10))
for idx, attack in enumerate(attack_types):
    ax = fig.add_subplot(2, 3, idx + 1, projection='3d')
    ax.set_title(f'{attack} Attack')
    ax.set_xlabel('n (Dimension)')
    ax.set_ylabel('q (Modulus)')
    ax.set_zlabel('t (CB Param)')
    norm = voxel_data[attack] / np.max(voxel_data[attack])
    ax.scatter(X.flatten(), Y.flatten(), Z.flatten(), c=norm.flatten(), cmap='coolwarm', marker='o', s=60)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "matplotlib_3d_summary.png"))
print("📊 Saved: matplotlib_3d_summary.png")

# --- Plotly HTMLs ---
for attack in attack_types:
    print(f"💾 Saving Plotly plot for {attack}...")
    fig = go.Figure(data=[go.Scatter3d(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        mode='markers',
        marker=dict(
            size=6,
            color=voxel_data[attack].flatten(),
            colorscale='Viridis',
            colorbar=dict(title='Bit Security'),
            opacity=0.85
        )
    )])
    fig.update_layout(
        title=f'{attack} Attack - Interactive',
        scene=dict(
            xaxis_title='n (Dimension)',
            yaxis_title='q (Modulus)',
            zaxis_title='t (CB Param)'
        )
    )
    html_path = os.path.join(OUTPUT_DIR, f"{attack}_interactive.html")
    fig.write_html(html_path)
    print(f"✅ {html_path} saved.")
