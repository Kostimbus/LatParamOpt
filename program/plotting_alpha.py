from estimator import *

Logging.set_level(Logging.LEVEL0)

# params = LWE.Parameters(n=23, q=3371, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(4))
# params_bigger = LWE.Parameters(n=40, q=3371, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(4))
# output = LWE.primal_usvp(params)
# output_bigger = LWE.primal_usvp(params_bigger, red_shape_model="gsa")
# #output_all = LWE.estimate(params)
# params_bigger_q = LWE.Parameters(n=40, q=7477, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(4))
# output_bigger_q = LWE.primal_usvp(params_bigger_q)

# params_diff_noise = LWE.Parameters(n=40, q=7477, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(8))
# output_diff_noise = LWE.primal_usvp(params_diff_noise)

# output_all_bigger_q = LWE.estimate(params_bigger_q)
# output_all_diff_noise = LWE.estimate(params_diff_noise)

# print("\nn=23, usvp\n", output, "\n")

# print("\nn=40, q=3371, Xe=Bin(4)\n", output_bigger, "\n")

# print("\nn=40, q=7477, Xe=Bin(4)\n", output_bigger, "\n")

# print("\nn=40, q=7477, Xe=Bin(8)\n", output_diff_noise, "\n")
# #print(LWE.estimate(params_diff_noise))


# print("\n\n", output_all_bigger_q, "\n\n")
# print("\n\n", output_all_diff_noise, "\n\n")


# print("\n\n\n TEST \n")
# params_test = LWE.Parameters(n=40, q=257, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(10))
# output_test = LWE.primal_usvp(params_test, red_shape_model="gsa")
# params_test_2 = LWE.Parameters(n=40, q=7477, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(10))
# output_test_2 = LWE.primal_usvp(params_test_2, red_shape_model="gsa")
# print(output_test, "\n\n\n", output_test_2)




# params_list = {
#     "n=50, q=3371, 4": LWE.Parameters(n=50, q=3371, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(4)),
#     "n=50, q=100, 4": LWE.Parameters(n=50, q=100, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(4)),
#     "n=50, q=5407, 4": LWE.Parameters(n=50, q=5407, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(4)),
#     "n=60, q=5407, 4": LWE.Parameters(n=60, q=5407, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(4)),
#     "n=50, q=5500, 4": LWE.Parameters(n=50, q=5500, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(4)),
#     "n=50, q=9533, 4": LWE.Parameters(n=50, q=9533, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(4)),
#     "n=50, q=9533, 2": LWE.Parameters(n=50, q=9533, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(2)),
#     "n=50, q=9533, 16": LWE.Parameters(n=50, q=9533, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(16)),
#     "n=100, q=9533, 16": LWE.Parameters(n=100, q=9533, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(16)),
#     "n=50, q=95330, 16": LWE.Parameters(n=50, q=95330, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(16))
# }

# for name, params in params_list.items():
#     print(f"\n Estimating security for: {name}")
#     output_all = LWE.estimate(params)



# import matplotlib.pyplot as plt
# import numpy as np
# import seaborn as sns
# import re

# # Define ranges for parameters
# n_values = range(100, 401, 100)         # LWE dimension
# q_values = range(1000, 5001, 1000)      # Modulus
# t_values = [2, 4, 6]                    # CenteredBinomial parameters

# # Initialize result matrices for Dual and Dual-Hybrid attacks
# dual_security = np.zeros((len(n_values), len(q_values)))
# hybrid_security = np.zeros((len(n_values), len(q_values)))

# # Perform estimations
# for i, n in enumerate(n_values):
#     for j, q in enumerate(q_values):
#         for t in t_values:
#             print(f"Estimating for n={n}, q={q}, t={t}...")  # Progress tracking
            
#             # CenteredBinomial standard deviation
#             sigma = (t / 2) ** 0.5
#             alpha = sigma / q
            
#             # Set up LWE parameters with CenteredBinomial
#             params = LWE.Parameters(n=n, q=q, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(t))
#             output_all = LWE.estimate(params)
            
#             # Extracting security estimates for Dual and Dual-Hybrid
#             dual_str = output_all['dual'].str()
#             hybrid_str = output_all['dual_hybrid'].str()
            
#             # Parse ROP values
#             dual_rop = float(re.search(r'rop: ≈2\^([0-9.]+)', dual_str).group(1))
#             hybrid_rop = float(re.search(r'rop: ≈2\^([0-9.]+)', hybrid_str).group(1))
            
#             # Store results
#             dual_security[i, j] = dual_rop
#             hybrid_security[i, j] = hybrid_rop

# # Plot heatmaps
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
# sns.heatmap(dual_security, annot=True, fmt=".2f", cmap='coolwarm', xticklabels=q_values, yticklabels=n_values, ax=ax1)
# ax1.set_title('Dual Attack Security (bits)')
# ax1.set_xlabel('Modulus (q)')
# ax1.set_ylabel('Dimension (n)')

# sns.heatmap(hybrid_security, annot=True, fmt=".2f", cmap='coolwarm', xticklabels=q_values, yticklabels=n_values, ax=ax2)
# ax2.set_title('Dual-Hybrid Attack Security (bits)')
# ax2.set_xlabel('Modulus (q)')
# ax2.set_ylabel('Dimension (n)')

# plt.tight_layout()
# plt.show()




# import matplotlib.pyplot as plt
# import numpy as np
# import re
# from mpl_toolkits.mplot3d import Axes3D

# # Define ranges for parameters
# n_values = [100, 200, 300, 400]         # LWE dimensions
# q_values = [1000, 2000, 3000, 4000, 5000]  # Moduli
# t_values = [2, 4, 6]                    # CenteredBinomial parameters
# attack_types = ["arora-gb", "bkw", "usvp", "bdd", "dual", "dual_hybrid"]

# # Prepare the grid
# X, Y, Z = np.meshgrid(n_values, q_values, t_values, indexing='ij')
# voxel_data = {}

# # Initialize data storage for each attack type
# for attack in attack_types:
#     voxel_data[attack] = np.zeros(X.shape)

# # Perform estimations for each attack type
# for i, n in enumerate(n_values):
#     for j, q in enumerate(q_values):
#         for k, t in enumerate(t_values):
#             print(f"Estimating for n={n}, q={q}, t={t}...")  # Progress tracking

#             # Set up LWE parameters with CenteredBinomial
#             params = LWE.Parameters(n=n, q=q, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(t))
#             output_all = LWE.estimate(params)

#             # For each attack type, extract security estimation
#             for attack in attack_types:
#                 attack_str = output_all[attack].str()
#                 rop = float(re.search(r'rop: ≈2\^([0-9.]+)', attack_str).group(1))
#                 voxel_data[attack][i, j, k] = rop

# # --- Plotting ---
# fig = plt.figure(figsize=(18, 10))
# for idx, attack in enumerate(attack_types):
#     ax = fig.add_subplot(2, 3, idx + 1, projection='3d')
#     ax.set_title(f'{attack} Attack')
#     ax.set_xlabel('n (Dimension)')
#     ax.set_ylabel('q (Modulus)')
#     ax.set_zlabel('t (Binomial Parameter)')
    
#     # Voxel plotting
#     colors = voxel_data[attack] / voxel_data[attack].max()  # Normalize colors
#     ax.scatter(X.flatten(), Y.flatten(), Z.flatten(), c=colors.flatten(), cmap='coolwarm', marker='o')
#     plt.tight_layout()

# plt.show()




import matplotlib.pyplot as plt
import numpy as np
import re
import pickle
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go

# Define ranges for parameters
n_values = [100, 200, 300, 400, 500]         # LWE dimensions
q_values = [1000, 2000, 3000, 4000, 5000, 6000]  # Moduli
t_values = [2, 4, 6]                    # CenteredBinomial parameters
attack_types = ["arora-gb", "bkw", "usvp", "bdd", "dual", "dual_hybrid"]

# Prepare the grid
X, Y, Z = np.meshgrid(n_values, q_values, t_values, indexing='ij')
voxel_data = {}

# Initialize data storage for each attack type
for attack in attack_types:
    voxel_data[attack] = np.zeros(X.shape)

# Perform estimations for each attack type
for i, n in enumerate(n_values):
    for j, q in enumerate(q_values):
        for k, t in enumerate(t_values):
            print(f"Estimating for n={n}, q={q}, t={t}...")  # Progress tracking

            # Set up LWE parameters with CenteredBinomial
            params = LWE.Parameters(n=n, q=q, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(t))
            output_all = LWE.estimate(params)

            # For each attack type, extract security estimation
            for attack in attack_types:
                attack_str = output_all[attack].str()
                rop = float(re.search(r'rop: ≈2\^([0-9.]+)', attack_str).group(1))
                voxel_data[attack][i, j, k] = rop

# Save the computed voxel data for later use
with open('voxel_data.pkl', 'wb') as f:
    pickle.dump(voxel_data, f)
print("Data successfully saved to voxel_data.pkl")

# --- Plotting with Matplotlib ---
fig = plt.figure(figsize=(18, 10))
for idx, attack in enumerate(attack_types):
    ax = fig.add_subplot(2, 3, idx + 1, projection='3d')
    ax.set_title(f'{attack} Attack')
    ax.set_xlabel('n (Dimension)')
    ax.set_ylabel('q (Modulus)')
    ax.set_zlabel('t (Binomial Parameter)')
    
    # Voxel plotting with bigger dots
    colors = voxel_data[attack] / voxel_data[attack].max()  # Normalize colors
    ax.scatter(X.flatten(), Y.flatten(), Z.flatten(), c=colors.flatten(), cmap='coolwarm', marker='o', s=60)
    plt.tight_layout()
plt.show()

# --- Interactive Plotly Visualization ---
for attack in attack_types:
    print(f"Saving interactive plot for {attack}...")
    fig = go.Figure(data=[go.Scatter3d(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        mode='markers',
        marker=dict(
            size=6,  # Increased size for Plotly
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
            zaxis_title='t (Binomial Parameter)'
        )
    )
    
    # Save interactive HTML
    html_filename = f"{attack}_interactive.html"
    fig.write_html(html_filename)
    print(f"{html_filename} saved successfully.")
