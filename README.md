# LatParamOpt: Automated Parameter Selection for Lattice-Based Cryptography

## Overview
**LatParamOpt** (Lattice Parameter Optimizer) is an automated command-line tool designed to compute secure and efficient parameters for lattice-based cryptographic schemes. It supports fundamental problems including **LWE, RLWE, MLWE, SIS, RSIS, and MSIS**.

Selecting concrete parameters for lattice-based schemes is a complex task involving trade-offs between security against various quantum/classical attacks and computational efficiency. This tool automates the process by bridging the **Z3 SMT solver** with Albrecht's **Lattice Estimator**. It allows users to declare target bit security and structural constraints, and systematically discovers optimized parameter sets.

This project was developed as a Bachelor's Thesis: *"Computing Secure and Efficient Parameters for Lattice-Based Cryptographic Schemes"* by Konstantin Gubenko.

## Key Features
* **Multi-Problem Support:** Computes parameters for LWE, RLWE, MLWE, SIS, RSIS, and MSIS.
* **Declarative Constraints:** Users can specify complex structural constraints (e.g., modulo conditions, relations between parameters) directly via the CLI.
* **Non-Linear Constraint Handling:** Built-in support for solver-friendly encodings of non-linear constraints like primality checks (`is_prime(q)`), powers of two (`is_power_of_two(n)`), and logarithmic relations (`log2(x)`).
* **Smart Search Strategy:** Utilizes a "windowed search" strategy with adaptive step sizes to efficiently navigate huge parameter spaces without brute-forcing.
* **Optimization & Refinement:** Supports setting optimization goals (Lexicographic, Pareto, Box) to maximize efficiency. Includes a post-processing refinement step to eliminate parameter overkill and shrink parameters closer to the exact target security.
* **Big-Number Mode:** Supports massive domains required by some schemes (e.g., modulus $q$ up to $2^{500}$).
* **Interactive 3D Visualizations:** Includes Python scripts (`plotting_*.py`) that generate interactive 3D HTML plots (using Plotly) to visualize estimated bit security across different parameter dimensions (e.g., dimension $n$, modulus $q$, and error term $t$). These plots help intuitively understand how different attacks (Primal, Dual, BKW, Arora-GB) react to parameter changes.

## Installation Requirements
To run the program, the following software prerequisites must be installed on your system:

1. **Python** (recommended: 3.10) and **git**.
2. **SymPy** (Python package).
3. **Z3 Python bindings** (`z3-solver`).
4. **SageMath** (provides the `sage` command and the `sageall` Python module).
5. **lattice-estimator** (GitHub checkout): Clone the estimator from `https://github.com/malb/lattice-estimator/tree/main`.
6. **Project Files**: The source files and folders from this repository (contents of the `program` folder) must be copied into the root directory of the cloned `lattice-estimator`.

## Quick Start
The main entry point for the tool is `finder.py`. You pass the problem name as the first positional argument, followed by the target bit security and the variables you want to define.

### Example 1: Basic LWE Search
Find parameters for LWE targeting 40-bit security:
```bash
python3 finder.py LWE -bit 40 -n n -q q -t t
```

### Example 2: Adding Constraints and Optimization (MLWE)
Find parameters for MLWE with specific constraints and optimize for the module rank `d`:
```bash
python3 finder.py MLWE -bit 40 -N N -n n -d d -q q -t t -c "n % 2 == 0" -c "d > 1" -c "d < 10" -opt d
```

### Example 3: Big Mode and Refinement
Enable big-number mode for large domains and refine the results to reduce overkill:
```bash
python3 finder.py LWE -bit 50 -n n -q q -t t --enable-big --enable-refine --bit-sec-margin 8
```

## Interactive Plots
The `interactive_lwe_plots` directory contains pre-generated interactive 3D HTML plots and the raw voxel data used to generate them. You can open the `.html` files directly in any modern web browser to explore how the parameter space impacts the hardness of various attacks.

## Scientific Contribution
Foundational hardness results (like worst-case to average-case reductions) prove that solving lattice problems is computationally hard asymptotically. However, they do not tell engineers which exact values to choose for a specific security level in practice. This tool bridges the gap between theoretical cryptography and real-world deployment. By modeling cryptographic design constraints symbolically and leveraging SMT solvers, it provides a reproducible, auditable, and automated path from security requirements to concrete deployment parameters.
