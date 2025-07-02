# IQ Fit Puzzle - MPI Parallel Solver

This project implements a **parallel IQ Fit Puzzle solver** using **MPI (Message Passing Interface)** and supports multi-core execution.

---

## 🧱 Project Overview

- **Language:** C++11
- **Parallelization:** MPI (e.g., OpenMPI or MS-MPI)
- **Target Binary:** `iqfit_mpi`
- **Board:** 11x5 IQ Fit board
- **Pieces:** 12 unique pieces with all rotations/flips
- **Goal:** Find all valid puzzle solutions using distributed processing

---

## 🔧 Build Instructions

### 📦 Dependencies

- C++ compiler with MPI support (e.g., `mpic++`)
- MPI runtime (e.g., OpenMPI or MS-MPI)
- Unix-like terminal (macOS/Linux or WSL on Windows)

### ⚙️ Build with Makefile

To compile the solver:

```bash
make
```

This builds the `iqfit_mpi` executable from `iqfit_mpi.cpp`.

---

## 🚀 Run Instructions

### ✅ Run via Makefile Targets

Choose a core count and run using:

```bash
make run1   # Run with 1 core
make run2   # Run with 2 cores
make run4   # Run with 4 cores
make run8   # Run with 8 cores
make run12  # Run with 12 cores
```

📁 Each run stores console output in the `log/` directory:

- Example: `log/run4.txt`

### ▶️ Run Manually

Alternatively, run directly with `mpirun`:

```bash
mpirun -np 4 ./iqfit_mpi
```

Replace `4` with your desired number of MPI processes.

---

## 📂 Output

- All valid solutions are written to `solutions.txt`
- Terminal output (for performance evaluation, progress, etc.) is saved to `log/runX.txt`

---

## 🛉 Clean Up

To remove compiled binaries, logs, and solution files:

```bash
make clean
```

This deletes:

- `iqfit_mpi` binary
- `solutions.txt`
- `log/` folder

---

🌐 GUI Visualization

You can view the first 100 puzzle solutions interactively on the following web interface:

https://the-iq-fit-solver-gui-asur.vercel.app/

## 📌 Notes

- This solver is optimized for multi-core environments and was tested on up to 12 MPI processes.
- Ensure your system allows the number of MPI processes you specify (use `--oversubscribe` on some systems if needed).
