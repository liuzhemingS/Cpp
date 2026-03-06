# Cpp

## Electron physics simulation (relativistic, many-body)

This repository now contains a C++ simulation of an electron plasma cloud using a physically grounded model:

- Relativistic particle motion (`u = gamma * v`) integrated with a **relativistic Boris pusher**.
- Pairwise electron-electron **Coulomb repulsion** (`O(N^2)` direct summation).
- Idealized external fields: uniform magnetic field + quadrupole electric trap.
- Diagnostic output in CSV format for time-evolution analysis.

### Build

```bash
g++ -std=c++17 -O2 -Wall -Wextra -pedantic main.cpp -o electron_sim
```

### Run

```bash
./electron_sim
```

The run writes diagnostics to `electron_diagnostics.csv` with columns:

- `step`
- `time_s`
- `mean_kinetic_eV`
- `total_kinetic_eV`
- `total_pair_potential_eV`
- `mean_radius_um`

## Scientific scope and limits

The simulation is realistic in the sense of classical + relativistic charged-particle dynamics, but (as in most tractable particle simulations) it still makes approximations:

- No full quantum many-body wavefunction.
- No spin dynamics or exchange-correlation treatment.
- No retarded fields (instantaneous Coulomb interaction used).
- No radiation reaction / synchrotron losses.

These effects can matter in extreme conditions, but this model is a robust baseline for physically meaningful electron dynamics.
