# Numerical Analysis on Equilibrium Composition of Air Plasma

**Author:** Wataru Furutaka

A Fortran-based program to calculate the chemical equilibrium composition of high-temperature air plasma using Gibbs free energy minimization and a damped Newton-Raphson method with backtracking line search.

## Introduction

This project simulates the chemical equilibrium composition of high-temperature air plasma (from 298.15 K up to 20,000 K) under varying pressures. It models the thermodynamic properties of 11 chemical species and solves the non-linear system of chemical equilibrium and mass/charge conservation equations.

The solver utilizes a log-transformed damped Newton-Raphson method with a backtracking line search (Armijo condition) to guarantee convergence even at lower temperatures where equilibrium constants span many orders of magnitude. Linear algebra system solving at each iteration is offloaded to LAPACK via Intel MKL.

## Thermodynamic Modeling

### NASA 9-Coefficients Polynomials

Thermodynamic properties (specific heat $C_p^\circ$, enthalpy $H^\circ$, and entropy $S^\circ$) for each species are evaluated using the NASA 9-coefficients polynomial fit. The coefficients are partitioned into three temperature intervals (200K–1000K, 1000K–6000K, and 6000K–20000K) sourced from the NASA CEA database (derived from `thermo.inp`).

The non-dimensional equations are:

- **Specific Heat**

  $$
  \frac{C_p^\circ}{R} = a_1 T^{-2} + a_2 T^{-1} + a_3 + a_4 T + a_5 T^2 + a_6 T^3 + a_7 T^4
  $$

- **Enthalpy**

  $$
  \frac{H^\circ}{RT} = -a_1 T^{-2} + a_2 T^{-1} \ln T + a_3 + \frac{a_4}{2} T + \frac{a_5}{3} T^2 + \frac{a_6}{4} T^3 + \frac{a_7}{5} T^4 + \frac{b_1}{T}
  $$

- **Entropy**

  $$
  \frac{S^\circ}{R} = -\frac{a_1}{2} T^{-2} - a_2 T^{-1} + a_3 \ln T + a_4 T + \frac{a_5}{2} T^2 + \frac{a_6}{3} T^3 + \frac{a_7}{4} T^4 + b_2
  $$

- **Gibbs Free Energy**

  $$
  G^\circ = RT \left( \frac{H^\circ}{RT} - \frac{S^\circ}{R} \right)
  $$

Where $a_1 \dots a_7$ are the polynomial coefficients, and $b_1, b_2$ are integration constants. These parameters are stored in `src/mod_constants.f90`.

## Chemical Equilibrium System

The plasma is assumed to consist of 11 species:

| Index | Species Name | Type | Description |
|-------|--------------|------|-------------|
| 1 | O | Atom | Oxygen Atom |
| 2 | N | Atom | Nitrogen Atom |
| 3 | NO | Compound | Nitric Oxide |
| 4 | O2 | Molecule | Oxygen Molecule |
| 5 | N2 | Molecule | Nitrogen Molecule |
| 6 | O+ | Ion | Oxygen Ion |
| 7 | N+ | Ion | Nitrogen Ion |
| 8 | O2+ | Ionized Mol | Oxygen Molecular Ion |
| 9 | N2+ | Ionized Mol | Nitrogen Molecular Ion |
| 10 | NO+ | Ion | Nitric Oxide Ion |
| 11 | e- | Atom | Electron |

### Governing Equations

To find the equilibrium composition (partial pressures $P_i$), we solve a system of 11 equations consisting of 8 chemical equilibrium constraints, 1 total pressure constraint, 1 charge neutrality constraint, and 1 mass conservation constraint for the Nitrogen-to-Oxygen ratio.

1. **Oxygen Cleavage:** $O_2 \rightleftharpoons 2O$

   $$
   K_{p, O_2} P_{O_2} = P_O^2
   $$

2. **Nitrogen Cleavage:** $N_2 \rightleftharpoons 2N$

   $$
   K_{p, N_2} P_{N_2} = P_N^2
   $$

3. **NO Dissociation:** $NO \rightleftharpoons N + O$

   $$
   K_{p, NO} P_{NO} = P_N P_O
   $$

4. **Oxygen Ionization:** $O^+ + e^- \rightleftharpoons O$

   $$
   K_{p, O^+} P_{O^+} P_{e^-} = P_O
   $$

5. **Nitrogen Ionization:** $N^+ + e^- \rightleftharpoons N$

   $$
   K_{p, N^+} P_{N^+} P_{e^-} = P_N
   $$

6. **Oxygen Molecule Ionization:** $O_2^+ + e^- \rightleftharpoons 2O$

   $$
   K_{p, O_2^+} P_{O_2^+} P_{e^-} = P_O^2
   $$

7. **Nitrogen Molecule Ionization:** $N_2^+ + e^- \rightleftharpoons 2N$

   $$
   K_{p, N_2^+} P_{N_2^+} P_{e^-} = P_N^2
   $$

8. **NO Ionization:** $NO^+ + e^- \rightleftharpoons N + O$

   $$
   K_{p, NO^+} P_{NO^+} P_{e^-} = P_N P_O
   $$

9. **Total Pressure Constraint**

   $$
   P_{\text{tot}} = \sum_{i=1}^{11} P_i
   $$

10. **Charge Neutrality**

   $$
   P_{e^-} = P_{O^+} + P_{N^+} + P_{O_2^+} + P_{N_2^+} + P_{NO^+}
   $$

11. **Nitrogen-to-Oxygen Ratio** (mass conservation, approx. 78:21)

   $$
   \frac{N_N}{N_O} = \frac{78}{21} \implies 78 \cdot N_O - 21 \cdot N_N = 0
   $$

   Where $N_O$ and $N_N$ are the total abundance of oxygen and nitrogen atoms across all species.

## Numerical Solver

To handle the large variation in partial pressures (which can drop to $10^{-30}$ atm or lower), the solver computes the log-transformed pressures:

$$
x_i = \ln\left(P_i / P_{\text{atm}}\right)
$$

This transformation ensures that $P_i$ remains strictly positive and improves the conditioning of the Jacobian matrix.

### Continuation Method & OpenMP Parallelization

To drastically accelerate convergence, the solver employs a continuation method: after the first temperature step, the converged solution from the previous temperature is used as the initial guess for the next. This reduces the required Newton-Raphson iterations to just 1-3 per temperature step.

Furthermore, the outer loop over varying atmospheric pressures is completely parallelized using OpenMP, distributing the independent computation paths across all available CPU cores.

### Damped Newton-Raphson Method

At each iteration $k$, the linear system is solved for the correction vector $\Delta x$:

$$
J(x^{(k)}) \Delta x = -f(x^{(k)})
$$

Using the LAPACK double-precision general solver (`dgesv`).

### Backtracking Line Search (Armijo Rule)

To prevent divergence when far from the solution, a backtracking line search is performed:

$$
x^{(k+1)} = x^{(k)} + \alpha \Delta x
$$

where the step size $\alpha$ starts at $1.0$ and is successively halved ($\alpha \leftarrow 0.5 \alpha$) until the Residual Sum of Squares (RSS) decreases sufficiently:

$$
\text{RSS}(x^{(k)} + \alpha \Delta x) \le \text{RSS}(x^{(k)}) - c_1 \alpha \|f(x^{(k)})\|^2
$$

with $c_1 = 10^{-4}$.

## Directory Structure

```text
Gibbs/
├── Makefile                # Multi-configuration compiler settings
├── src/                    # Fortran source files
│   ├── main.f90            # Main runner calculating constants and compositions
│   ├── mod_constants.f90   # Physical constants and NASA polynomial coefficients
│   ├── mod_types.f90       # OOP definition of chemical species and reaction types
│   ├── mod_types_functions.f90    # Submodule for thermodynamic calculations (H, S, G)
│   ├── mod_types_subroutines.f90 # Submodule for species initialization
│   └── mod_newton_raphson.f90    # Backtracking, function definitions, NR solver
├── test/                   # Unit tests & verification programs
│   ├── gibbses.f90         # Tests individual species Gibbs energy outputs
│   └── newton-test/        # Experimental & mock tests for solver validation
├── experiments/            # Visualizations and scratchpads
│   └── gnuplot/
│       └── test.gp         # Gnuplot script for plotting mole fractions
└── output/                 # Destination for generated data tables (*.dat)
```

## Building and Running

### Prerequisites

- Intel Fortran Compiler (`ifx`)
- Intel oneAPI Math Kernel Library (MKL) for LAPACK support

### Compiling

To build the executable in release mode (default):

```bash
make
```

To build with debug checks and symbols:

```bash
make BUILD=debug
```

### Running the Simulation

Run the compiled binary:

```bash
./bin/calc_pressure
```

This will calculate:

1. The equilibrium constants for all species at different temperatures (outputted to `output/k_test.dat`).
2. The partial pressures (mole fractions) for atmospheric pressures of 1.0, 0.1, and 0.01 atm across temperatures 298.15 K to 20,000 K. The results are automatically saved to dynamically named, separated output files incorporating the pressure and gas ratio (e.g. `output/1p00atm78v21.dat`, `output/0p10atm78v21.dat`).

### Plotting Results

You can plot the mole fraction distribution as a function of temperature using the provided Gnuplot script:

```bash
cd output
gnuplot ../experiments/gnuplot/test.gp
```

This plots the species distribution over the temperature range [0:20000] K.

## Planned Features & Future Enhancements

### Object-Oriented Refactoring & Encapsulation

Refactor the Newton-Raphson solver and backtracking line search routines into clean Fortran modules. Target complete encapsulation of solver parameters, Jacobians, and residual state variables in object-oriented structures (derived types) to reduce global/module variables.

### Solver Customization

- Implement configurable convergence tolerance and maximum iteration limits.
- Allow swapping or selecting between different iterative schemes (e.g. standard Newton-Raphson vs. Damped Newton-Raphson with backtracking line search).

### Dynamic Gas Composition

- Allow custom gas mixture ratios (e.g. arbitrary $N_2$ to $O_2$ concentrations) to be passed via preprocessor definitions (e.g. `-Dratio_N=...`) at compile time.

### Output Formatting

- Add descriptive column headers (chemical species names) to the first line of the output data files and adapt the Gnuplot script to read columns dynamically.

### Logging and Diagnostics

- Route status/convergence warnings and LAPACK `dgesv` solver errors directly to `stderr` or a dedicated run log file for better monitoring.
