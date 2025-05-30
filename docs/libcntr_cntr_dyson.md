## Overview

The Dyson equation is a cornerstone of quantum many-body theory, providing a way to calculate the full Green's function \f$ G \f$ of an interacting system from a non-interacting Green's function \f$ G_0 \f$ (or a simpler reference) and the self-energy \f$ \Sigma \f$, which encapsulates all interaction effects. This module within the `libcntr` library provides a comprehensive suite of functions for solving Dyson equations on the Kadanoff-Baym contour.

The primary forms of the Dyson equation addressed are:
1.  **Integral form (typically for Matsubara components):** \f$ G = G_0 + G_0 \Sigma G \f$. This is often solved iteratively or by matrix inversion in frequency-momentum space.
2.  **Differential form (typically for real-time components):** \f$ [ i\frac{d}{dt} + \mu - H(t) ] G(t,t') - [\Sigma \ast G](t,t') = \delta(t,t') \f$. This is an integro-differential equation that can be solved by time-stepping.

The solvers are designed to work with `cntr::herm_matrix<T>` objects for Green's functions and self-energies, and `cntr::function<T>` for time-dependent Hamiltonians. Both serial and OpenMP-parallelized versions of key routines are provided.

## Key Components

The Dyson equation solvers are organized into functions for handling the Matsubara (imaginary time) part and the real-time part, with global functions orchestrating the complete solution.

**1. Matsubara Dyson Equation Solvers:**
   - **`dyson_mat(...)`**:
     *   **Purpose:** Solves the Dyson equation for the Matsubara Green's function \f$ G^M \f$.
     *   **Equation solved:** Effectively \f$ ( (G_0^M)^{-1} - \Sigma^M ) G^M = 1 \f$, where \f$ (G_0^M)^{-1} \f$ is related to \f$ (i\omega_n + \mu - H_0) \f$.
     *   **Methods (`method` parameter):**
         *   `CNTR_MAT_FOURIER` (default for some overloads): Uses discrete Fourier transforms. \f$ G(i\omega_n) = (i\omega_n + \mu - H_0 - \Sigma(i\omega_n))^{-1} \f$.
         *   `CNTR_MAT_FIXPOINT` (default for others): An iterative fixed-point solver, likely for \f$ G = ( (G_0)^{-1} - \Sigma_{corr} )^{-1} \f$ where \f$ G_0 \f$ might include a mean-field part of \f$ \Sigma \f$. Based on Volterra Integral Equation of the 2nd kind (VIE2) solvers.
         *   `CNTR_MAT_CG` (steepest descent / conjugate gradient): An iterative method for solving the matrix equation. Also based on VIE2 solvers.
     *   **Key Parameters:**
         *   `G`: (Output) `herm_matrix<T>` for the solution \f$ G^M \f$.
         *   `Sigma`: Input `herm_matrix<T>` for the self-energy \f$ \Sigma^M \f$.
         *   `mu`: Chemical potential.
         *   `H`: Input `cntr::function<T>` for the static Hamiltonian part \f$ H_0 \f$ (uses value at `t=-1`).
         *   `SigmaMF` (optional): `cntr::function<T>` for a mean-field part of the self-energy, which can be absorbed into \f$ G_0 \f$.
         *   `beta`: Inverse temperature.
         *   `SolveOrder` / `I` (`integration::Integrator<T>`): Integration order for internal convolutions if needed by iterative methods.
         *   `force_hermitian`: Boolean to enforce \f$ G^M(\tau) = (G^M(\beta-\tau))^\dagger \f$.
     *   No direct OpenMP version for `dyson_mat` itself is listed, but the underlying convolutions or VIE solvers might be parallelized.

**2. Real-Time Dyson Equation Solvers:**
   These functions solve the integro-differential equation: \f$ [ i\frac{d}{dt} + \mu - H(t) ] G(t,t') - \int_{\mathcal{C}} d\bar{t} \Sigma(t,\bar{t})G(\bar{t},t') = \delta(t,t') \f$.

   - **`dyson_start(...)`**:
     *   **Purpose:** Solves for \f$ G(t,t') \f$ for the initial `SolveOrder` real time steps. This is a specialized start-up procedure because the standard time-stepping formulas require a history of `SolveOrder` points. It assumes \f$ G^M \f$ is already computed and uses it for boundary conditions.
     *   **Key Parameters:** `G` (output), `mu`, `H` (time-dependent), `Sigma`, `beta`, `h` (real time step), `SolveOrder`.
     *   Internally calls `dyson_start_ret`, `dyson_start_tv`, `dyson_start_les`.

   - **`dyson_timestep(...)`**:
     *   **Purpose:** Advances the solution for \f$ G(t,t') \f$ by one real time-step `n`, assuming solutions for previous time steps up to `n-1` (and the initial `SolveOrder` steps) are known.
     *   **Key Parameters:** `n` (current time step), `G` (input/output), `mu`, `H`, `Sigma`, `beta`, `h`, `SolveOrder`.
     *   Internally calls `dyson_timestep_ret`, `dyson_timestep_tv`, `dyson_timestep_les`.
     *   **OpenMP Version:** `dyson_timestep_omp(...)` is available, which parallelizes the computation of retarded, lesser, and time-vertical components using functions like `dyson_timestep_ret_omp`, etc. The parallelization strategy often involves distributing calculations for different matrix elements or the second time argument across threads.

   - **`dyson(...)`**:
     *   **Purpose:** Global solver that orchestrates the full solution on the contour.
     *   **Steps:**
         1.  Calls `dyson_mat` to solve for the Matsubara Green's function.
         2.  If real-time steps are required (`G.nt() >= 0`), calls `dyson_start` to initialize the first `SolveOrder` real-time steps.
         3.  Iteratively calls `dyson_timestep` for `n` from `SolveOrder + 1` to `G.nt()`.
     *   **Key Parameters:** `G` (output), `mu`, `H`, `Sigma`, `beta`, `h`, `SolveOrder`, `matsubara_method`, `force_hermitian`.
     *   No direct OpenMP version of the global `dyson` is declared, but if an OpenMP-enabled `dyson_timestep` (or custom loop calling `dyson_timestep_omp`) is used internally by a user, the time-stepping part can be parallel.

**3. Specialized OpenMP Implementations:**
   - `cntr_dyson_omp_decl.hpp` and `cntr_dyson_omp_impl.hpp` define OpenMP versions for the core component calculations:
     *   `dyson_timestep_ret_omp`
     *   `dyson_timestep_tv_omp`
     *   `dyson_timestep_les_omp`
     *   `pseudodyson_timestep_tv_omp` and `pseudodyson_timestep_omp` for `herm_pseudo<T>` objects, suggesting a similar Dyson-like solver for a different matrix type.

## Important Variables/Constants

-   **`CNTR_DYSON_DECL_H`, `CNTR_DYSON_IMPL_H`, etc.**: Include guards.
-   **`CNTR_USE_OMP`**: Preprocessor macro, enables OpenMP versions.
-   **`CNTR_MAT_FOURIER`, `CNTR_MAT_FIXPOINT`, `CNTR_MAT_CG`**: Integer constants defining the solution method for the Matsubara Dyson equation.
-   **`MAX_SOLVE_ORDER`**: Default integration/differentiation order for simplified interfaces.
-   **Input Parameters:**
    *   `mu`: Chemical potential.
    *   `beta`: Inverse temperature.
    *   `h`: Real-time discretization step.
    *   `SolveOrder`: Numerical integration/differentiation order.
    *   `matsubara_method`: Method selector for `dyson_mat`.
    *   `force_hermitian`: Boolean flag for Matsubara solver.
    *   `omp_num_threads`: For OpenMP versions, specifies thread count.

## Usage Examples

The Dyson solvers are high-level tools for complex simulations.

**Conceptual Usage:**

```cpp
#include "cntr/cntr_dyson_decl.hpp" // Or main cntr.hpp include
#include "cntr/cntr_herm_matrix_decl.hpp"
#include "cntr/cntr_function_decl.hpp"
#ifdef CNTR_USE_OMP
#include "cntr/cntr_dyson_omp_decl.hpp"
#endif

using GMatrix = cntr::herm_matrix<T>;
using CFunction = cntr::function<T>;

// Initialize Green's functions, self-energy, Hamiltonian, and parameters
GMatrix G_solution, Sigma_input;
CFunction H_hamiltonian;
T mu_chemical_potential = 0.0;
T beta_inverse_temp = 10.0;
T h_real_timestep = 0.02;
int num_real_timesteps = 100;
int num_matsubara_timesteps = 128;
int orbital_size = 1; // Example for a 1x1 matrix
int integration_order = 3;

// Setup G_solution, Sigma_input, H_hamiltonian with dimensions and initial data
// G_solution: (num_real_timesteps, num_matsubara_timesteps, orbital_size, fermi_flag)
// Sigma_input: same as G_solution
// H_hamiltonian: (num_real_timesteps, orbital_size)

// Example 1: Solve the full Dyson equation (Matsubara + real-time)
cntr::dyson(G_solution, mu_chemical_potential, H_hamiltonian, Sigma_input,
            beta_inverse_temp, h_real_timestep, integration_order,
            CNTR_MAT_FIXPOINT, true);

// Example 2: Solve only the Matsubara part using Fourier method
GMatrix G_matsubara_only( -1, num_matsubara_timesteps, orbital_size, G_solution.sig() );
cntr::dyson_mat(G_matsubara_only, mu_chemical_potential, H_hamiltonian, Sigma_input,
                beta_inverse_temp, integration_order, CNTR_MAT_FOURIER, true);

#ifdef CNTR_USE_OMP
// Example 3: Time-stepping part using OpenMP (conceptual, assumes G_matsubara is already computed and set in G_solution_omp)
GMatrix G_solution_omp = G_solution; // Copy initial state if needed
cntr::dyson_start(G_solution_omp, mu_chemical_potential, H_hamiltonian, Sigma_input, beta_inverse_temp, h_real_timestep, integration_order);
int num_threads = 4;
for (int n = integration_order + 1; n <= num_real_timesteps; ++n) {
    cntr::dyson_timestep_omp(num_threads, n, G_solution_omp, mu_chemical_potential, H_hamiltonian, Sigma_input,
                             beta_inverse_temp, h_real_timestep, integration_order);
}
#endif
```

## Dependencies and Interactions

-   **`cntr_herm_matrix_decl.hpp`, `cntr_function_decl.hpp`**: For the primary data structures (`herm_matrix<T>`, `function<T>`).
-   **`cntr_convolution_impl.hpp` (and `_decl`)**: Essential for the \f$ \Sigma \ast G \f$ terms in the Dyson equation. The real-time Dyson solvers heavily rely on these convolution routines for different components.
-   **`integration.hpp`**: Provides the `integration::Integrator<T>` class and numerical rules (weights) for both integration (in convolutions) and differentiation (implicitly in the \f$ i d/dt \f$ term).
-   **`cntr_matsubara_impl.hpp`**: Used by `dyson_mat_fourier` for Discrete Fourier Transforms.
-   **`cntr_equilibrium_decl.hpp`**: `green_from_H` is used to obtain the non-interacting \f$ G_0 \f$ in some Matsubara solution methods (e.g., fixpoint, steep).
-   **`cntr_vie2_decl.hpp`**: Provides solvers for Volterra Integral Equations of the 2nd kind, used in `dyson_mat_fixpoint` and `dyson_mat_steep`.
-   **`cntr_elements.hpp`, `cntr_utilities_decl.hpp`**: For low-level matrix/element operations and utilities like forcing hermiticity.
-   **OpenMP**: If `CNTR_USE_OMP == 1`, the OpenMP runtime is a dependency for parallel versions.
-   **`cntr_herm_pseudo_decl.hpp`, `cntr_pseudo_convolution_decl.hpp`, `cntr_pseudodyson_decl.hpp`**: Indicated by includes in `cntr_dyson_omp_impl.hpp`, suggesting similar Dyson solvers might exist for `herm_pseudo<T>` objects, handled by `pseudodyson_timestep_omp`.

The Dyson equation solvers are among the most complex and high-level components in `libcntr`, tying together functionality from convolution, numerical integration/differentiation, and specialized matrix equation solvers. They are central to simulating interacting quantum systems.
