## Overview

The `cntr_equilibrium_decl.hpp` and `cntr_equilibrium_impl.hpp` files in the `libcntr` library provide essential tools for calculating and constructing Green's functions and related quantities in thermal equilibrium. These functionalities are critical for initializing non-equilibrium simulations (where a system often starts from an equilibrium state) or for performing purely equilibrium studies.

The module offers several ways to define and compute equilibrium Green's functions:
1.  From a given Density of States (DOS), using Hilbert transforms (via Fourier transforms) combined with Fermi-Dirac or Bose-Einstein statistics.
2.  From a given quadratic Hamiltonian (either time-independent or time-dependent), by solving the corresponding non-interacting Schrödinger equation and populating states according to equilibrium distributions.
3.  Specific functions for simple bosonic systems.

Standard Fermi-Dirac and Bose-Einstein distribution functions are also provided.

## Key Components

**1. Distribution Functions:**
   - **`fermi(beta, omega)`**:
     *   Calculates the Fermi-Dirac distribution \f$ n_F(\omega) = 1 / (1 + e^{\beta \omega}) \f$.
     *   Available for single `omega` or `dvector` of energies.
   - **`bose(beta, omega)`**:
     *   Calculates the Bose-Einstein distribution \f$ n_B(\omega) = 1 / (e^{\beta \omega} - 1) \f$.
     *   Available for single `omega` or `dvector` of energies.
   - **`fermi_exp(beta, tau, omega)` / `bose_exp(beta, tau, omega)` (@private):**
     *   Calculate terms like \f$ e^{\omega \tau} n_F(\omega) \f$ or \f$ e^{\omega \tau} n_B(\omega) \f$, which are building blocks for Matsubara Green's functions.

**2. Green's Functions from Density of States (DOS):**
   - **`green_equilibrium(G, dos, beta, h, mu, limit, nn)`**:
     *   **Purpose:** Computes the full equilibrium contour Green's function (`herm_matrix<T>& G`) given a spectral density of states (`dos_function& dos`).
     *   **Method:** Uses Fourier transforms (specifically `fourier::adft_func`) of the DOS combined with the appropriate statistical distribution (Fermi/Bose, determined by `G.sig()`) and phase factors to construct each component (Matsubara, retarded, lesser, time-vertical).
     *   The `dos_function` object must provide `operator()(double omega)` and energy range `lo_`, `hi_`.
     *   **Parameters:**
         *   `G`: Output `herm_matrix<T>`.
         *   `dos`: User-supplied DOS object.
         *   `beta`: Inverse temperature.
         *   `h`: Real-time step.
         *   `mu`: Chemical potential (shifts energies in DOS: \f$ \omega \rightarrow \omega - \mu \f$).
         *   `limit`, `nn`: Parameters for the adaptive Fourier transform.
   - **`green_equilibrium_mat(G, dos, beta, limit, nn, mu)`**:
     *   Calculates only the Matsubara component of the Green's function from the DOS.
   - **Predefined DOS classes:**
     *   `bethedos`: Semicircular density of states. Convenience wrappers `green_equilibrium_bethe` and `green_equilibrium_mat_bethe` use this.
     *   `ohmic`: Ohmic bath DOS: \f$ \omega^2 e^{-\omega/\omega_c} \f$.
     *   `smooth_box`: A box-like DOS smoothed with Fermi functions.

**3. Green's Functions from Hamiltonian (`green_from_H`):**
   - **Purpose:** Calculates the non-interacting Green's function \f$ G_0 \f$ for a system described by a quadratic Hamiltonian \f$ H(t) - \mu N \f$.
   - **Overloads for constant Hamiltonian (`cdmatrix& eps` for \f$H_0\f$):**
     *   `green_from_H(G, mu, eps, beta, h)`: For full `herm_matrix<T>`.
     *   `green_from_H(tstp, G, mu, eps, beta, h)`: For `herm_matrix_timestep<T>` or to set a specific timestep of a `herm_matrix<T>`.
     *   **Method:** Diagonalizes \f$ H_0 - \mu I \f$ to get eigenvalues \f$ \epsilon_k \f$ and eigenvectors \f$ \phi_k \f$.
         *   Matsubara: \f$ G^M(\tau) = -\sum_k \phi_k \phi_k^\dagger e^{-\epsilon_k \tau} / (1 + e^{-\beta \epsilon_k}) \f$ (fermions).
         *   Real-time: Uses propagators \f$ U(t) = e^{-i(H_0-\mu)t} \f$ and equilibrium occupations.
   - **Overloads for time-dependent Hamiltonian (`cntr::function<T>& eps` for \f$H(t)\f$):**
     *   `green_from_H(G, mu, eps, beta, h, SolveOrder, cf_order)`: For full `herm_matrix<T>`.
     *   `green_from_H(tstp, G, mu, eps, beta, h, fixHam, SolveOrder, cf_order)`: For `herm_matrix_timestep<T>` or a specific `tstp` of `herm_matrix<T>`.
     *   **Method:** Uses high-order commutator-free (CF) Magnus-type integrators (e.g., CF2:1, CF4:2 from Alvermann & Fehske, JCP 230, 5930 (2011)) to compute the time evolution operator \f$ U(t) \f$ by solving \f$ i \partial_t U(t,\bar{t}) = (H(t)-\mu) U(t,\bar{t}) \f$.
         *   The equilibrium density matrix is determined by \f$ H(t=0)-\mu \f$.
         *   Green's functions are then constructed as \f$ G_0(t,t') = -i U(t,0) \langle T_C c(0) c^\dagger(0) \rangle_{H(0)} U(0,t') \f$.
     *   **Parameters:**
         *   `SolveOrder`: Order for interpolation/extrapolation of \f$ H(t) \f$ if needed by the CF integrator.
         *   `cf_order`: Order of the commutator-free expansion (2 or 4).
         *   `fixHam`: If true, assumes \f$ H(t) \f$ is known for all required times (no extrapolation).

**4. Bosonic Green's Functions for Single Pole:**
   - **`green_single_pole_XX(D0, w, beta, h)`**:
     *   Calculates the equilibrium Green's function \f$ D_0(t,t') = -i \langle T_C X(t) X(t') \rangle \f$ for the displacement operator \f$ X = (b+b^\dagger)/\sqrt{2} \f$ of a single boson mode with frequency `w`.
   - **`green_single_pole_XX_timestep(tstp, D0, w, beta, h)`**:
     *   Calculates a specific timestep `tstp` for the above bosonic Green's function.
     *   Internal helpers `green_single_pole_XX_mat` (for Matsubara) and `green_single_pole_XX_timestep` (for real-time components) perform the actual calculations using Bose distribution functions and complex exponentials.

## Important Variables/Constants

-   **`CNTR_EQUILIBRIUM_DECL_H` / `CNTR_EQUILIBRIUM_IMPL_H`**: Include guards.
-   **Input Parameters common to many functions:**
    *   `beta`: Inverse temperature.
    *   `mu`: Chemical potential.
    *   `h`: Real-time step size (for functions constructing real-time components).
-   **Numerical Parameters:**
    *   `limit`, `nn`: Control accuracy and sampling for Fourier transforms in DOS-based calculations.
    *   `SolveOrder`: Order for polynomial interpolation/extrapolation of time-dependent Hamiltonians.
    *   `cf_order`: Order (2 or 4) of the commutator-free Magnus expansion for `green_from_H` with time-dependent \f$ H(t) \f$.

## Usage Examples

This module is typically used to set up the initial state of a system before applying a non-equilibrium perturbation or to study equilibrium properties.

**Conceptual Usage:**

```cpp
#include "cntr/cntr_equilibrium_decl.hpp"
#include "cntr/cntr_herm_matrix_decl.hpp"
#include "cntr/cntr_function_decl.hpp" // For time-dependent H

using GMatrix = cntr::herm_matrix<T>;
using CFunction = cntr::function<T>;

// Parameters
T beta_val = 10.0;
T h_step = 0.01;
T mu_val = 0.0;
int orbital_dim = 1;
int n_real_steps = 50;
int n_matsubara_steps = 64;

// Example 1: Equilibrium G from Bethe DOS
GMatrix G_bethe(n_real_steps, n_matsubara_steps, orbital_dim, cntr::FERMI);
cntr::green_equilibrium_bethe(G_bethe, beta_val, h_step, 100, 20, mu_val);

// Example 2: Equilibrium G0 for a constant Hamiltonian H0
GMatrix G0_const_H(n_real_steps, n_matsubara_steps, orbital_dim, cntr::FERMI);
cntr::cdmatrix H_matrix(orbital_dim, orbital_dim);
// ... (initialize H_matrix, e.g., H_matrix(0,0) = some_energy_offset)
cntr::green_from_H(G0_const_H, mu_val, H_matrix, beta_val, h_step);

// Example 3: Equilibrium G0 for a time-dependent Hamiltonian H(t)
GMatrix G0_td_H(n_real_steps, n_matsubara_steps, orbital_dim, cntr::FERMI);
CFunction H_t_func(n_real_steps, orbital_dim);
// ... (initialize H_t_func, e.g., H_t_func.set_value(t, H_matrix_at_t) )
// H_t_func.set_value(-1, H_matrix_at_tau) for Matsubara Hamiltonian (usually H(0))
int solve_order_interp = 3;
int cf_expansion_order = 2; // Commutator-Free expansion order
cntr::green_from_H(G0_td_H, mu_val, H_t_func, beta_val, h_step, solve_order_interp, cf_expansion_order);

// Example 4: Bosonic equilibrium Green's function for X operator
GMatrix D0_boson(n_real_steps, n_matsubara_steps, 1, cntr::BOSE);
T boson_freq = 1.0;
cntr::green_single_pole_XX(D0_boson, boson_freq, beta_val, h_step);
```

## Dependencies and Interactions

-   **`cntr_global_settings.hpp`**: Global type definitions and constants.
-   **`fourier.hpp`**: Provides `fourier::adft_func` used in DOS-based Green's function calculations.
-   **`integration.hpp`**: Provides `integration::Integrator<T>` for interpolation/extrapolation weights in `green_from_H` with time-dependent Hamiltonians.
-   **`cntr_elements.hpp`**: For low-level element-wise operations.
-   **`cntr_herm_matrix_decl.hpp`, `cntr_herm_matrix_timestep_decl.hpp`, `cntr_function_decl.hpp`**: Defines the primary data structures.
-   **`cntr_utilities_decl.hpp`**: May be used for utilities like `force_matsubara_hermitian` if applicable, though not directly visible in equilibrium setup.
-   **Eigen library**: Used for `cdmatrix` (likely `Eigen::MatrixXcd`) and its linear algebra operations (e.g., `SelfAdjointEigenSolver`, `.exp()`).

The equilibrium module provides foundational Green's functions that can then be used as \f$ G_0 \f$ in Dyson equation solvers or as initial conditions for time-propagation of Kadanoff-Baym equations.
