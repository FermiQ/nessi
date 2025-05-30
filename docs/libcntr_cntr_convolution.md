## Overview

The files `cntr_convolution_decl.hpp` and `cntr_convolution_impl.hpp` are core components of the `libcntr` library, providing extensive functionality for computing convolutions of contour Green's functions (or similar matrix-like objects). The general form of the convolution is \f$ C(t,t') = \int_{\mathcal{C}} d\bar{t} A(t,\bar{t}) B(\bar{t},t') \f$, often referred to as \f$ C = A \ast B \f$. The library also supports convolutions involving an additional function \f$ F \f$, as \f$ C(t,t') = \int_{\mathcal{C}} d\bar{t} A(t,\bar{t}) F(\bar{t}) B(\bar{t},t') \f$ (\f$ C = A \ast F \ast B \f$).

These operations are fundamental in solving many-body problems, particularly within the Kadanoff-Baym formalism, for calculating self-energies, polarizations, and other physical quantities. The implementation handles the complexities of the contour, including Matsubara (imaginary time) and real-time (retarded, lesser, advanced, greater, time-vertical/left-mixing) components, using numerical integration techniques.

## Key Components

The convolution functionality is exposed through a large set of templated functions, primarily operating on `cntr::herm_matrix<T>` objects and optionally `cntr::function<T>` objects. Many functions have overloads for handling general matrices, Hermitian matrices, and versions with/without an intermediate function \f$ F \f$. OpenMP parallelized versions are also available for performance enhancement.

**Main Function Groups:**

1.  **`convolution_timestep` and `convolution`:**
    *   **Purpose:** Compute the full contour convolution \f$ C = A \ast B \f$ or \f$ C = A \ast F \ast B \f$.
    *   `convolution_timestep(int n, ...)`: Computes the convolution for a specific real time-step `n` (or Matsubara part if `n == -1`).
    *   `convolution(...)`: Computes the full convolution across all time steps (Matsubara and real-time by iterating `convolution_timestep`).
    *   **Key Parameters:**
        *   `C`: Output `herm_matrix<T>` for the result.
        *   `A`, `B`: Input `herm_matrix<T>` objects.
        *   `Acc`, `Bcc`: Optional `herm_matrix<T>` objects representing the adjoints of `A` and `B`. If not provided, `A` and `B` are assumed Hermitian.
        *   `ft` (optional): A `cntr::function<T>` object for the \f$ F \f$ term in \f$ A \ast F \ast B \f$. Raw pointers `f0` (Matsubara part of F) and `ft` (real-time part of F) are used in some internal/private versions.
        *   `beta`: Inverse temperature.
        *   `h`: Real-time step size.
        *   `SolveOrder` (or `integration::Integrator<T> I`): Specifies the order/method for numerical integration.
    *   **Forms:**
        *   `C = A \ast B`: `convolution_timestep(n, C, A, Acc, B, Bcc, ...)`
        *   `C = A \ast B` (A, B Hermitian): `convolution_timestep(n, C, A, B, ...)`
        *   `C = A \ast F \ast B`: `convolution_timestep(n, C, A, Acc, ft, B, Bcc, ...)`
        *   `C = A \ast F \ast B` (A, B Hermitian): `convolution_timestep(n, C, A, ft, B, ...)`

2.  **`convolution_density_matrix`:**
    *   **Purpose:** Computes specific components relevant to the density matrix, typically \f$ \rho = -i (A \ast B)^< \f$ or \f$ \rho = -i (A \ast F \ast B)^< \f$ at a specific time-step `tstp`.
    *   **Key Parameters:** Similar to `convolution_timestep`, but output is `cdmatrix& rho` (likely a complex Eigen matrix for the density matrix) or a raw pointer `std::complex<T>* rho`.
    *   **Forms:** Similar variations for with/without function `ft` and with/without explicit adjoints `Acc, Bcc`.

3.  **OpenMP Parallelized Versions:**
    *   Functions like `convolution_omp`, `convolution_timestep_omp`, `convolution_matsubara_omp` mirror their serial counterparts but include an `omp_num_threads` parameter.
    *   These aim to speed up calculations by parallelizing loops over time steps or matrix elements.

4.  **Internal Helper and Specialized Functions (mostly `@private`):**
    *   `convolution_matsubara(...)`: Computes the Matsubara component of the convolution.
    *   `convolution_timestep_ret(...)`: Computes the retarded component.
    *   `convolution_timestep_tv(...)`: Computes the time-vertical (left-mixing) component.
    *   `convolution_timestep_les(...)`: Computes the lesser component by combining terms like \f$A^R B^< \f$, \f$A^< B^A \f$, and \f$A^\rceil B^\lceil \f$.
        *   Further specialized helpers like `convolution_timestep_les_tvvt`, `convolution_timestep_les_lesadv`, `convolution_timestep_les_retles` calculate these individual terms.
    *   `matsubara_integral_1`, `matsubara_integral_2`: Low-level functions for performing the Matsubara frequency/time integrals using an `integration::Integrator`.
    *   `incr_convolution_*`: Incremental versions that add the result of a convolution (multiplied by a scalar `alpha`) to an existing `C` matrix. These are used by the OpenMP versions to accumulate partial results.
    *   Many functions are templated with `SIZE1` (e.g., `convolution_matsubara_dispatch<T, GG, SIZE1>`) to specialize for scalar (SIZE1=1) or general matrix (SIZE1=LARGESIZE) operations.

**Integration Scheme:**
The convolutions rely on an `integration::Integrator<T>` object (or a `SolveOrder` parameter which likely constructs one). This object provides weights for numerical integration (e.g., Gregory weights) and handles boundary conditions, crucial for accuracy.

## Important Variables/Constants

-   **`CNTR_CONVOLUTION_DECL_H` / `CNTR_CONVOLUTION_IMPL_H`**: Include guards.
-   **`CNTR_USE_OMP`**: Preprocessor macro to enable/disable OpenMP parallelized versions.
-   **`LARGESIZE`**: A placeholder (likely a large integer or -1) used in template specializations to denote general matrix sizes as opposed to fixed small sizes (e.g., SIZE1=1).
-   **`MAX_SOLVE_ORDER`**: Default integration order for simplified interfaces.
-   **Input Parameters:**
    *   `beta`: Inverse temperature.
    *   `h`: Real-time discretization step.
    *   `n` or `tstp`: The specific real time-step index for `_timestep` functions. `n=-1` usually triggers Matsubara calculation.
    *   `SolveOrder`: Integer defining the order of the numerical integration method.
    *   `omp_num_threads`: Number of threads for OpenMP versions.

## Usage Examples

The convolution functions are intended for use in complex physics simulations where Green's functions are manipulated.

**Conceptual Usage:**

```cpp
#include "cntr/cntr_convolution_decl.hpp" // Or the main cntr.hpp
#include "cntr/cntr_herm_matrix_decl.hpp"
#include "cntr/cntr_function_decl.hpp"

// Assuming T is float or double
using GMatrix = cntr::herm_matrix<T>;
using CFunction = cntr::function<T>;

// Initialize Green's functions and parameters
GMatrix G_A, G_B, G_C;
GMatrix G_A_adj, G_B_adj; // Adjoints, if needed
CFunction F_func;         // Optional intermediate function

T beta = 1.0;    // Inverse temperature
T h_real = 0.01; // Real time step
int nt = 100;    // Number of real time steps
int ntau = 100;  // Number of Matsubara steps
int solve_order = 3; // Integration order

// ... (Initialize G_A, G_B, F_func with appropriate dimensions and data)
// ... (G_C should be initialized with same parameters as A, B)

// Example 1: Full convolution C = A * B (A, B Hermitian)
// Uses default integration order MAX_SOLVE_ORDER if not specified
cntr::convolution(G_C, G_A, G_B, beta, h_real, solve_order);

// Example 2: Convolution C = A * F * B at a specific time step `current_step`
int current_step = 10;
// Assuming A, B are not necessarily Hermitian, provide adjoints
// cntr::convolution_timestep(current_step, G_C, G_A, G_A_adj, F_func, G_B, G_B_adj, beta, h_real, solve_order);
// Simplified interface (if A, B are hermitian)
cntr::convolution_timestep(current_step, G_C, G_A, F_func, G_B, beta, h_real, solve_order);


// Example 3: Matsubara part of C = A * B
cntr::convolution_timestep(-1, G_C, G_A, G_B, beta, h_real, solve_order); // h_real is not used for n=-1

// Example 4: Calculate density matrix rho = -i * (A*B)^< at tstp
cntr::cdmatrix rho(G_A.size1(), G_A.size1()); // Eigen complex matrix
int t_idx_for_rho = 5;
cntr::convolution_density_matrix(t_idx_for_rho, rho, G_A, G_B, beta, h_real, solve_order);

#ifdef CNTR_USE_OMP
// Example 5: Parallel convolution C = A * B
int num_threads = 4;
cntr::convolution_omp(num_threads, G_C, G_A, G_B, beta, h_real, solve_order);
#endif
```

The user chooses the appropriate function based on:
- Whether the full contour or a specific time-step is needed.
- Whether an intermediate function \f$ F \f$ is part of the convolution.
- Whether the input matrices `A` and `B` are Hermitian (simplified calls) or require explicit adjoints (`Acc`, `Bcc`).
- Whether OpenMP parallelism is desired.

## Dependencies and Interactions

-   **`cntr_global_settings.hpp`**: Provides global settings, possibly including `LARGESIZE`, `MAX_SOLVE_ORDER`, `CNTR_USE_OMP`.
-   **`integration.hpp`**: Crucial for providing the `integration::Integrator<T>` class and related integration weights/methods (e.g., `gregory_omega`, `gregory_weights`, `poly_integration`).
-   **`cntr_elements.hpp`**: Likely contains low-level element-wise operations for matrices (e.g., `element_incr`, `element_mult`, `element_conj`).
-   **`cntr_herm_matrix_decl.hpp`**: Defines `herm_matrix<T>`, the primary data structure these functions operate on.
-   **`cntr_function_decl.hpp`**: Defines `function<T>`, used for the intermediate function \f$ F \f$ in \f$ A \ast F \ast B \f$ convolutions.
-   **`eigen_map.hpp`**: Used for `cdmatrix` type (likely `Eigen::MatrixXcd`) in `convolution_density_matrix`.
-   **C++ Standard Library**: `<vector>`, `<complex>`, `<cassert>`.
-   **OpenMP**: If `CNTR_USE_OMP == 1`, it depends on the OpenMP runtime library.

The convolution routines are high-level operations that build upon more fundamental matrix operations, element manipulations, and numerical integration schemes provided by other parts of the `libcntr` library and its dependencies. They are essential for constructing and solving integral equations like the Dyson or Bethe-Salpeter equations in quantum many-body theory.
