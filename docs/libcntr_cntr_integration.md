## Overview

The `integration.hpp` file (located at `libcntr/cntr/integration.hpp`) is a core utility module within the `libcntr` library. Its primary purpose is to provide a consistent and efficient way to access pre-calculated weights for various numerical integration and differentiation schemes. These schemes are based on polynomial approximations of functions sampled at equally spaced points and are fundamental for many algorithms in `libcntr`, such as solving time evolution equations (e.g., Kadanoff-Baym equations) and computing convolution integrals.

The module defines the `integration::Integrator<T>` class, which encapsulates these weights for a given order of accuracy.

## Key Components

**Namespace: `integration`**

**1. `Integrator<T>` Class:**
   - **Template Parameter `T`**: Typically `double` or `float`, specifying the numerical precision of the weights.
   - **Purpose**: Acts as a repository and provider for various numerical quadrature and differentiation coefficients. An `Integrator` object is initialized for a specific order `k`.
   - **Constructor**: `Integrator(int k = 0)`
     *   Initializes the integrator for a given order `k`. The maximum order for Gregory integration is `GREGORY_KMAX = 8`, though other polynomial methods might be practically limited to lower orders (e.g., up to 5, as suggested by comments).
     *   Loads pre-calculated weights by calling internal `read_*` functions.
   - **Key Methods for Accessing Weights:**
     *   `get_k()` or `k()`: Returns the order `k` for which the integrator instance was initialized.
     *   `gregory_weights(int n, int j)`: Returns Gregory quadrature weights \f$w_{n,j}\f$ for integrating \f$\int_0^{x_n} f(x) dx = h \sum_j w_{n,j} f(x_j)\f$. Handles boundary conditions for small `n` and uses `gregory_omega` for larger `n`.
     *   `gregory_omega(int j)`: Returns the specific Gregory weights \f$\omega_j\f$ used at the start and end of the integration interval when the interval is sufficiently large.
     *   `bd_weights(int l)`: Returns coefficients \f$c_l\f$ for the (k+1)-th order Backward Differentiation Formula (BDF): \f$ f'(x_i) \approx \frac{1}{h} \sum_{l=0}^{k+1} c_l f(x_{i-l}) \f$.
     *   `poly_differentiation(int i, int l)`: Returns coefficients \f$c_{i,l}\f$ for polynomial differentiation: \f$ f'(x_i) \approx \frac{1}{h} \sum_{l=0}^k c_{i,l} f(x_l) \f$, where the polynomial passes through \f$f(x_0), ..., f(x_k)\f$.
     *   `poly_integration(int i, int j, int l)`: Returns weights \f$w_{i,j,l}\f$ for integrating a k-th order polynomial: \f$ \int_{x_i}^{x_j} f(x) dx = h \sum_{l=0}^k w_{i,j,l} f(x_l) \f$. Used for small intervals or specific integration schemes like Volterra equation solvers.
     *   `poly_interpolation(int alpha, int l)`: Returns coefficients for constructing an interpolating polynomial \f$ P(x) = \sum_{\alpha,l} (x/h)^\alpha c_{\alpha,l} f(x_l) \f$. (Primarily for completeness, as stated in comments).
     *   `rcorr(int m, int j, int l)`: Returns special weights \f$w_{m,j,l}\f$ for boundary corrections in Matsubara convolutions of the form \f$ \int_0^{\tau_m} A(\tau_m-\tau)B(\tau) d\tau \approx \Delta\tau \sum_{j,l=0}^k w_{m,j,l} A(\tau_j)B(\tau_l) \f$, when \f$m < k\f$.
   - **Example Integration Method:**
     *   `integrate(const std::vector<M>& f, int n)`: A template method to compute \f$\int_0^{x_n} f(x) dx\f$ using Gregory quadrature for a vector `f` of a generic type `M` (where `M=0` is defined).

**2. Static Factory Function `I(int k)`:**
   - `template<typename T> Integrator<T>& I(int k)`:
     *   Provides convenient access to cached (static) `Integrator` instances of a specific order `k`. This avoids repeated initialization of weights if the same order is used multiple times.
     *   Supports orders from 1 to 5 for the cache.

**3. Weight Initialization Functions (Internal):**
   - `read_poly_interpolation(int k, T* P)`
   - `read_poly_differentiation(int k, T* P)`
   - `read_poly_integration(int k, T* P)`
   - `read_gregory_weights(int k, T* w)`
   - `read_bd_weights(int k, T* w)`
   - `read_rcorr(int k, T* w)`
   - These functions are declared and are responsible for loading the hardcoded or pre-calculated numerical coefficients into the `Integrator` instance. Their definitions are likely in a corresponding `.cpp` file (e.g., `integration.cpp` if it were present and contained these, or compiled from another source).

## Important Variables/Constants

-   **`GREGORY_KMAX`**: A preprocessor define, set to 8. This is the maximum order supported for Gregory integration.
-   **`MAX_SOLVE_ORDER` (from `cntr_global_settings.hpp`)**: Often used as the default order `k` when an `Integrator` is requested or used by other modules (e.g., Dyson solvers, convolution routines). This value is typically between 3 and 5.

## Usage Examples

The `Integrator` class is not usually instantiated directly by the end-user of the `libcntr` library. Instead, it's used internally by higher-level functions that perform numerical time-stepping, differentiation, or convolution.

**Conceptual Internal Usage by other `libcntr` modules:**

```cpp
// Within a function in, for example, cntr_convolution_impl.hpp
// Get an integrator of a specific order (e.g., SolveOrder defined elsewhere)
integration::Integrator<double>& integrator = integration::I<double>(SolveOrder);

// Using Gregory weights for a convolution integral over 'm' steps:
// C[i] = sum_j w_ij * A[i-j] * B[j]
// where w_ij might be integrator.gregory_weights(m, j) or integrator.gregory_omega(j)
// depending on the part of the sum.

// Using BDF weights for a time derivative in a Dyson solver:
// dA_dt_at_n = (1.0/h) * sum_l integrator.bd_weights(l) * A[n-l];

// Using polynomial integration for a specific part of a Volterra equation:
// integral_val = sum_l integrator.poly_integration(idx1, idx2, l) * F[l];
```
The `SolveOrder` parameter, often seen in `dyson`, `convolution`, and `differentiation` functions, dictates the order `k` of the `Integrator` instance used internally.

## Dependencies and Interactions

-   **Standard C++ Library**: `<cmath>`, `<cassert>`, `<iostream>`, `<complex>`, `<vector>`, `<stdlib.h>`.
-   **`cntr_global_settings.hpp`**: Likely for `MAX_SOLVE_ORDER` and other global typedefs or constants.
-   **Other `libcntr` Modules**: The `Integrator` is a low-level utility. Its coefficients are crucial for:
    *   `cntr_convolution_*.hpp`: For performing numerical convolutions (time integrals).
    *   `cntr_differentiation_*.hpp`: For numerical differentiation of contour functions.
    *   `cntr_dyson_*.hpp`: For time-stepping in solving Dyson's equations, which involves both differentiation and convolution.
    *   `cntr_equilibrium_*.hpp`: Potentially for interpolating Hamiltonians in `green_from_H` if `SolveOrder` is passed.
-   The actual definitions of the `read_*` functions that load the numerical weights are expected to be in a compiled source file (e.g., `integration.cpp` if it were part of the core library source, or linked from a precompiled part of `libcntr`). The header only declares them.

The `integration::Integrator` class centralizes the provision of numerical coefficients, ensuring consistency and reusability across the `libcntr` library for tasks requiring finite difference or polynomial-based numerical methods.
