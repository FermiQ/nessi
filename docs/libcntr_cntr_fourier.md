## Overview

The Fourier transform module, primarily defined in `libcntr/cntr/fourier.hpp` and `libcntr/cntr/fourier.cpp`, provides tools for performing accurate numerical Fourier transforms. These are essential for converting functions between time (or imaginary time \f$\tau\f$) and frequency (\f$\omega\f$) domains, a common operation in many-body physics, particularly when dealing with Green's functions and self-energies.

The module focuses on Discrete Fourier Transforms (DFT) with endpoint corrections (linear and cubic) to improve accuracy when approximating continuous Fourier integrals. It also features an adaptive DFT class (`adft_func`) that can refine the integration grid to achieve a desired level of precision.

## Key Components

All functionalities are within the `fourier` namespace.

**1. DFT Correction Functions (defined in `fourier.cpp`):**
   These functions calculate correction factors needed to improve the accuracy of simple DFTs when approximating a continuous Fourier integral \f$I(\omega) = \int_a^b dt\, e^{i \omega t} f(t)\f$. The corrections account for the behavior of the function at the endpoints of the integration interval. (See Press, Teukolsky, Vetterling, Flannery, "Numerical Recipes", Ch. 13.9).

   - **`get_dftcorr_linear(double th, double* corfac, std::complex<double>* endcor)`**:
     *   Calculates correction factors for linearly corrected DFT.
     *   `th`: \f$\theta = \omega h\f$ (product of frequency and grid spacing).
     *   `corfac`: (Output) The main correction factor \f$W(\theta)\f$.
     *   `endcor`: (Output) Boundary correction term \f$\alpha_0(\theta)\f$.
   - **`get_dftcorr_cubic(double th, double* corfac, std::complex<double>* endcor)`**:
     *   Calculates correction factors for cubically corrected DFT.
     *   `th`: \f$\theta = \omega h\f$.
     *   `corfac`: (Output) The main correction factor \f$W(\theta)\f$.
     *   `endcor`: (Output) Array of boundary correction terms \f$\alpha_j(\theta)\f$ for j=0,1,2,3.
   - **`complex_dftcor_cubic(double w, double delta, double a, double b, std::complex<double>* endpts, std::complex<double>* endcor, double* corfac)`**:
     *   Combines the cubic correction factors with the actual endpoint values of the function being transformed.
     *   `endpts`: Array containing function values at the start (\f$f_0, ..., f_3\f$) and end (\f$f_N, ..., f_{N-3}\f$) of the interval.
     *   `endcor`: (Output) The total boundary correction sum.
     *   `corfac`: (Output) The main correction factor \f$W(\theta=\omega\delta)\f$.
   - **`dft_cplx(double w, int n, double a, double b, std::complex<double>* f_values, std::complex<double>& result, std::complex<double>& error_estimate)`**:
     *   Computes the Fourier integral \f$\int_a^b dt\, e^{i \omega t} f(t)\f$ using cubically corrected DFT.
     *   `f_values`: Array of `n+1` function values on an equidistant grid over \f$[a,b]\f$.
     *   `result`: The calculated value of the Fourier integral.
     *   `error_estimate`: An estimate of the error, obtained by comparing the result with one computed using half the number of points.

**2. `adft_func` Class (Adaptive DFT):**
   - **Purpose**: Provides a class to compute Fourier transforms with adaptive refinement of the integration interval to achieve higher accuracy.
   - **Constructor**: `adft_func(int size, int nmax)` initializes internal buffers. `size` is the total number of sample points available, `nmax` is the maximum number of subintervals.
   - **Key Public-like Method**:
     *   `sample(double w, double a, double b, Function& fz, int n_pts_per_interval, int limit_max_intervals)`:
         *   This is the main routine to prepare the `adft_func` object for a transform at frequency `w`.
         *   It adaptively samples the function `fz` (which must be callable, e.g., `fz(x)`) over the interval \f$[a,b]\f$.
         *   The interval is recursively split if the error estimate from `dft_cplx` on a subinterval is too large.
         *   `n_pts_per_interval`: Number of points used for DFT on the finest subintervals (must be >= `ADFT_MINPTS` and even).
         *   `limit_max_intervals`: Maximum number of subintervals the original interval can be divided into.
   - **Key Method for Calculation**:
     *   `dft(double w, std::complex<double>& result, std::complex<double>& err)`:
         *   After `sample()` has been called for a specific `w`, this method computes the Fourier integral by summing the results of `dft_cplx` over all the adaptively determined subintervals.
         *   `result`: The final value of the Fourier integral.
         *   `err`: The sum of error estimates from subintervals.
   - **Internal Methods** (marked `@private` but crucial for operation):
     *   `init(...)`: Initializes the adaptive process with a single interval.
     *   `split(...)`: Splits a subinterval into two, re-sampling the function.
     *   `interval_data(int i)`: Returns a pointer to the sampled function values for subinterval `i`.

## Important Variables/Constants

-   **`PI`**: Defined as `3.14159265358979323846`.
-   **`ADFT_MINPTS`**: Defined as `16`. The minimum number of points an interval in `adft_func` can have.

## Usage Examples

The Fourier transform tools are likely used when transforming data between time/imaginary-time and frequency/Matsubara-frequency domains. The `adft_func` is particularly useful for functions that are not smooth or require high precision.

**Conceptual Usage of `adft_func`:**

```cpp
#include "cntr/fourier.hpp"
#include <complex>
#include <vector>
#include <cmath>

// Example function to be Fourier transformed
struct MyFunction {
    std::complex<double> operator()(double t) const {
        return std::exp(-t * t) * std::cos(10.0 * t); // Example: a Gaussian wavepacket
    }
};

int main() {
    fourier::adft_func adaptive_ft;
    MyFunction my_func;

    double omega_target = 10.0; // Frequency at which to compute the transform
    double t_start = -5.0;
    double t_end = 5.0;
    int points_per_interval = 32; // Must be >= 16 and even
    int max_intervals = 100;      // Max adaptive subdivisions

    // Prepare the adaptive grid for the target frequency
    adaptive_ft.sample(omega_target, t_start, t_end, my_func, points_per_interval, max_intervals);

    // Compute the Fourier transform at omega_target
    std::complex<double> ft_result, ft_error_est;
    adaptive_ft.dft(omega_target, ft_result, ft_error_est);

    // std::cout << "FT at omega = " << omega_target << " is " << ft_result 
    //           << " (error est: " << ft_error_est << ")" << std::endl;

    return 0;
}
```

The non-adaptive `dft_cplx` could be used directly if the function is smooth and a fixed grid is sufficient. The correction factor functions (`get_dftcorr_*`) are used internally by `dft_cplx` and `adft_func`, and also by `cntr_matsubara_impl.hpp` for its `matsubara_ft` function.

## Dependencies and Interactions

-   **Standard C++ Library**: `<cmath>`, `<cassert>`, `<iostream>`, `<complex>`, `<vector>`, `<stdlib.h>`.
-   **`cntr_matsubara_impl.hpp`**: The Matsubara Fourier transform utilities in `cntr_matsubara` (specifically `matsubara_ft`) use the `get_dftcorr_linear` and `get_dftcorr_cubic` functions from this Fourier module to apply endpoint corrections to their DFTs.
-   **No direct external library dependencies like FFTW are apparent in these files.** The DFT performed is a direct summation, potentially made accurate by endpoint corrections and adaptive interval splitting, rather than a Fast Fourier Transform algorithm.

This module provides robust tools for accurate DFT calculations, which are essential for spectral analysis and transformations in the context of `libcntr`'s focus on contour Green's functions.
