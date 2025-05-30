## Overview

The linear algebra module within `libcntr`, primarily defined by `libcntr/cntr/eigen_typedef.h` and `libcntr/cntr/linalg_eigen.cpp`, provides a set of C-style wrapper functions for common linear algebra operations. These functions internally leverage the Eigen C++ template library for numerical computation. The module also defines convenient typedefs for frequently used Eigen matrix and vector types.

The main purpose of this module is to offer a simplified, procedural interface to Eigen's capabilities, potentially for easier integration within the broader `libcntr` codebase or to abstract direct Eigen dependencies from other parts of the library. Operations include solving linear systems, matrix inversion, and eigenvalue problems for Hermitian matrices.

## Key Components

**1. Eigen Typedefs (`eigen_typedef.h`):**
   This header includes `<eigen3/Eigen/Dense>` and defines several type aliases for common Eigen dynamic-size matrix and vector types:
   -   **Vectors:**
        *   `ivector`: Alias for `Eigen::VectorXi` (vector of integers).
        *   `dvector`: Alias for `Eigen::VectorXd` (vector of doubles).
        *   `cdvector`: Alias for `Eigen::VectorXcd` (vector of `std::complex<double>`).
   -   **Matrices:**
        *   `imatrix`: Alias for `Eigen::MatrixXi` (matrix of integers).
        *   `dmatrix`: Alias for `Eigen::MatrixXd` (matrix of doubles).
        *   `cdmatrix`: Alias for `Eigen::MatrixXcd` (matrix of `std::complex<double>`).
   These typedefs are also brought into the global scope within `eigen_typedef.h` using `using Eigen::...`.

**2. Linear Algebra Functions (`linalg_eigen.cpp` within `namespace linalg`):**
   These functions typically take raw C-style pointers (`void*`, `double*`, `std::complex<double>*`) and dimensions as input, map them to Eigen objects, perform the computation using Eigen, and then copy results back if necessary.

   -   **Data Conversion Functions:**
        *   `set_cdmatrix(int n, void* a, cdmatrix& A)`: Maps raw complex data `a` to an Eigen matrix `A` (\f$n \times n\f$). Uses `Eigen::Map` and transposes.
        *   `set_dmatrix(int n, void* a, dmatrix& A)`: Maps raw double data `a` to an Eigen matrix `A` (\f$n \times n\f$). Uses `Eigen::Map` and transposes.
        *   `set_cdvector(int n, void* a, cdvector& A)`: Maps raw complex data `a` to an Eigen vector `A` (length `n`).
        *   `set_dvector(int n, void* a, dvector& A)`: Maps raw double data `a` to an Eigen vector `A` (length `n`).
        *   `get_cdmatrix(int n, void* a, cdmatrix& A)`: Copies data from Eigen matrix `A` to raw pointer `a`.
        *   `get_dmatrix(int n, void* a, const dmatrix& A)`: Copies data from Eigen matrix `A` to raw pointer `a`.
        *   `get_cdvector(int n, void* a, cdvector& A)`: Copies data from Eigen vector `A` to raw pointer `a`.
        *   `get_dvector(int n, void* a, dvector& A)`: Copies data from Eigen vector `A` to raw pointer `a`.

   -   **Linear System Solvers:**
        *   `cplx_sq_solve(void* a, void* b, void* x, int n)`: Solves \f$ Ax=b \f$ for complex square matrices \f$A\f$ (\f$n \times n\f$). Uses `Eigen::FullPivLU`.
        *   `cplx_sq_solve_many(void* a, void* b, void* x, int n, int d)`: Solves `d` systems \f$ Ax_i=b_i \f$ where \f$A\f$ is \f$n \times n\f$.
        *   `real_sq_solve(double* a, double* b, double* x, int n)`: Solves \f$ Ax=b \f$ for real square matrices.

   -   **Matrix Inversion:**
        *   `cplx_matrix_inverse(void* a, void* x, int n)`: Computes the inverse of a complex \f$n \times n\f$ matrix \f$A\f$, storing it in \f$X\f$. Uses `Eigen::FullPivLU`.
        *   `linalg_matrix_inverse(double* a, double* x, int n)`: Computes the inverse of a real \f$n \times n\f$ matrix. (Likely intended as `real_matrix_inverse`).

   -   **Eigenvalue Problems:**
        *   `eigen_hermv(int n, std::complex<double>* A_ptr, double* eval_ptr, std::complex<double>* evec_ptr)`: Computes eigenvalues (`eval_ptr`) and eigenvectors (`evec_ptr`) of a complex Hermitian \f$n \times n\f$ matrix `A_ptr`. Uses `Eigen::SelfAdjointEigenSolver`.

   -   **Commented-out Functionality:**
        *   `QR_decomposition`: A QR decomposition function was commented out in `linalg_eigen.cpp`.

   *(Note: The declarations for these functions would reside in `cntr_linalg_decl.hpp` or a similar header, which was not found but is implied by usage.)*

## Important Variables/Constants

-   **Type Aliases (from `eigen_typedef.h`):** `cdmatrix`, `dmatrix`, `cdvector`, `dvector`, `imatrix`, `ivector` are important type definitions used throughout `libcntr` for Eigen matrices and vectors.
-   **`cdouble` (from `cntr_global_settings.hpp` via `linalg_eigen.cpp`):** Likely a typedef for `std::complex<double>`.

## Usage Examples

These linear algebra functions are likely used internally by other modules in `libcntr` that require solving equations or matrix manipulations. For instance, `cntr::herm_matrix` methods might use these for operations on their underlying data, or Dyson solvers might use `cplx_sq_solve`.

**Conceptual Internal Usage:**

```cpp
// Assume data for matrix A and vector b are in raw C-style arrays
// double* a_data;
// std::complex<double>* b_data_cplx;
// std::complex<double>* x_solution_cplx;
// int n_dim = /* matrix dimension */;

// Example: Solve a complex linear system Ax = b
// linalg::cplx_sq_solve(a_data_cplx, b_data_cplx, x_solution_cplx, n_dim);

// Example: Invert a real matrix
// double* inv_a_data;
// linalg::linalg_matrix_inverse(a_data, inv_a_data, n_dim);

// Example: Eigenvalue decomposition of a Hermitian matrix
// std::complex<double>* herm_matrix_data;
// double* eigenvalues;
// std::complex<double>* eigenvectors_data;
// linalg::eigen_hermv(n_dim, herm_matrix_data, eigenvalues, eigenvectors_data);

// Using the typedefs for Eigen matrices directly in other C++ code:
// #include "eigen_typedef.h" // Or "cntr_linalg_decl.hpp" which should include it
// ...
// cdmatrix my_eigen_matrix(rows, cols);
// my_eigen_matrix(0,0) = std::complex<double>(1.0, 2.0);
// dvector my_eigen_vector(rows);
```

The `set_*` and `get_*` functions facilitate bridging between raw pointer data management (perhaps used in other parts of `libcntr` or for performance) and Eigen's object-oriented matrix/vector types.

## Dependencies and Interactions

-   **Eigen3 Library**: This module is fundamentally dependent on the Eigen3 library (`<eigen3/Eigen/Dense>`). All core computations are delegated to Eigen.
-   **`eigen_typedef.h`**: Provides the essential typedefs for Eigen matrices and vectors. Any file using these typedefs or the `linalg` functions would need to include the relevant header (presumably `cntr_linalg_decl.hpp` which in turn includes `eigen_typedef.h`).
-   **`cntr_global_settings.hpp`**: Included in `linalg_eigen.cpp`, likely for global type definitions such as `cdouble`.
-   **Other `libcntr` Modules**:
    *   The typedefs (`cdmatrix`, `dvector`, etc.) are widely used in other `libcntr` classes like `cntr::function`, `cntr::herm_matrix`, and within various computation routines (e.g., `cntr_equilibrium_impl.hpp` uses `cdmatrix` for Hamiltonian representations).
    *   The solver functions (`cplx_sq_solve`, `eigen_hermv`) are likely called by higher-level algorithms within the library that require solving linear systems or diagonalizing matrices (e.g., in Dyson solvers or certain equilibrium calculations).

This module serves as a C++ wrapper around Eigen, providing basic linear algebra capabilities through a C-style functional API and defining standard matrix/vector type aliases for use throughout the `libcntr` project. The use of `Eigen::Map` in `set_*` functions is an efficient way to handle raw data with Eigen, while `get_*` functions currently involve manual copying.
