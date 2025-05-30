## Overview

The files `cntr_herm_matrix_decl.hpp` and `cntr_herm_matrix_impl.hpp` define the `cntr::herm_matrix<T>` class, a central data structure in the `libcntr` library. This class is meticulously designed to represent two-time physical quantities on the Kadanoff-Baym contour, such as Green's functions \f$ G(t,t') \f$ or self-energies \f$ \Sigma(t,t') \f$, that possess Hermitian symmetry.

The Kadanoff-Baym contour includes an imaginary time branch (Matsubara) and real-time branches (Keldysh). `herm_matrix` efficiently stores only the non-redundant components of these two-time functions by exploiting their inherent Hermitian and contour symmetries. It provides a comprehensive suite of methods for construction, access, manipulation, and I/O of these complex objects.

## Key Components

**`cntr::herm_matrix<T>` Class:**

-   **Template Parameter `T`:** Typically `double` or `float`, indicating the precision of the `std::complex<T>` data.

-   **Contour Components Storage:**
    The class stores the following essential components of a two-time function \f$ C(t,t') \f$:
    *   **Matsubara (`mat_`)**: \f$ C^M(\tau_k) \f$ for discrete imaginary times \f$ \tau_k \f$. Stored as `mat_[(ntau_+1)*element_size_]`.
    *   **Retarded (`ret_`)**: \f$ C^R(t_i, t_j) \f$ for real times \f$ t_j \le t_i \f$. Stored in a packed format for the lower triangle (including diagonal) as `ret_[((nt_+1)*(nt_+2))/2*element_size_]`.
    *   **Lesser (`les_`)**: \f$ C^<(t_i, t_j) \f$ for real times \f$ t_i \le t_j \f$. Stored similarly to the retarded component for the upper triangle (including diagonal) as `les_[((nt_+1)*(nt_+2))/2*element_size_]`.
    *   **Time-Vertical / Left-Mixing (`tv_`)**: \f$ C^\rceil(t_i, \tau_k) \f$ for real time \f$ t_i \f$ and imaginary time \f$ \tau_k \f$. Stored as `tv_[(nt_+1)*(ntau_+1)*element_size_]`.

-   **Symmetry Utilization:**
    Other Keldysh components (e.g., Advanced \f$C^A\f$, Greater \f$C^>\f$, Right-Mixing \f$C^\lceil\f$) are derived from these stored components using established symmetry relations, such as:
    *   \f$ C^A(t_j, t_i) = (C^R(t_i, t_j))^\ddagger \f$
    *   \f$ C^>(t_i, t_j) = C^R(t_i, t_j) - C^<(t_i, t_j) + C^A(t_i, t_j) \f$ (Langreth rules) or other direct relations.
    *   \f$ C^<(t_j, t_i) = -\eta (C^<(t_i, t_j))^\ddagger \f$ (for \f$i \ne j\f$, specific to Hermitian definition)
    *   \f$ C^\lceil(\tau_k, t_i) = -\eta (C^\rceil(t_i, \beta - \tau_k))^\ddagger \f$

-   **Key Member Variables:**
    *   `les_`, `ret_`, `tv_`, `mat_`: Raw `std::complex<T>*` pointers to the data arrays for each stored component.
    *   `nt_`: Integer, number of real-time steps.
    *   `ntau_`: Integer, number of imaginary-time (Matsubara) steps.
    *   `size1_`, `size2_`: Integers, matrix dimensions at each (t,t') point. Typically `size1_ == size2_`.
    *   `element_size_`: Integer, `size1_ * size2_`.
    *   `sig_`: Integer, particle statistics flag: +1 for bosons, -1 for fermions.

-   **Constructors & Destructor:**
    *   Default constructor, parameterized constructors (`nt, ntau, size1, sig` or `nt, ntau, size1, size2, sig`), copy constructor, and C++11 move constructor/assignment.
    *   The destructor handles deallocation of all component data arrays.

-   **Resizing and Clearing:**
    *   `resize_discard(nt, ntau, size1)`: Resizes and discards existing data.
    *   `resize_nt(nt)`: Resizes only the real-time dimension, attempting to preserve existing data where possible.
    *   `resize(nt, ntau, size1)`: General resize, may discard data if dimensions change incompatibly.
    *   `clear()`: Sets all stored data to zero.

-   **Data Access (Getters):**
    *   `get_les(i, j, Matrix& M) const`, `get_ret(i, j, Matrix& M) const`, etc.: Retrieve specific components into a provided matrix object `M` (e.g., an Eigen matrix). These methods correctly apply symmetry relations to construct any requested component from the stored ones.
    *   Scalar overloads `get_les(i, j, cplx& x) const`, etc., for 1x1 matrices.
    *   `density_matrix(int tstp, Matrix& M) const`: Computes the density matrix \f$\rho(t) = i \eta C^<(t,t)\f$ or \f$\rho = -C^M(\beta)\f$.

-   **Data Modification (Setters):**
    *   `set_les(i, j, Matrix& M)`, `set_ret(i, j, Matrix& M)`, etc.: Set specific components from a matrix `M`.
    *   `set_mat_herm()`: Enforces \f$C^M(\tau) = (C^M(\tau))^\ddagger\f$ for all \f$\tau\f$.

-   **Time-Step Operations:**
    *   `set_timestep_zero(int tstp)`: Sets all components for a given real time-step `tstp` (or Matsubara if `tstp == -1`) to zero.
    *   `set_timestep(int tstp, herm_matrix& g1)` / `set_timestep(int tstp, herm_matrix_timestep<T>& ts)`: Copies data for a specific `tstp` from another `herm_matrix` or a `herm_matrix_timestep` object.
    *   `get_timestep(int tstp, herm_matrix_timestep<T>& ts) const`: Extracts data for `tstp` into a `herm_matrix_timestep` object.
    *   `incr_timestep(...)`: Adds data from another `herm_matrix` or `herm_matrix_timestep` to `tstp` of this object, optionally with a scaling factor.

-   **Matrix and Element Operations:**
    *   `set_matrixelement(...)`: Modifies individual matrix elements within the contour object for specific times or all times.
    *   `set_submatrix(...)`: Modifies a block of matrix elements.
    *   `left_multiply(int tstp, function<T>& ft, ...)` / `right_multiply(...)`: Multiplies a time-slice of the `herm_matrix` by a `cntr::function<T>` object from the left/right. E.g., \f$C(t_{tstp}, t') \leftarrow F(t_{tstp}) C(t_{tstp}, t')\f$.
    *   `left_multiply_hermconj(...)` / `right_multiply_hermconj(...)`: Similar, but multiplies with the Hermitian conjugate of the `cntr::function<T>`.
    *   `smul(int tstp, T weight)` / `smul(int tstp, cplx weight)`: Scales a time-slice by a scalar.

-   **File I/O:**
    *   Text format: `print_to_file()`, `read_from_file()`.
    *   HDF5 format (if `CNTR_USE_HDF5 == 1`): `write_to_hdf5()`, `read_from_hdf5()`. Includes specialized HDF5 I/O for "slices" and "average-relative time" (Wigner) coordinates (`write_to_hdf5_slices`, `write_to_hdf5_tavtrel`).

-   **MPI Support (if `CNTR_USE_MPI == 1`):**
    *   `Reduce_timestep()`, `Bcast_timestep()`, `Send_timestep()`, `Recv_timestep()`: For parallel communication of time-slice data.

**`herm_matrix_timestep_view<T>` and `herm_matrix_timestep<T>`:**
Although `herm_matrix_timestep_view_impl.hpp` is included, indicating its use, `herm_matrix_timestep<T>` is the type more directly visible in `herm_matrix`'s public interface for timestep-specific operations. `herm_matrix_timestep<T>` likely represents or provides a view of all components (retarded, lesser, tv) associated with a single real time index `t` (or the Matsubara component if `t=-1`), facilitating efficient operations on individual time slices without manipulating the entire `herm_matrix` object directly or copying large amounts of data unnecessarily.

## Important Variables/Constants

-   **`CNTR_HERM_MATRIX_DECL_H` / `CNTR_HERM_MATRIX_IMPL_H`**: Include guards.
-   **`CNTR_USE_HDF5` / `CNTR_USE_MPI`**: Preprocessor macros enabling HDF5 and MPI features.
-   **Member `sig_`**: Defines particle statistics (+1 Bose, -1 Fermi), crucial for symmetry relations.

## Usage Examples

`cntr::herm_matrix<T>` is the primary container for Green's functions and self-energies in `libcntr`.

**Conceptual Usage:**

```cpp
#include "cntr/cntr_herm_matrix_decl.hpp"
#include "cntr/cntr_function_decl.hpp" // For multiplication examples
// Assuming Eigen is available if using Matrix template arguments directly
// #include <Eigen/Dense>

using GMatrix = cntr::herm_matrix<double>;
using CFunction = cntr::function<double>;

// Parameters
int nt_real = 10;
int ntau_matsubara = 64;
int matrix_dimension = 1; // Scalar Green's function
int statistics = cntr::FERMI; // or cntr::BOSE

// Construction
GMatrix G(nt_real, ntau_matsubara, matrix_dimension, statistics);
G.clear(); // Initialize to zero

// Setting a Matsubara value (e.g., from an analytical formula or previous calculation)
Eigen::MatrixXcd mat_val(matrix_dimension, matrix_dimension);
mat_val(0,0) = std::complex<double>(1.0, 0.5);
if (matrix_dimension == 1) G.set_mat(5, mat_val(0,0)); // Scalar version
else G.set_mat(5, mat_val);                         // Matrix version

// Getting a retarded component value
Eigen::MatrixXcd ret_val_check(matrix_dimension, matrix_dimension);
G.get_ret(3, 2, ret_val_check); // Get G_R(t=3*h, t'=2*h)

// Example of time-step operation using herm_matrix_timestep
cntr::herm_matrix_timestep<double> G_slice_t5;
G.get_timestep(5, G_slice_t5); // Get all components for t=5
// ... operate on G_slice_t5 ...
// G.set_timestep(5, G_slice_t5); // Put it back if modified

// Example: Multiply by a time-dependent function H(t)
CFunction H_func(nt_real, matrix_dimension);
// ... initialize H_func ...
for (int t = -1; t <= nt_real; ++t) {
    G.left_multiply(t, H_func); // G(t,t') = H(t)G(t,t')
}

// File I/O
G.print_to_file("green_function.dat");
GMatrix G_loaded;
G_loaded.read_from_file("green_function.dat");
```

## Dependencies and Interactions

-   **`cntr_global_settings.hpp`**: Provides global constants and type definitions.
-   **`cntr_elements.hpp`**: For low-level, element-wise matrix operations (`element_set`, `element_mult`, `element_conj`, etc.).
-   **`cntr_function_decl.hpp`**: Defines `cntr::function<T>`, used in multiplication operations.
-   **`cntr_herm_matrix_timestep_decl.hpp` / `cntr_herm_matrix_timestep_view_impl.hpp`**: Define and implement views or dedicated objects for single time-slices, used by `herm_matrix` for efficient timestep-wise operations.
-   **Eigen Library**: While not a direct include in the header, the `get_X`/`set_X` methods taking `Matrix& M` arguments strongly imply that `M` is often an Eigen matrix type (like `cdmatrix`).
-   **HDF5 and MPI libraries**: Conditional dependencies based on `CNTR_USE_HDF5` and `CNTR_USE_MPI`.

The `herm_matrix` class is fundamental for representing the central objects (Green's functions, self-energies) in simulations based on Kadanoff-Baym equations or related contour-ordered formalisms. Higher-level modules, such as Dyson equation solvers or convolution routines, operate extensively on `herm_matrix` objects.
