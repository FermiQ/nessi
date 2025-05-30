## Overview

The files `cntr_function_decl.hpp` and `cntr_function_impl.hpp` define the `cntr::function<T>` class, a versatile C++ template class within the `libcntr` library. This class is designed to represent generic time-dependent functions, \f$ F(t) \f$, which can be scalar or matrix-valued. The time argument \f$ t \f$ can be real, or a special index `t = -1` is used to denote data associated with the Matsubara (imaginary time) axis, typically representing equilibrium or initial state values.

This `function` class serves as a fundamental building block for representing various physical quantities that depend on a single time variable, such as time-dependent Hamiltonians, local potentials, fields, or time-local components of self-energies. It provides storage, accessors, and basic mathematical operations for these time-dependent objects.

## Key Components

**`cntr::function<T>` Class:**

-   **Template Parameter `T`:** Typically `double` or `float`, defining the precision of the complex data stored.

-   **Storage and Structure:**
    *   Internally stores data as a contiguous array of `std::complex<T>`.
    *   The data represents a series of "time slices". Each slice corresponds to a specific time point.
    *   Time points include `t = -1` (for Matsubara/initial data) and `t = 0, 1, ..., nt_` (for real-time data).
    *   At each time point, the function can be a scalar or a matrix of dimensions `size1_` x `size2_`.

-   **Key Member Variables:**
    *   `data_`: A raw pointer `std::complex<T>*` to the allocated memory holding the function data.
    *   `nt_`: Integer representing the maximum index for real-time steps. The total number of time slices is `nt_ + 2`.
    *   `size1_`: Integer, number of rows in the matrix representation at each time point.
    *   `size2_`: Integer, number of columns in the matrix representation at each time point. For scalar functions, `size1_ = size2_ = 1`.
    *   `element_size_`: Integer, `size1_ * size2_`, the number of complex values per time slice.
    *   `total_size_`: Integer, total number of complex values stored: `(nt_ + 2) * element_size_`.

-   **Constructors & Destructor:**
    *   `function()`: Default constructor (data is null, sizes are zero).
    *   `function(int nt, int size1 = 1)`: Constructs a function for `nt` real time steps, with each time point holding a square matrix of size `size1`x`size1`.
    *   `function(int nt, int size1, int size2)`: Constructs for general `size1`x`size2` matrices.
    *   `function(const function& f)`: Copy constructor.
    *   `function(function&& f) noexcept`: Move constructor (C++11).
    *   `~function()`: Destructor, deallocates `data_`.

-   **Assignment Operators:**
    *   `operator=(const function& f)`: Copy assignment.
    *   `operator=(function&& f) noexcept`: Move assignment (C++11).

-   **Accessors and Mutators:**
    *   `ptr(int t)`: Returns a raw pointer `std::complex<T>*` to the beginning of the data for time slice `t`.
    *   `resize(int nt, int size1)`: Resizes the function (square matrix assumed). Data is lost.
    *   `set_zero()`: Sets all elements of the function at all times to zero.
    *   `set_constant(std::complex<T>* f0)`: Sets the function at all time points to the constant matrix provided by `f0`.
    *   `set_constant(EigenMatrix& M)`: Sets the function at all time points to the constant Eigen matrix `M`.
    *   `set_value(int tstp, EigenMatrix& M)`: Sets the matrix value at time `tstp` from an Eigen matrix `M`.
    *   `set_value(int tstp, std::complex<T> x)`: Sets the scalar value at time `tstp` (for 1x1 functions).
    *   `get_value(int tstp, EigenMatrix& M) const`: Retrieves the matrix value at time `tstp` into an Eigen matrix `M`.
    *   `operator[](int t)`: Direct access to the complex value at time `t`; only safe for scalar functions (1x1).

-   **Element-wise (in time) Mathematical Operations:**
    *   `smul(T weight)`: Multiplies the entire function (all time points, all matrix elements) by a scalar `weight`.
    *   `left_multiply(std::complex<T>* f, T weight = 1.0)` or `left_multiply(function<T>& f, T weight = 1.0)`: Performs \f$ F_{this}(t) \leftarrow M(t) \cdot F_{this}(t) \f$ for each time `t`. \f$ M(t) \f$ is provided by the argument function `f`.
    *   `right_multiply(std::complex<T>* f, T weight = 1.0)` or `right_multiply(function<T>& f, T weight = 1.0)`: Performs \f$ F_{this}(t) \leftarrow F_{this}(t) \cdot M(t) \f$.
    *   `incr(function<T>& g, T weight = 1.0)`: Adds function `g` (scaled by `weight`) to this function: \f$ F_{this}(t) \leftarrow F_{this}(t) + \text{weight} \cdot g(t) \f$.

-   **Matrix Element Manipulation:**
    *   `set_matrixelement(int tstp, int i1, int i2, EigenMatrix& M, int j1, int j2)`: Sets the element `(i1,i2)` of the function at time `tstp` to element `(j1,j2)` of Eigen matrix `M`.
    *   `set_matrixelement(int i1, int i2, function& g, int j1, int j2)`: Sets element `(i1,i2)` of this function (at all times) from element `(j1,j2)` of function `g` (at corresponding times).
    *   `get_matrixelement(int i1, int i2, function<T>& g)`: Extracts element `(i1,i2)` of this function (at all times) into function `g` (which becomes a scalar function).

-   **File I/O:**
    *   `print_to_file(const char* file, int precision = 16) const`: Writes the function data to a human-readable text file.
    *   `read_from_file(const char* file)` / `read_from_file(int nt1, const char* file)`: Reads function data from a text file, optionally up to `nt1` time steps.
    *   HDF5 Support (if `CNTR_USE_HDF5 == 1`):
        *   `write_to_hdf5(...)`: Various overloads to write data to an HDF5 file/group.
        *   `read_from_hdf5(...)`: Various overloads to read data from an HDF5 file/group.

-   **MPI Support (if `CNTR_USE_MPI == 1`):**
    *   `Bcast_timestep(int tstp, int root)`: Broadcasts the data for a specific time step `tstp` from a root process to all other MPI processes.

The class `herm_function` was mentioned in the subtask description but is not defined in the provided header files. It's possible it's a type alias or a derived class intended to enforce Hermiticity, but documentation is based on the provided `cntr::function` class.

## Important Variables/Constants

-   **`CNTR_FUNCTION_DECL_H` / `CNTR_FUNCTION_IMPL_H`**: Include guards.
-   **`CNTR_USE_HDF5`**: Preprocessor macro to enable HDF5 I/O.
-   **`CNTR_USE_MPI`**: Preprocessor macro to enable MPI functionalities.

## Usage Examples

The `cntr::function<T>` class is used to represent various single-time physical quantities in simulations.

**Conceptual Usage:**

```cpp
#include "cntr/cntr_function_decl.hpp"
// Assuming Eigen is available for matrix operations if needed
// #include <Eigen/Dense> // For cdmatrix example

using CFunction = cntr::function<double>; // double-precision complex function

// Example 1: Scalar time-dependent field
int num_timesteps = 100;
CFunction scalar_field(num_timesteps, 1); // size1=1 for scalar

scalar_field.set_value( -1, std::complex<double>(1.0, 0.0) ); // Matsubara/initial value
for (int t = 0; t <= num_timesteps; ++t) {
    scalar_field.set_value(t, std::complex<double>(cos(0.1 * t), sin(0.1 * t)));
}
// Access value (for scalar function)
std::complex<double> val_at_t5 = scalar_field[5];


// Example 2: Matrix-valued time-dependent Hamiltonian H(t)
int matrix_dim = 2;
CFunction H_t(num_timesteps, matrix_dim);

Eigen::MatrixXcd H_matrix_t(matrix_dim, matrix_dim);
H_matrix_t << 1.0, 0.5,
              0.5, 2.0;
H_t.set_value(-1, H_matrix_t); // Initial Hamiltonian for Matsubara reference

for (int t = 0; t <= num_timesteps; ++t) {
    H_matrix_t(0,1) = 0.5 * cos(0.05 * t);
    H_matrix_t(1,0) = 0.5 * cos(0.05 * t); // Assuming Hermitian
    H_t.set_value(t, H_matrix_t);
}

// Example 3: Scale the Hamiltonian
H_t.smul(0.5);

// Example 4: Add another function (e.g., a perturbation)
CFunction V_t(num_timesteps, matrix_dim);
// ... initialize V_t ...
H_t.incr(V_t, 1.0); // H_t = H_t + V_t

// Example 5: File I/O
H_t.print_to_file("hamiltonian_data.txt");
CFunction H_t_loaded;
H_t_loaded.read_from_file("hamiltonian_data.txt");

#ifdef CNTR_USE_HDF5
// H_t.write_to_hdf5("sim_output.h5", "dynamic_hamiltonian");
#endif
```

## Dependencies and Interactions

-   **`cntr_global_settings.hpp`**: Provides global settings like `CNTR_USE_HDF5`, `CNTR_USE_MPI`.
-   **`linalg.hpp`**: Likely provides `cdmatrix` (presumably an Eigen matrix typedef) and possibly other linear algebra utilities used in methods like `set_value`, `get_value`.
-   **`cntr_elements.hpp`**: Provides low-level element-wise operations (`element_set`, `element_mult`, `element_smul`) used in the implementation of mathematical operations for the `function` class.
-   **Standard C++ Library**: `<complex>`, `<vector>` (implicitly), `<fstream>`, `<string>`, `<cstring>` (for `memcpy`, `memset`), `<stdexcept>`/`<cstdlib>` (for `assert`, `abort`).
-   **HDF5 Library**: If `CNTR_USE_HDF5 == 1`, it links against the HDF5 C library for HDF5 I/O operations. Header `hdf5_io.hpp` (not shown but implied) would provide helper functions like `store_int_attribute_to_hid`.
-   **MPI Library**: If `CNTR_USE_MPI == 1`, it links against an MPI library for parallel communication.

The `cntr::function<T>` class is a foundational data structure in `libcntr`, used extensively by higher-level modules like Dyson equation solvers (`cntr_dyson_*.hpp`) or equilibrium Green's function constructors (`cntr_equilibrium_*.hpp`) to represent time-dependent Hamiltonians or local self-energy contributions.
