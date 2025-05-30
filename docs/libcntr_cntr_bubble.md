## Overview

The files `cntr_bubble_decl.hpp` and `cntr_bubble_impl.hpp` provide functionality within the `libcntr` library for calculating "bubble diagrams." These diagrams are fundamental in many-body physics and quantum field theory, typically representing convolutions of two Green's functions. This module supplies templated functions to compute two specific types of bubble diagrams, `Bubble1` and `Bubble2`, which differ in their time argument conventions for the constituent functions. The calculations are performed on the Kadanoff-Baym contour, handling retarded, lesser, advanced, greater, and Matsubara components appropriately.

## Key Components

The primary components are two sets of overloaded templated functions, `Bubble1` and `Bubble2`, declared in `cntr_bubble_decl.hpp` and implemented in `cntr_bubble_impl.hpp`.

**1. `Bubble1` functions:**
   - **Purpose:** Evaluate a bubble diagram \f$ C \f$ from two-time contour functions \f$ A, B \f$ according to the formula:
     \f$ C_{c_1,c_2}(t_1,t_2) = i \cdot A_{a_1,a_2}(t_1,t_2) \cdot B_{b_2,b_1}(t_2,t_1) \f$.
   - The calculation is performed at a given `tstp` (time step) for all contour components (retarded, lesser, left-mixing/time-vertical, and Matsubara).
   - **Key Template Parameters:**
     - `GGC`: Type of the resulting Green's function \f$ C \f$.
     - `GGA`: Type of the Green's function \f$ A \f$.
     - `GGB`: Type of the Green's function \f$ B \f$.
     (These types can be `herm_matrix`, `herm_matrix_timestep`, `herm_matrix_timestep_view`, or compatible types.)
   - **Main Overloads:**
     - `Bubble1(int tstp, GGC &C, int c1, int c2, GGA &A, GGA &Acc, int a1, int a2, GGB &B, GGB &Bcc, int b1, int b2)`:
       General version taking matrix indices (`c1,c2` for `C`; `a1,a2` for `A`; `b1,b2` for `B`) and explicit adjoint objects (`Acc` for `A`, `Bcc` for `B`).
     - `Bubble1(int tstp, GGC &C, int c1, int c2, GGA &A, int a1, int a2, GGB &B, int b1, int b2)`:
       Assumes \f$ A \f$ and \f$ B \f$ have Hermitian symmetry (i.e., their adjoints are themselves).
     - `Bubble1(int tstp, GGC &C, GGA &A, GGA &Acc, GGB &B, GGB &Bcc)`:
       Assumes \f$ A, B, C \f$ are \f$ 1 \times 1 \f$ matrices (scalar-like). Takes explicit adjoints.
     - `Bubble1(int tstp, GGC &C, GGA &A, GGB &B)`:
       Assumes \f$ A, B, C \f$ are \f$ 1 \times 1 \f$ matrices and \f$ A, B \f$ have Hermitian symmetry.
   - **Implementation Details (`cntr_bubble_impl.hpp`):**
     - Uses helper function `get_bubble_1_mat` for Matsubara components:
       \f$ C_{c_1,c_2}(\tau) = - A_{a_1,a_2}(\tau) \cdot B_{b_2,b_1}(-\tau) \f$
     - Uses helper function `get_bubble_1_timestep` for real-time components (retarded, lesser, tv). This involves careful handling of different Green's function components (gtr, les, ret, etc.) on the contour.
     - Wraps Green's function objects with `herm_matrix_timestep_view` internally if they are not already of that type, allowing generic application.

**2. `Bubble2` functions:**
   - **Purpose:** Evaluate a bubble diagram \f$ C \f$ from two-time contour functions \f$ A, B \f$ according to the formula:
     \f$ C_{c_1,c_2}(t_1,t_2) = i \cdot A_{a_1,a_2}(t_1,t_2) \cdot B_{b_1,b_2}(t_1,t_2) \f$.
   - Note the difference in indices for \f$ B \f$ compared to `Bubble1`.
   - Calculation details and template parameters are analogous to `Bubble1`.
   - **Main Overloads:** Similar to `Bubble1`, with versions for general matrix indices, Hermitian symmetry, and \f$ 1 \times 1 \f$ matrices.
     - `Bubble2(int tstp, GGC &C, int c1, int c2, GGA &A, GGA &Acc, int a1, int a2, GGB &B, GGB &Bcc, int b1, int b2)`
     - `Bubble2(int tstp, GGC &C, int c1, int c2, GGA &A, int a1, int a2, GGB &B, int b1, int b2)`
     - `Bubble2(int tstp, GGC &C, GGA &A, GGA &Acc, GGB &B, GGB &Bcc)`
     - `Bubble2(int tstp, GGC &C, GGA &A, GGB &B)`
   - **Implementation Details (`cntr_bubble_impl.hpp`):**
     - Uses helper function `get_bubble_2_mat` for Matsubara components:
       \f$ C_{c_1,c_2}(\tau) = - A_{a_1,a_2}(\tau) \cdot B_{b_1,b_2}(\tau) \f$
     - Uses helper function `get_bubble_2_timestep` for real-time components.
     - Also uses `herm_matrix_timestep_view` for genericity.

**Helper Functions (private, in `cntr_bubble_impl.hpp`):**
- `get_bubble_1_mat<T>(...)`: Calculates Matsubara component for `Bubble1`.
- `get_bubble_1_timestep<T>(...)`: Calculates real-time components for `Bubble1`.
- `get_bubble_2_mat<T>(...)`: Calculates Matsubara component for `Bubble2`.
- `get_bubble_2_timestep<T>(...)`: Calculates real-time components for `Bubble2`.
  These functions take raw pointers to complex data, sizes, indices, and other parameters to perform the low-level calculations.

## Important Variables/Constants

- **`CNTR_BUBBLE_DECL_H` / `CNTR_BUBBLE_IMPL_H`**: Include guards for the respective header files.
- **`ii` (std::complex<T>(0, 1.0))**: Represents the imaginary unit \f$ i \f$, used in the bubble calculation formulas. This is defined locally within implementation functions.
- **`tstp`**: An integer parameter passed to all `Bubble1` and `Bubble2` functions, representing the current time step at which the calculation is performed. If `tstp == -1`, only the Matsubara component is calculated.
- **`ntau`**: Number of Matsubara time steps, obtained from the Green's function objects.
- **`sigb` / `B.sig()`**: A sign factor associated with the \f$ B \f$ Green's function, used in Matsubara calculations for `Bubble1`.

## Usage Examples

The `Bubble1` and `Bubble2` functions are intended to be used within the `libcntr` framework where contour Green's function objects (of types like `herm_matrix_timestep` or compatible structures) are already set up.

A typical usage scenario would involve:
1. Initializing Green's function objects `A`, `B`, and `C` with appropriate sizes, time steps (`nt`, `ntau`), and other parameters.
2. Populating `A` and `B` with data from previous calculations or model definitions.
3. Calling one of the `Bubble` functions to compute `C`.

**Conceptual Example (Illustrative):**
```cpp
// Assuming G_type is a Green's function type compatible with the templates
// (e.g., cntr::herm_matrix_timestep<double>)

// Initialize Green's functions
G_type G_A, G_B, G_C; 
// ... (setup dimensions, time parameters for G_A, G_B, G_C)
// ... (fill G_A and G_B with data)

int current_timestep = 10; // Example time step
int a1=0, a2=0, b1=0, b2=0, c1=0, c2=0; // Example matrix indices

// Calculate Bubble1: C(t1,t2) = i * A(t1,t2) * B(t2,t1)
// Assuming G_A and G_B are Hermitian and 1x1
cntr::Bubble1(current_timestep, G_C, G_A, G_B);

// Or, for matrix versions with explicit adjoints:
// G_type G_A_adj, G_B_adj;
// ... (populate G_A_adj, G_B_adj)
// cntr::Bubble1(current_timestep, G_C, c1, c2, G_A, G_A_adj, a1, a2, G_B, G_B_adj, b1, b2);

// Calculate Bubble2: C(t1,t2) = i * A(t1,t2) * B(t1,t2)
// cntr::Bubble2(current_timestep, G_C, G_A, G_B); 
```
The user selects the appropriate overload based on whether the Green's functions are matrix-valued (requiring indices) or scalar-like (1x1), and whether explicit adjoint objects are provided or Hermitian symmetry is assumed. The `tstp` parameter controls the point up to which the time-dependent parts of the bubble are calculated.

## Dependencies and Interactions

- **`cntr_global_settings.hpp`**: Included in `cntr_bubble_decl.hpp`, suggesting it might provide global type definitions or constants used by the bubble functions.
- **`cntr_herm_matrix_timestep_view_decl.hpp`**: Included in `cntr_bubble_impl.hpp`. The `herm_matrix_timestep_view` class is extensively used in the implementation of `Bubble1` and `Bubble2` functions to provide a common interface/view over different Green's function data structures. This implies a strong dependency on this component of `libcntr`.
- **Green's Function Objects**: The bubble functions operate on generic Green's function objects (`GGC`, `GGA`, `GGB`). These objects are expected to provide methods like `ntau()`, `tstp_`, `size1()`, `matptr()`, `retptr()`, `tvptr()`, `lesptr()`, `sig()`, and a `scalar_type` typedef. This indicates a tight coupling with the Green's function representations within `libcntr`.
- **C++ Standard Library**: Uses `<complex>` for complex number arithmetic and `<cassert>` (implied by `assert`) for runtime checks.
- **Namespace `cntr`**: All declarations and implementations are within the `cntr` namespace.

These bubble calculation routines are likely used in higher-level algorithms within the Kaba-Neco library, such as solving Dyson equations, calculating self-energies, or computing response functions where such diagrammatic terms appear.
