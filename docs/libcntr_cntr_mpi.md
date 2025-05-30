## Overview

The `libcntr` library incorporates Message Passing Interface (MPI) functionality to enable parallel computations, particularly for distributing data and workloads across multiple processes. Unlike some libraries that might have a dedicated MPI abstraction module (e.g., `cntr_mpi_decl.hpp`), `libcntr` integrates MPI capabilities more directly.

The primary configuration for MPI is managed within `libcntr/cntr/cntr_global_settings.hpp`. This file handles the conditional inclusion of MPI headers and defines a macro for the master process rank. Actual MPI communication calls (like broadcasts, reductions, send/receive) are typically embedded as member functions within the core data structure classes of the library, such as `cntr::herm_matrix` and `cntr::function`. These MPI-enabled methods are compiled if the `CNTR_USE_MPI` macro is defined.

## Key Components

The "MPI functionality" in `libcntr` is not a standalone module with its own set of wrapper classes or extensive utility functions. Instead, it's characterized by:

**1. Configuration via `cntr_global_settings.hpp`:**
   - **`CNTR_USE_MPI` Macro**:
     *   This preprocessor macro controls the inclusion of MPI features. If `CNTR_USE_MPI` is defined (e.g., via a compiler flag like `-DCNTR_USE_MPI`), it is then defined as `1` within `cntr_global_settings.hpp`.
     *   Code sections using MPI are typically guarded by `#if CNTR_USE_MPI == 1`.
   - **MPI Header Inclusion**:
     *   If `CNTR_USE_MPI == 1`, the standard MPI header `<mpi.h>` is included. This makes all standard MPI functions, types, and constants available to the `libcntr` source files.
   - **`MPI_MASTER` Macro**:
     *   Defined as `#define MPI_MASTER 0`. This constant provides a convenient way to refer to the rank of the primary MPI process (root process for collectives, etc.).

**2. MPI Operations within Core Classes:**
   - As observed in the documentation for classes like `cntr::herm_matrix` and `cntr::function`, MPI operations are implemented as member functions of these classes. Examples include:
     *   `herm_matrix::Bcast_timestep(int tstp, int root)`
     *   `herm_matrix::Reduce_timestep(int tstp, int root)`
     *   `herm_matrix::Send_timestep(int tstp, int dest, int tag)`
     *   `herm_matrix::Recv_timestep(int tstp, int root, int tag)`
     *   `function::Bcast_timestep(int tstp, int root)`
   - These methods internally use standard MPI calls (e.g., `MPI_Bcast`, `MPI_Reduce`, `MPI_Send`, `MPI_Recv`) to communicate data relevant to the specific object instance (e.g., a particular time-slice of a Green's function).

**3. No Dedicated MPI Wrapper Module:**
   - The specified files `libcntr/cntr/cntr_mpi_decl.hpp` and `libcntr/cntr/cntr_mpi_impl.hpp` were not found. This indicates that `libcntr` does not provide a separate abstraction layer or utility suite specifically for MPI operations beyond the basic setup in `cntr_global_settings.hpp`. Users of the library interact with MPI through the parallel methods of the main data classes.

## Important Variables/Constants

-   **`CNTR_USE_MPI`**: Preprocessor macro. Its definition enables MPI code paths.
-   **`MPI_MASTER`**: Defined as `0`, conventionally the rank of the master process in MPI communications.
-   Standard MPI constants and types from `<mpi.h>` become available when `CNTR_USE_MPI` is active (e.g., `MPI_COMM_WORLD`, `MPI_DOUBLE_COMPLEX`, `MPI_INT`).

## Usage Examples

Since there isn't a separate MPI module, usage examples involve calling MPI-enabled methods of other `libcntr` classes.

**Conceptual Usage (within a hypothetical parallel application using `libcntr`):**

```cpp
#include "cntr/cntr_global_settings.hpp" // This will include mpi.h if CNTR_USE_MPI is set
#include "cntr/cntr_herm_matrix_decl.hpp"

#if CNTR_USE_MPI == 1
void run_parallel_simulation() {
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    cntr::herm_matrix<double> G_global;
    int current_timestep = 5;

    if (world_rank == MPI_MASTER) {
        // Initialize or compute G_global on the master process
        // ...
    }

    // Broadcast a specific timestep of G_global from master to all other processes
    // Assuming G_global is properly sized on all processes before the call.
    // The Bcast_timestep method itself handles the MPI_Bcast call.
    G_global.Bcast_timestep(current_timestep, MPI_MASTER);

    // ... (Each process now has G_global's data for current_timestep) ...

    // Example: All processes compute some partial result and reduce it
    cntr::herm_matrix<double> G_local_sum_contrib = G_global; // just for structure
    // ... perform local computation into G_local_sum_contrib for current_timestep ...
    
    // Reduce (sum) the timestep data from all processes to the master
    // The Reduce_timestep method would internally call MPI_Reduce
    G_local_sum_contrib.Reduce_timestep(current_timestep, MPI_MASTER);

    if (world_rank == MPI_MASTER) {
        // G_local_sum_contrib on master now holds the sum for current_timestep
        // ...
    }
}
#endif

int main(int argc, char* argv[]) {
#if CNTR_USE_MPI == 1
    MPI_Init(&argc, &argv);
#endif

    // ... application logic ...
#if CNTR_USE_MPI == 1
    // run_parallel_simulation();
#endif

#if CNTR_USE_MPI == 1
    MPI_Finalize();
#endif
    return 0;
}
```
The application itself is responsible for initializing and finalizing MPI. `libcntr` provides MPI-aware methods within its classes that can then be used in a parallel context.

## Dependencies and Interactions

-   **MPI Library**: The core dependency. The system must have an MPI implementation installed (e.g., OpenMPI, MPICH), and `libcntr` must be compiled with `CNTR_USE_MPI` defined and linked against the MPI library.
-   **`cntr_global_settings.hpp`**: Central point for enabling MPI via the `CNTR_USE_MPI` macro and including `<mpi.h>`.
-   **Core Data Structure Classes (e.g., `cntr::herm_matrix`, `cntr::function`)**: These classes contain member functions that implement MPI communication patterns for their data (e.g., broadcasting or reducing a time-slice).
-   **Build System**: The build system (e.g., CMake) would be responsible for detecting MPI and setting the `CNTR_USE_MPI` definition when compiling `libcntr` and any applications using it.

In summary, MPI support in `libcntr` is achieved by conditional compilation of MPI-specific code within its main data handling classes, rather than through a separate MPI abstraction layer module. Configuration is minimal and centralized in `cntr_global_settings.hpp`.
