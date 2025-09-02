# GEMINI.md - HypergeometricODEs

## Project Overview

This project, `HypergeometricODEs`, is a Julia-based research project focused on the numerical computation of hypergeometric functions, particularly the $_2F_1$ function. It explores and implements different methods for this purpose, including:

*   **Taylor Series Method:** A core method implemented in `Taylor.jl` for solving the underlying ordinary differential equation (ODE) of the hypergeometric function.
*   **Conformal Mapping:** An alternative approach, found in `Conformal2F1.jl`, which is likely used to improve convergence or accuracy in specific regions of the complex plane.
*   **Transformations:** The project utilizes various mathematical transformations of the $_2F_1$ function to handle different argument values effectively.

The project also includes extensive code for testing and validation, comparing the results of its methods against established libraries like `ArbNumerics` and even proprietary software like Mathematica. The presence of plotting libraries (`CairoMakie`, `GLMakie`) indicates that visualization of the results is also a key part of the research.

## Key Files

*   **`Project.toml` and `Manifest.toml`:** These files define the project's dependencies and structure, confirming it as a Julia project.
*   **`pFq.jl`:** This appears to be the main entry point, providing the user-facing functions `_2f1` and `_pfq` for computing hypergeometric functions.
*   **`Taylor.jl`:** Contains the implementation of the Taylor series method for solving the $_2F_1$ ODE.
*   **`Conformal2F1.jl`:** Implements the conformal mapping approach for computing $_2F_1$.
*   **`Initialization.jl`:** Handles the initial value setup for the ODE solver.
*   **`Transformations.jl`:** Contains the various transformation formulas for the $_2F_1$ function.
*   **`Notes.md`:** Provides valuable insights into the research, including TODOs, questions, and observations about the methods being developed.

## Building and Running

As a Julia project, the typical workflow would be to use the Julia REPL.

1.  **Open the Julia REPL.**
2.  **Navigate to the project directory:**
    ```julia
    cd("/home/merlin/Documents/Julia/Research/HypergeometricODEs")
    ```
3.  **Activate the project environment:**
    ```julia
    using Pkg
    Pkg.activate(".")
    ```
4.  **Instantiate the project (to install dependencies):**
    ```julia
    Pkg.instantiate()
    ```
5.  **Run the code:**
    To use the functions, you would include the main file and then call the desired functions. For example:
    ```julia
    include("pFq.jl")
    _2f1(1, 2, 3, 0.5)
    ```

## Development Conventions

*   The code is well-structured, with different functionalities separated into different files.
*   The use of `BenchmarkTools` suggests that performance is a consideration.
*   The `Notes.md` file indicates an active research and development process, with a clear set of goals and open questions.
*   The presence of a `Depreciated` folder suggests that the code is actively being refactored and improved.
