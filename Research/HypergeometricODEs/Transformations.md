# Initialization Ideas
We can choose an initialization via:
* Get unbounded transformation radii.
* Check if transformed z is within any of the radii.
    - If z is in the one of the radii, initialize with the corresponding transformation.
* No direct evaluation found so continue to "abc" sign criteria to obtain a valid class of transformations.
    - If $\mathrm{Re}(z) > 1$, choose the transformation that puts $z_0$ to the right of $z = 1$.
        - Choose initialization spot so that straight stepping avoids $z = 1$.
    - Choose transformation that has the largest radius of convergence.

