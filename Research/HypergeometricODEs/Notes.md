# Notes on the Taylor Method for ${_2}F_1$
> [!important] TODOs
> > [!caution] Top priority
> > - [ ] Fully implement conformal mapping approach
> > - [ ] Test conformal mapping
> > - [-] Eliminate poorly initialized transformations
> >     - z |-> 1/(1-z) is the culprit
> - [ ] Improve branch cut graphics.
>     - Separate sheet, reim plots.
>     - Use black and white wireframe plot
> - [ ] Add scaling factor for reliable radii. Do larger scaling factors fix these tests?

> [!note] Questions
> - Are there any common trends in the lower accuracy tests (not just failed tests)?
>     - Almost all of the lower accuracy tests have a dense ring of roots about z = 1 and/or z = 0.
>     - A few tests with c > 0 have roots but I am not sure if that is a full indication of the problem.
> - Are straight stepping paths more than enough?
> - Does the Wronskian indicate a preferred direction of travel?
>     - This may be more important for complex valued tests.
> - Can we separate the root and singular parts and stable parts of 2F1?
>     - Is the expansion about z = 1 into a regular and singular part the one we are looking for? No, the expansion is just a local form of the z |-> 1 - z transformation.
>     - What do the z |-> 1 - z 2F1 functions demonstrate the splitting?
>         - Not too much as far as the roots go. However, these graphics might be nice when looking into initialization errors.

> [!warning] Notes
> - The 1-2 digit loss in accuracy from some of the transformations is due to conditioning.
>     - For the tail error tests, there are regions of less than perfect conditioning.
> - Cécile's expansion about z = 1 is a local form of the z |-> 1 - z transformation formula.
> - The Wronskian IVP is slow an inaccurate when using a black box solver. Possibly more successful with our own IVP solver like the Taylor method.
> - The Wronskian integral is not the slowest but it doesn't improve accuracy at all in the "curious tests."
