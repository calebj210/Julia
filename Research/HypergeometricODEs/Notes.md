# Notes on the Taylor Method for ${_2}F_1$
> [!important] TODOs
> > [!caution] Top priority
> > - [-] Finish documenting purely conformal mapping results.
> >     - Pass July 8th, July 31st, and new handout to Cécile.
> - [ ] Use Taylor method with conformal mapping.
> - [ ] Eliminate poorly initialized transformations
>     - z |-> 1/(1-z) is the culprit.
>     - Specifically, the closer initialization point was on the opposite side of the origin. 
> - [ ] Improve branch cut graphics.
>     - Separate sheet, reim plots.
>     - Use black and white wireframe plot

> [!note] Questions
> - Are there any common trends in the lower accuracy tests (not just failed tests)?
>     - Almost all of the lower accuracy tests have a dense ring of roots about z = 1 and/or z = 0.
>     - A few tests with c > 0 have roots but I am not sure if that is a full indication of the problem.
> - Are straight stepping paths more than enough?
> - Does the Wronskian indicate a preferred direction of travel?
>     - This may be more important for complex valued tests.
> - Can we separate the root and singular parts and stable parts of 2F1?
> - What is the ideal choice for selecting an initialization transformation?
>     - Does the smallest mapped z-value indicate the "best" choice?
>     - Does the center and radius of a transformation work even better?

> [!warning] Notes
> - The 1-2 digit loss in accuracy from some of the transformations is due to conditioning.
>     - For the tail error tests, there are regions of less than perfect conditioning.
> - Cécile's expansion about z = 1 is a local form of the z |-> 1 - z transformation formula.
> - The Wronskian IVP is slow an inaccurate when using a black box solver. Possibly more successful with our own IVP solver like the Taylor method.
> - The Wronskian integral is not the slowest but it doesn't improve accuracy at all in the "curious tests."
> - The conformal mapping approach offers partly better convergence compared the regular Maclaurin series. There are a few cases however where the mapping is not better. 
