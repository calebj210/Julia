# Notes on the Taylor Method for ${_2}F_1$
## TODOs
- [ ] Can we separate the root and singularity part and stable parts of 2F1?
    - [ ] Is the expansion about z = 1 into a regular and singular part the one we are looking for?
    - [x] What do the z |-> 1 - z 2F1 functions demonstrate the splitting?
        - Not too much as far as the roots go. However, these graphics might be nice when looking into initialization errors.
- [ ] Remove tail errors.
    - [ ] Add scaling factor for reliable radii. Do larger scaling factors fix these tests?
- [ ] Eliminate 1-2 digit initialization errors.
    - [ ] Check the transformation 2F1 functions individually.
- [ ] Improve branch cut graphics.
    - [ ] Separate sheet, reim plots.

## Questions
- Are there any common trends in the lower accuracy tests (not just failed tests)?
    - Almost all of the lower accuracy tests have a dense ring of roots about z = 1 and/or z = 0.
    - A few tests with c > 0 have roots but I am not sure if that is a full indication of the problem.
- Are straight stepping paths more than enough?
- Does the Wronskian indicate a preferred direction of travel?
    - This may be more important for complex valued tests.
- Can we remove the 1-2 digits of accuracy loss present in the transformed initializations?
    - This appears most _prominently_ for 1/z and 1/(1-z). Still they do appear for 1 - z and 1 - 1/z.
