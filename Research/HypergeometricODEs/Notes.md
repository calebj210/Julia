# Notes on the Taylor Method for ${_2}F_1$
## TODOs
- [ ] Remove tail errors.
    - [ ] Improve initialization criteria to include line-circle intersections.
    - [ ] Look at inner tail parameters.
    - [ ] Factor in z/(z-1) for Taylor initialization.
- [ ] Improve branch cut graphics.
    - [ ] Separate sheet, reim plots.

## Questions
- Are there any common trends in the lower accuracy tests (not just failed tests)?
- Are straight stepping paths more than enough?
- Does the Wronskian indicate a preferred direction of travel?
    - This may be more important for complex valued tests.
- Can we remove the 1-2 digits of accuracy loss present in the transformed initializations?
