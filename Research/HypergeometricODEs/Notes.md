# Notes on the Taylor Method for ${_2}F_1$
## TODOs
- [ ] Remove tail errors.
- [ ] Improve branch cut graphics.
    - [ ] Pin down branch cuts.
    - [ ] Seamlessly stitch sheets.

## Questions
- Are there any common trends in the lower accuracy tests (not just failed tests)?
- Are straight stepping paths more than enough?
- Does the Wronskian indicate a preferred direction of travel?
- Can we remove the 1-2 digits of accuracy loss present in the transformed initializations?

## Draft Comments
- [x] 6:  “… for it’s many parameters” reminds me of a story about the famous English mathematician G.H. Hardy, who started a lecture by “Let us consider a large integer, say 2”. There are precisely 4 parameters, so “many” is a bit misleading. We could say something like “… for any value of its argument and its three additional parameters”.
- [x] 21: I would hesitate to say “the greatest”, maybe simpler just “a”.
- [x] 24: “method” -> “study”.
- [x] 26: “and z” -> “and argument z” .
- [x] 35: “where there are” -> “which has”.
- [x] 38: I suggest removing “invite you to”.
- [x] 45: “coursed into” – you probably mean “coerced into”
- [x] 48: “many” unnecessary.
- [x] 51: “immense” -> “extensive”
- [x] 77: “computation has required a slew of techniques” -> “computations have been implemented in many different ways,”
- [x] 79: “use intricate” -> “use an intricate”. In the following Figure 1, the legends are much too small (exponents on the vertical axes for Figures 4, 5, 6 are also on the small side).
- [x] 117: “often fails to provide accuracy nor speed” -> “is computationally inefficient” (we can’t say it fails to provide accuracy, as it will do so under uneconomical refinement).
- [x] 118: “A nice and simple improvement of Euler’s method is the so called” -> “A significant improvement is provided by the”
- [x] 143: “it is more accurate to step towards the origin rather than away in floating point arithmetic” -> “it is in floating point arithmetic more accurate to step towards the origin rather than away”
- [x] 147: “gives us the” -> “gives the”
- [x] 163: “Johansson” -> “the one used in [18]” (not essential, but it might be nicer to avoid using names when criticizing).
- [x] 169: “requires” -> “often requires”
- [x] 191-192: Table 1 caption: “Ideal transformations to use” -> “Transformations used” (words like “ideal” and “optimal” have been heavily misused by many authors, so some readers might be a bit put off seeing them). It is very clear anyway from the caption that what is listed are our preferences.
- [x] 202: “section” -> “Example”.
- [x] 259: “another” -> “the second”.
- [x] 261: “an alternate” ->  “the second”.
- [x] 263: “Strikingly, the alternate branch picks up a singularity at the origin.” -> “Additional branches will in general feature a singularity also at the origin (the form of the expansion (1.1) prevents this from appearing on the primary sheet).”
- [x] 281: You can after this line add a footnote to the effect that in these tests, the result each z-value was calculated separately, i.e., we have not utilized that calculations across a grid can for the Taylor methos be vastly speeded up by following path strategies as used in the above cited Painlevé calculations.
- [x] 286: “We wish” -> “We now wish”
- [x] 288: What does “uniformly sampled” mean? Grid-based, or (presumably) random with statistically uniform sampling density?
- [x] 307: “2f1” ->  “2F1”
- [x] 323: “computation” -> “Computation”
- [x] 327: “Mccoy” -> “McCoy”
- [x] 356” “love” ->  “Love” and 357 “schwarzschild” -> “Schwarzschild”
- [x] 365: “odes” -> “ODEs” and 366 ”taylor” -> “Taylor”
- [x] 376: “”auxiliary” -> “Auxiliary”

- [ ] 260: Add a second figure that makes more clear how the two sheets join (we can safely assume that lots of readers are quite unfamiliar with the visual aspects of multiple sheets).
- [ ] 296: As we have discussed, we need to get totally rid of this tail. Scrutinize the worst cases to see just what goes wrong, and either update the general method or put in some test to identify and then treat separately these cases.
- [ ] 301: We probably want to replace “4.3. Remarks” with “5. Conclusions”. This Chapter needs then to be expanded a bit (and, in particular, not end mid-sentence). Do you think that we can also say something along the lines that the key ideas of the Taylor approach can be expected to generalize to the confluent 1F1 case and also to p+1Fp cases (p = 2, 3, …), however then with the complication that fewer functional identities are then available.  
