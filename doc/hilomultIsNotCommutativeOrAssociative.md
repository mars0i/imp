hilo multiplication is neither associative nor commutative
===


hilo multiplication is not commutative:

```ocaml
let pmat, qmat = S.make_wf_interval 200 [W.{w11=1.0; w12=0.3; w22=0.7}; W.{w11=1.0; w12=0.5; w22=0.2}];;
let hl2 = S.next_bounds_mats pmat qmat psum qsum (pmat, qmat);;
let hl3 = S.next_bounds_mats pmat qmat psum qsum hl2;;
let hl3' = let p', q' = hl2 in S.next_bounds_mats p' q' psum qsum (pmat, qmat);;
let lodif, hidif = M.(sum ((fst hl3) - (fst hl3')), sum ((snd hl3) - (snd hl3')));;

val lodif : M.mat =
       C0
R 5.52946

val hidif : M.mat =
        C0
R -120.498
```

Since `hl2` was just the hilo product of `(pmat, qmat)` with itself,
this also shows that hilo multiplication is not associative.
