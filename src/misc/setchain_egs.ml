(** Examples from Hartfiel's _Markov Set-Chains_, Springer 1998.
    For use with module Models.Setchains. *)


(************************************************************)
(** Example 2.10 *)

let mp  = M.of_array [|0.35; 0.55; 0.25; 0.65|] 2 2;;
let mq  = M.of_array [|0.45; 0.65; 0.35; 0.75|] 2 2;;
let s0p = M.of_array [|0.4; 0.5|] 1 2;;
let s0q = M.of_array [|0.5; 0.6|] 1 2;;
let m = [mp; mq];;
let s0 = [s0p; s0q];;

let m_verts = mat_vertices ~uniq:true ~digits:3 mp mq;;
let s0_verts = mat_vertices ~uniq:true ~digits:3 s0p s0q;;

let w1 = Pm.cross_apply M.dot s0_verts m_verts;;  (* dot = ( *@ )  *)
let s1_verts = [(L.last w1); (L.first w1)];;  (* (0.29, 0.71), (0.4, 0.6); see Hartfiel *)

let w2 = Pm.cross_apply M.dot s1_verts m_verts;;

(************************************************************)
(** Example 2.13 *)

let p = M.of_array [|0.25; 0.0; 0.25|] 1 3;;
let q = M.of_array [|0.75; 0.0; 0.75|] 1 3;;

let l = M.of_array [|0.25; 0.25; 0.5|] 3 1;;
let h = M.of_array [|0.75; 0.75; 0.5|] 3 1;;

(************************************************************)
(** Example 2.14 (same as Ex. 2.11).
    See HartfielSetChainsErratat.txt. *)

(* Note: These are already tight. *)

let pmat = M.of_array [|0.0;  0.25; 0.25;
                        0.25; 0.50; 0.25;
                        0.25; 0.25; 0.0|]
                      3 3;;

let qmat = M.of_array [|0.0;  0.75; 0.75;
                        0.25; 0.50; 0.25;
                        0.75; 0.75; 0.0|]
                      3 3;;

(** Modifies a matrix to introduce Hartfiel's typo, for easier 
    comparison with the text.  See doc/HartfielMarkovSetChainsErrata.txt *)
let swap12 m = 
  let m' = M.copy m in
  M.swap_rows m' 1 2;
  m';;

let swap12both (lo, hi) = swap12 lo, swap12 hi;;

(************************************************************)
(** Example 2.15.
    See HartfielSetChainsErratat.txt for notes on how my results do
    and don't match what's in the book. *)

(* Note: These are already tight. *)

let pmat = M.of_array [|0.423; 0.459; 0.043;
                        0.029; 0.677; 0.222;
                        0.0  ; 0.478; 0.461|]
                      3 3;;

let qmat = M.of_array [|0.473; 0.509; 0.093;
                        0.079; 0.724; 0.272;
                        0.036; 0.528; 0.511|]
                      3 3;;

let p = M.of_array [|0.4; 0.1; 0.2|] 1 3;;
let q = M.of_array [|0.6; 0.3; 0.4|] 1 3;;

(************************************************************)
(* example using make_wf_interval: *)

let n = 100;;
let p', q' = W.(make_wf_interval n {w11=1.0; w12=0.3; w22=0.1} {w11=1.0; w12=0.9; w22=0.5});;
let p, q = tighten_mat_interval p' q';;


(************************************************************)
(* ad hoc example *)

let p, q = 
  let size = 6 in 
  let x = 1. /. (float size) in
  let p' = M.(x $- ((uniform 1 size) *$ 0.1)) in
  let q' = M.(p' + ((uniform 1 size) *$ 0.2)) in
  tighten_vec_interval p' q';;
