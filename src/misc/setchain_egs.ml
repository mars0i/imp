(** Examples from Hartfiel's _Markov Set-Chains_, Springer 1998. *)


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

(************************************************************)
(** Example 2.15.
    See HartfielSetChainsErratat.txt. *)

(* Note: These are already tight. *)

let pmat = M.of_array [|0.423; 0.459; 0.043;
                        0.029; 0.677; 0.222;
                        0.0  ; 0.478; 0.461|]
                      3 3;;

let qmat = M.of_array [|0.473; 0.509; 0.093;
                        0.079; 0.724; 0.272;
                        0.036; 0.528; 0.511|]
                      3 3;;
