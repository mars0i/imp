

(* NOT equation (3) in Skulj ?? *)
let flowers max_n min_alg max_alg =
  let flower (idxs, lower) = 
    let complement_upper = (L.assoc (list_complement max_n idxs) max_alg)
    in (idxs, (max lower (1. -. complement_upper)))
  in L.map flower min_alg

(* NOT equation (4) in Skulj ?? *)
let fuppers max_n min_alg max_alg =
  let fupper (idxs, upper) = 
    let complement_lower = (L.assoc (list_complement max_n idxs) min_alg)
    in (idxs, (min upper (1. -. complement_lower)))
  in L.map fupper max_alg


let rec subtract_list xs ys =
  match xs with 
  | [] -> []
  | x::more_xs -> match ys with
    | [] -> xs
    | y::more_ys -> 
      if x = y 
      then subtract_list more_xs more_ys
      else x::(subtract_list more_xs ys)

(* Note on next_intsets: An alternative would be to skip the append 
 * and just have each element contain the new additions to the power 
 * set.  Then to get a power set, you have to concat all sets up until 
 * and including that one.  This is more space efficient but means 
 * re-appending again if you want another powerset. *)

(* Note on algebra_probs: Could use Sets instead of lists as keys, but they don't 
 * display their contents by default, which makes playing around at the repl 
 * inconvenient. *)

(* Perhaps it's a bit fragile, but since I'm representing sets of
 * atoms as lists of integers in decreasing order, we can use that
 * to compute complements.  Since this is for exploratory
 * experimentation, I prefer this representation to using Batteries or Core
 * Sets, which don't display their contents by default.
*)


(** OLD DEFINITIONS OF pri_f_field_lowers/uppers  **)

(** OLD Calculate a lower, L value from a probability interval.
    Equation (3) p. 1317 in Skulj 2009 (funny capitalization evokes the paper) *)
let pri_F_field_Lower  omega_max  atom_mins  atom_maxs  subset_idxs = 
  let mins_sum = prob_sum atom_mins subset_idxs in
  let maxs_comp_sum = prob_sum atom_maxs (list_complement omega_max subset_idxs) in
  max mins_sum (1. -. maxs_comp_sum)

(** OLD Calculate an upper, U value from a probability interval.
    Equation (4) p. 1317 in Skulj 2009 *)
let pri_F_field_Upper  omega_max  atom_mins  atom_maxs  subset_idxs = 
  let maxs_sum = prob_sum atom_maxs subset_idxs in
  let mins_comp_sum = prob_sum atom_mins (list_complement omega_max subset_idxs) in
  min maxs_sum (1. -. mins_comp_sum)

(** OLD Calculate L values for all members of the algebra and return an
    (atoms, L-value) alist *)
let pri_F_field_Lowers omega_max atom_mins atom_maxs =
  let algebra_idx_sets = LL.at algebra_sets omega_max in
  let lower idx_set =
    (idx_set, pri_F_field_Lower omega_max atom_mins atom_maxs idx_set) in
  L.map lower algebra_idx_sets

(** OLD Calculate U values for all members of the algebra and return an
    (atoms, U-value) alist *)
let pri_F_field_Uppers omega_max atom_mins atom_maxs =
  let algebra_idx_sets = LL.at algebra_sets omega_max in
  let upper idx_set = 
    (idx_set, pri_F_field_Upper omega_max atom_mins atom_maxs idx_set) in
  L.map upper algebra_idx_sets

(** Given *two* algebra_probs lists, return a similar alist in which
    values are the maximum of the two corresponding probabilities. *)
let algebra_maxs = algebra_extrema max


let compare m1 m2 =
  if m1 = m2 then 0 else
    if forall2 (<=) m1 m2 then -1 else 
      if forall2 (>=) m1 m2 then 1 
      else raise (Failure "incomparable matrices")
(* TODO This is inefficient.  It potentially loops through all pairs
 * of elements three times.  
 * A new strategy:
 * Look at this element.  If it's =, then keep looking for anything.
 * If this element is not zero, then if it's <, all subsequent pairs
 * must be < or =; else raise an exception.
 * If this element is >, all subsequent elements must be > or =; else
 * raise an exception.
 * If we are still in the = state when we get done, then return 0.
 * If we are in the < or > state, then return -1 or 1, respectively.
 * So maybe have two functions or branchs:
 * The equality function/branch, which is open-minded, but that will
 * choose the other function/branch if it encounters a non-equal pair.
 * To this other fn/branch is passed a relation, <= or >=, which is
 * uses repeatedly. *)

(* TODO? An alternative would be to return 1 rather than raising an
 * exception when the matrices are incomparable.  This is an odd choice,
 * conceptually, but since when 1, i.e. when m1 > m2, we have no interval--
 * i.e. the interval creation function will return None, we could just
 * return None in the incomparable case, too, by returning 1.  i.e. in
 * general 1 essentially means that you'll get a non-interval, so why not
 * just piggy-back on this behavior? *)
(* i.e. how about this function (which is still inefficient, though less so)? *)
let compare' m1 m2 =
  if m1 = m2 then 0 else
    if forall2 (<=) m1 m2 then -1 else 1
(* or even: *)
let compare'' m1 m2 =
  if forall2 (<=) m1 m2 then -1 else 1
(* That misrepresents when m1 = m1.  Is that a problem? *)


let make_3D_pdfs ?(altitude=45.) ?(azimuth=125.) basename start_gen last_gen distlists =
  let make_pdf i dists =  (* i = t-1; dists = prob dists at t *)
    let gen = i + start_gen in
    let filename = basename ^ (Printf.sprintf "%03d" gen) ^ ".pdf" in 
    let xs, ys, zs = make_coords (sort_dists dists) in
    let h = Pl.create filename in
      Pl.set_background_color h 255 255 255;
      Pl.set_foreground_color h 150 150 150; (* grid lines *)
      Pl.set_ylabel h "frequency of A allele";
      Pl.set_xlabel h "possible distributions";
      Pl.set_zlabel h "probability";
      Pl.set_altitude h altitude;
      Pl.set_azimuth h azimuth;
      Pl.mesh ~h xs ys zs;
      Pl.output h;
      Printf.printf "%s\n%!" filename
  in LL.iteri make_pdf (sub_lazy_list start_gen last_gen distlists)