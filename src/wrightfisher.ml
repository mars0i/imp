
(** Wrightfisher:
    Calculate Wright-Fisher transition probabilities with selection and
    create plots of probabilities of population frequencies.

    Based on Ewens _Mathematical Population Genetics I, 2nd ed, 
    equations 1.58, 1.59, and 1.25, though similar formulas can be 
    found in many places. *)

module Mat = Owl.Mat
module Math = Owl.Maths (* note British->US translation *)
module Pl = Owl.Plot
module L = Batteries.List
module A = Batteries.Array
module LL = Batteries.LazyList

module U = Utils

let ( *@ ) = Mat.( *@ )  (* = dot: matrix multiplication *)

let is_odd n = n mod 2 <> 0
let is_even n = n mod 2 = 0


(** Return a lazy list that's a sublist of the argument, from element start 
    (zero-based) to element finish, inclusive. *)
let sub_lazy_list start finish ll =
  LL.take (finish - start + 1) (LL.drop start ll)

let make_init_dist allele_popsize a1count =
  let m = Mat.zeros 1 (allele_popsize + 1) in
  Mat.set m 0 a1count 1.0;
  m

type fitnesses = {w11 : float; w12 : float; w22 : float}

(** 1.59 in Ewens *)
let weight_i {w11; w12; w22} allele_popsize freq  =
  let i, i' = float freq, float (allele_popsize - freq) in
  let a_hom = w11 *. i *. i in
  let het   = w12 *. i *. i' in
  let b_hom = w22 *. i' *. i' in
  (a_hom +. het) /. (a_hom +. 2. *. het +. b_hom)

(** Wright-Fisher transition probability from frequency prev_freq (row index)
    to frequency next_freq (column index). *)
let prob_ij fitns allele_popsize prev_freq next_freq =
  let wt = weight_i fitns allele_popsize prev_freq in
  let other_wt = 1. -. wt in
  let j = float next_freq in
  let j' = float (allele_popsize - next_freq) in
  let comb = Math.combination_float allele_popsize next_freq in
  comb  *.  wt**j  *.  other_wt**j'

(** Make a transition matrix from fitnesses *)
let make_tranmat allele_popsize fitns =
  (* prob_ij with an extra ignored argument: *)
  let prob_ij_ fitns allele_popsize prev_freq next_freq _ =
    prob_ij fitns allele_popsize prev_freq next_freq
  in
  let dim = allele_popsize + 1 in
  let m = Mat.empty dim dim  in
  Mat.mapi (prob_ij_ fitns allele_popsize) m

(** Return pair of old distribution vector (i.e. the second argument) and
    distribution vector (i.e. the product of the two arguments), for use 
    as an element in a Batteries.LazyList *)
let next_dist tranmat dist = 
  (dist, dist *@ tranmat)

(** Returns a LazyList of probability distributions that are the result of 
    a stationary Markov chain with initial distribution (usually with 1 at 
    one entry and zero elsewhere, representing that the population has that 
    initial frequency) and a transition matrix. The initial distribution is 
    included in the list.*)
let make_dists tranmat init_dist =
  LL.from_loop init_dist (next_dist tranmat)

(** Convenience function: Takes elements from start to finish, inclusive, 
    from a LazyList, and convert the result to a List. *)
let take_to_list start finish ll = 
  LL.to_list (sub_lazy_list start finish ll)

(** Return second dimension of a matrix or vector. *)
let length m = snd (Mat.shape m)

(* Make a series of n 2D plot pdfs from dists using basename. *)
let make_2D_pdfs basename dists n =
  let dist_length = length (LL.at dists 0) in
  let xs = Mat.sequential 1 dist_length in (* vector of x-axis indices *)
  let make_pdf i dist =
    let filename = basename ^ (Printf.sprintf "%03d" i) ^ ".pdf" in
    let h = Pl.create filename in 
    Pl.set_yrange h 0.0 0.25;
    Pl.scatter ~h xs dist; 
    Pl.output h
  in LL.iteri make_pdf (LL.take n dists)

(** Return a triple containing x-coord, y-coord, and z-coord matrices.
    dist_list is a list of row vectors representing prob dists that will
    be concatenated into z coords.  Note that the x and y coord matrices
    will have the same shape as each other, but with values increasing
    in different directions.  This shape will be transposed/rotated wrt the
    z coord matrix that results. That's what Owl.Plot.{mesh,surf} need.
    For Wright-Fisher pdfs, I use y as frequency; x indexes probability
    distributions.
    *)
let make_coords ?(every=1) dist_list =
  let dist_list' = L.map (U.subsample_in_rows every) dist_list in (* identical if every=1 *)
  let (_, width) = Mat.shape (L.hd dist_list') in
  let height = L.length dist_list' in
  let everyf = float every in
  let xs = Mat.repeat ~axis:0 (Mat.sequential 1 height) width in
  let ys = Mat.repeat ~axis:1 (Mat.sequential ~step:everyf width 1) height in (* step so freqs match z vals if every>1 *)
  let zs = L.reduce Mat.concat_vertical dist_list' in
  (xs, ys, zs)


(** Given a list of transition matrices, mats, and a list of
    probability distributions, dists, returns a new list of
    probability distributions produced by multiplying all
    distributions by all matrices. *)
let next_dists tranmats dists =
  (dists, L.concat (L.map (fun dist -> L.map (Mat.dot dist) tranmats)
                          dists))

(** Given a list of transition matrices and a list of initial distributions
    (often one distribution with all weight on one frequency),
    return a lazy list of lists of probability distributions.
    Note this function does not not drop the first element. That way, the 
    number of dists in the nth distlist = (length tranmats)**n for one 
    initial distribution.  e.g. with
    two transition matrices and one initial distribution,
       LL.at distlists 1
    will produce 2**1 = 2 dists.  Or with more initial distributions, the
    number of dists at n is (length init_dists) * (length tranmats)**n . *)
let make_distlists_from_mats tranmats init_dists =
  LL.from_loop init_dists (next_dists tranmats)

(** Like make_distlists_from_mats, but uses basic parameters to generate the transition
    matrices and initial distributions that are arguments to make_distslists.
    Arg 1, size is the number of alleles, i.e. 2N; arg 2, init_freqs, is a 
    list of all initial frequencies for the population (usually there is only
    one, so the list will have only one element); arg 3, fitn_list is a list
    of fitness structures. *)
let make_distlists size init_freqs fitn_list =
  let init_dists = L.map (make_init_dist size) init_freqs in
  let tranmats = L.map (make_tranmat size) fitn_list in
  make_distlists_from_mats tranmats init_dists

(* Given a list of probability distribution vectors, tries to sort them
   so that similar lists are close in the order. *) 
let l2_sort_dists dists = L.sort U.l2_compare dists
let simple_sort_dists dists = L.sort U.difference_compare dists
let abs_sort_dists dists = L.sort U.absdiff_compare dists

(** Given a list of float fitness values, which should be in the order
       w11, w12, w22, w11, w12, w22, ...
    eat them in groups of three, using each three to create a
    fitness record and return a list of these records in order. 
    This can be used for commandline processing. *)
let group_fitns fitn_float_list =
  let rec loop l acc =
    match l with
    | [] -> acc
    | w11::w12::w22::tl -> loop tl ({w11=w11; w12=w12; w22=w22}::acc)
    | _ -> raise (Failure "Missing/extra fitness(es)")
  in L.rev (loop fitn_float_list [])

let plot_color = (120, 40, 0)

(** Add a single 3D plot to handle h. To be used with make_3D_pdfs.  *)
let add_3D_plot h altitude azimuth xs ys zs =
  Pl.set_altitude h altitude;
  Pl.set_azimuth h azimuth;
  Pl.set_ylabel h "freq of A allele";
  Pl.set_xlabel h "poss distributions";
  Pl.set_zlabel h "probability";
  Pl.plots2d3d ~h ~color:plot_color xs ys zs

(** Add a single 2D plot to handle h. To be used with make_pdfs.  *)
let add_2D_plot h ys zs =
  Pl.set_xlabel h "freq of A allele";
  Pl.set_ylabel h "probability";
  let _, n = Mat.shape ys in
  for i=0 to (n - 1) do 
    Pl.plot ~h ~color:plot_color (Mat.col ys i) (Mat.row zs i)
  done

type pdfdims = TwoD | ThreeD | BothDs

let make_page_groups pdfdim plots_per_page finite_lazy_distlists =
  let distlists = LL.to_list finite_lazy_distlists in
  let distlists' = 
    match pdfdim with (* if BothDs, we'll make 2 plots for each generation, so duplicate 'em *)
    | BothDs -> L.concat (L.fold_right (fun e acc -> [e;e]::acc) 
                                       distlists [])
    | _ -> distlists
  in
  let page_group_lists = L.ntake plots_per_page distlists'
  in L.map A.of_list page_group_lists  (* make sublists into arrays for easy indexing *)


(** Make a series of n 3D plot pdfs from distlists using basename.
    Example:
    let distlists = make_distlists 500 [200] 
                  [{w11=1.0; w12=0.8; w22=0.7}; {w11=1.0; w12=0.3; w22=0.7}];;
    make_pdfs "distsN=500init=200w11=1w22=0.7w12=0.8or0.3gen" distlists 9;; *)
(* Turned into spaghetti when I tried to add option of two different plots.  needs redoing. *)
let make_pdfs ?(pdfdim=ThreeD) ?(rows=1) ?(cols=1) ?(altitude=30.) ?(azimuth=300.) ?(every=1)
                         basename start_gen last_gen distlists =
  let plots_per_page = rows * cols in
  let max_row, max_col = rows - 1, cols - 1 in
  let gens_per_page = match pdfdim with  (* BothDs means two different plots per generation *)
                      | BothDs -> if is_odd plots_per_page 
                                  then raise (Invalid_argument "Two plots per generation (pdfdim=BothDs) with odd plots per page.")
                                  else plots_per_page / 2        (* Rarely makes sense and I don't want to handle the case: *)
                      | _ -> plots_per_page
  in
  (* Next convert distlists--an infinite lazylist of lists of vectors--into a 
   * list of arrays of lists of vectors, where the elements of each rows*cols
   * length array are lists of vectors for a plot on rowsXcols sized page of plots. *)
  let finite_lazy_distlists = sub_lazy_list start_gen last_gen distlists in     (* lazy list of only those gens we want *)
  let page_groups = make_page_groups pdfdim plots_per_page finite_lazy_distlists in
  (* fn to be applied to each array of lists of vectors to create a page: *)
  let make_pdf group_idx page_group = 
    (* Construct filename from basename and generation numbers: *)
    let group_len = A.length page_group in  (* differs if last group is short *)
    let group_start = start_gen + group_idx * gens_per_page in
    let group_last  = group_start + gens_per_page - 1 in
    let filename = basename ^ 
                   (Printf.sprintf "%02dto%02d" group_start group_last ) ^
                   ".pdf"
    in
    let h = Pl.create ~m:rows ~n:cols filename in
    Pl.set_background_color h 255 255 255; (* applies to all subplots *)
    let first_of_two = ref true in
    for row = 0 to max_row do
      for col = 0 to max_col do
        let idx = (row * cols) + col in  (* row not rows *)
        Pl.subplot h row col;
        if idx < group_len then  (* don't index past end of a short group *)
          (Pl.set_foreground_color h 170 170 170; (* grid color *)
           let xs, ys, zs = make_coords ~every (l2_sort_dists page_group.(idx))
           in
           match pdfdim, !first_of_two with
           | BothDs, true -> add_3D_plot h altitude azimuth xs ys zs;
                             first_of_two := false
           | BothDs, false -> add_2D_plot h ys zs;
                              first_of_two := true;
           | TwoD, _ -> add_2D_plot h ys zs;
           | ThreeD, _ -> add_3D_plot h altitude azimuth xs ys zs;
         )
        else (* short group *)
          (* Dummy plot to prevent plplot from leaving a spurious border: *)
          (Pl.set_foreground_color h 255 255 255; (* s/b same color as bg *)
           Pl.plot_fun ~h (fun x -> x) 0. 1.;)     (* y values must change *)
      done
    done;
    Pl.output h;
    Printf.printf "%s\n%!" filename
  in L.iteri make_pdf page_groups

(* Use this if Owl starts adding titles to 3D plots (write now that doesn't work):
           let gen = start_gen + group_idx + idx in
           Pl.set_title h (Printf.sprintf "generation %d" gen);
 *) 
