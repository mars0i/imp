(** Functions for creating files from sequences of distributions of credal sets. *)

module Mat = Owl.Mat
module Pl = Owl.Plot
module OU = Owl.Utils
module L = Batteries.List
module A = Batteries.Array
module G = Utils.Genl
module T = Tdists

let (%) f g = (fun x -> f (g x))

(** Data file creation functions *)

let tdists_marshal_ext = "mltds"
(* let datafile_extension = ".mld" *)

(*
let tdists_list_to_idx_cols finite_tdists_list =
  let open T in
  let cols_for_one_tdists td =
    let len = L.length td.dists in
    Mat.concat_horizontal (Mat.create len 1 (float td.gen))
                          (Mat.of_array (A.of_list td.dists) len 1)
  in
  let llist_of_mats = T.map cols_for_one_tdists finite_tdists_list in
  (* concat 'em here *)
  ()
*)


(*
let write_csv_tdists_list basename finite_tdists_list =
  let open T in
  let first_gen = (T.hd finite_tdists_list).gen in
  let last_gen = (T.last finite_tdists_list).gen in
  let filename = Printf.sprintf "%s%02dto%02d.%s" basename first_gen last_gen "csv" in
  (* make header row *)
  (* construct data rows *)
  (* Mat.save_txt  ...*)
  ()
*)
  

(** [write_tdists_finite_list basename finite_tdists_list].
    Lazy list [finite_tdists_list] must be finite! *)
let marshal_tdists_list basename finite_tdists_list =
  let open T in
  let first_gen = (T.hd finite_tdists_list).gen in
  let last_gen = (T.last finite_tdists_list).gen in
  let filename = Printf.sprintf "%s%02dto%02d.%s" 
                                 basename first_gen last_gen tdists_marshal_ext in
  OU.marshal_to_file finite_tdists_list filename

(** Lazy list must be finite! *)
let marshal_tdists_sublist basename start_gen last_gen tdists_list =
  marshal_tdists_list basename
                           (T.sub_lazy_list start_gen last_gen tdists_list)

let unmarshal_tdists_list filename =
  ((OU.marshal_from_file filename) : T.tdistslist)

(** PDF plot-writing functions *)

let default_fontsize = 3.25
let default_plot_color = Pl.RGB (160, 160, 160)
let twoD_x_margin = 1.5
let twoD_y_bottom = ~-.0.012
let twoD_line_width = 1.
let interval_fill_color = default_plot_color
(* let interval_fill_color = Pl.RGB (180, 180, 180) *)
(* let interval_shrink_increment = 0.05 (* amount to shrink interval fill vertically so boundary lines will be visible *) *)

(** Given a list of probability distribution vectors, tries to sort them
    so that similar lists are close in the order. *) 
let l2_sort_dists dists = L.sort G.l2_compare dists
let simple_sort_dists dists = L.sort G.difference_compare dists
let abs_sort_dists dists = L.sort G.absdiff_compare dists

type pdfdims = TwoD | ThreeD | BothDs

let make_page_groups pdfdim plots_per_page tdistlists =
  let nonlazytdistlists = T.to_list tdistlists in
  let nonlazytdistlists' = 
    match pdfdim with (* if BothDs, we'll make 2 plots for each generation, so duplicate 'em *)
    | BothDs -> L.concat (L.fold_right (fun e acc -> [e;e]::acc) 
                                       nonlazytdistlists [])
    | _ -> nonlazytdistlists
  in
  let page_group_lists = L.ntake plots_per_page nonlazytdistlists'
  in L.map A.of_list page_group_lists  (* make sublists into arrays for easy indexing *)

(** Return a triple containing x-coord, y-coord, and z-coord matrices.
    dist_list is a list of row vectors representing prob dists that will
    be concatenated into z coords.  Note that the x and y coord matrices
    will have the same shape as each other, but with values increasing
    in different directions.  This shape will be transposed/rotated wrt the
    z coord matrix that results. That's what Owl.Plot.{mesh,surf} need.
    *)
let make_coords sample_interval dist_list =
  let (_, width) = Mat.shape (L.hd dist_list) in
  let height = L.length dist_list in
  let widthf, heightf, sample_intervalf = float width, float height, float sample_interval in
  let xs, ys = Mat.meshgrid 0. (heightf -. 1.)  0. ((widthf *. sample_intervalf) -. sample_intervalf)  height width in
  let zs = Mat.transpose (L.reduce Mat.concat_vertical dist_list) in
  (xs, ys, zs)
(* QUESTION: does using meshgrid obviate the need for set_ydigits below?
   Are there other advantages/disadvantages of meshgrid?  (It's harder to understand.) *)

(** Add a single 3D plot to handle h. To be used with make_3D_pdfs.  *)
let add_3D_plot ?plot_max ?fontsize ?colors ?addl_3D_fn h altitude azimuth xs ys zs =
  let open Pl in
  let size = match fontsize with Some x -> x | None -> default_fontsize in
  let plot_color = L.hd (match colors with Some x -> x | None -> [default_plot_color]) in (* only first color is used for 3D *)
  set_font_size h size;
  set_ylabel h "freq of A allele";
  set_xlabel h "poss distributions";
  set_zlabel h "probability";
  mesh ~h ~spec:[plot_color; NoMagColor; ZLine Y; 
                 Altitude altitude; Azimuth azimuth] xs ys zs;
  (match addl_3D_fn with | Some f -> f h xs ys zs | None -> ()); (* bare matches should always be wrapped *)
  (match plot_max with                                           (* even at end of fn so you don't forget if move it *)
  | Some z -> Pl.set_zrange h 0. z
  | None -> ())

(* Kludge to allow running with vanilla Owl that doesn't include it: *)
(* let set_ydigits h n = () *)
let set_ydigits h n = Plplot.plsyax n 0


(** Add a single 2D plot to handle h. To be used with make_pdfs.  *)
let add_2D_plot ?plot_max ?fontsize ?colors ?addl_2D_fn h ys zs =  (* Note ys are x-coordinates, zs are y-coordinates. *)
  let open Pl in
  let size = match fontsize with | Some x -> x | None -> default_fontsize in
  let plot_colors = match colors with Some x -> x | None -> [default_plot_color] in
  let num_plot_colors = L.length plot_colors in
  set_font_size h size;
  set_xlabel h "freq of A allele";
  set_ylabel h "probability";
  set_ydigits h 50;
  let m, n = Mat.shape ys in
  Pl.set_xrange h (~-. twoD_x_margin) ((float m) +. twoD_x_margin);
  (match plot_max with                                        (* even at end of fn so you don't forget if move it *)
  | Some y -> Pl.set_yrange h twoD_y_bottom y
  | None -> ()); (* I'd like to apply the lower margin here, too, but needs Plplot guts to hack default margins process *)
  (* Add addl plot stuff if a function was passed: *)
  (match addl_2D_fn with | Some f -> f h ys zs | None -> ()); (* bare matches should always be wrapped *)
  (* Add one curve to a single plot area for each column in data matrices: *)
  for i=0 to (n - 1) do 
    let plot_color = L.at plot_colors (i mod num_plot_colors) in
    plot ~h ~spec:[plot_color; LineWidth twoD_line_width] (Mat.col ys i) (Mat.col zs i);
  done


(** Function that can be passed to add_2D_plot via the [~addl_2D_fn] argument.
    Fills the region between the upper and lower curves with
    interval_fill_color. *)
[@@@ warning "-8"] (* disable match warning https://stackoverflow.com/a/46006016/1455243 *)
let fill_bounds ?(spec=[interval_fill_color; FillPattern 0]) h ys zs =  (* args are modeled on add_2D_plot *)
  let [|x1; x2|] = Mat.to_cols ys in
  let [|y1; y2|] = Mat.to_cols zs in
  let xs = Mat.(concat_vertical (reverse x1) x2) in
  let ys = Mat.(concat_vertical (reverse y1) y2) in
  let open Pl in
  draw_polygon ~h ~spec xs ys
[@@@ warning "+8"]


(** Make a series of n plot pdfs from tdistlists using basename.
    Example:
      let distlists = make_distlists 500 [200] 
                    [{w11=1.0; w12=0.8; w22=0.7}; {w11=1.0; w12=0.3; w22=0.7}];;
      make_pdfs "foo" 4 5 distlists;; 
    
    leftright=true Lay out plots from left to right before down, vs down first
    pdfdim=ThreeD  Make 3D plots, vs. TwoD for 2D or BothDs for both 3D and 2D
    rows=1         Number of rows of plots in PDF
    cols=1         Number of columns of plots in PDF
    altitude=20.   Viewing position: altitude parameter to mesh, plmesh
    azimuth=300.   Viewing position: azimuth parameter to mesh, plmesh
    sample_interval=1  Sample data every k frequencies rather than all frequencies 
    plot_max       Maximum height displayed (default: let Owl.Plot decide).
    fontsize       Font size.  Of course
    colors         2D: list of Plot.RGB's to cycle through (vs default_plot_color)
    addl_2D_fn     Additional function to run on the file handle, data for 2D plots
    addl_3D_fn     Additional function to run on the file handl, data  for 3D plots

    Note will throw an error if you try to make 3D plots with only one set of 
    input fitnesses. *)
let make_pdfs ?(leftright=true) ?(pdfdim=ThreeD) ?(rows=1) ?(cols=1) 
              ?(altitude=20.) ?(azimuth=300.) ?(sample_interval=1)
              ?plot_max ?fontsize ?colors ?addl_2D_fn ?addl_3D_fn
              basename tdistlists = 
  let plots_per_page = rows * cols in
  let max_row, max_col = rows - 1, cols - 1 in
  (* Next convert distlists--a possibly infinite lazylist of tdists--into a 
   * list of arrays of lists of vectors, where the elements of each rows*cols
   * length array are lists of vectors for a plot on rowsXcols sized page of plots. *)
  let page_groups = make_page_groups pdfdim plots_per_page tdistlists in
  (* fn to be applied to each array of lists of vectors to create a page: *)
  let make_pdf group_idx page_group = 
    (* Construct filename from basename and generation numbers: *)
    let group_len = A.length page_group in  (* differs if last group is short *)
    let generations_string = String.concat "_" (L.map (string_of_int % T.gen) (A.to_list page_group)) in
    let filename = Printf.sprintf "%s%s.pdf" basename generations_string in
    let first_of_two = ref true in (* allows staying with one generation for BothDs *)
    let h = Pl.create ~m:rows ~n:cols filename in
    Pl.set_background_color h 255 255 255; (* applies to all subplots *)
    let max_i, max_j = if leftright then max_row, max_col else max_col, max_row in (* order left right vs up down *)
    (* main loop through plots *)
    for i = 0 to max_i do
      for j = 0 to max_j do
        let row, col = if leftright then i, j else j, i in (* order left right vs up down *)
        Pl.subplot h row col;
        let idx = (i * (max_j + 1)) + j in 
        if idx < group_len then  (* don't index past end of a short group *)
          (Pl.set_foreground_color h 0 0 0; (* grid and plot title color *)
           let T.{gen; dists}  = page_group.(idx) in
           let xs, ys, zs = make_coords sample_interval (simple_sort_dists dists) in
           (* pre_title: either a newline (for 3D) plots or an empty string, so that
            *  titles on 3D plots will be pushed down a bit.*)
           let pre_title = match pdfdim, !first_of_two with
                           | BothDs, true  -> add_2D_plot ?plot_max ?fontsize ?colors ?addl_2D_fn h ys zs;
                                              first_of_two := false;
                                              ""
                           | BothDs, false -> add_3D_plot ?plot_max ?fontsize ?colors ?addl_3D_fn h altitude azimuth xs ys zs;
                                              first_of_two := true;
                                              "\n"
                           | TwoD, _       -> add_2D_plot ?plot_max ?fontsize ?colors ?addl_2D_fn h ys zs;
                                               ""
                           | ThreeD, _     -> add_3D_plot ?plot_max ?fontsize ?colors ?addl_3D_fn h altitude azimuth xs ys zs;
                                              "\n"
           in Pl.set_title h (pre_title ^ (Printf.sprintf "Generation %d" gen))
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

(** Calls make_pdfs with arguments suitable for displaying setchain bounds
    created using Hartfiel's hi-lo algorithm.  Set [~pdfdim] to [TwoD],
    [~addl_2D_fn] to default value [fill_bounds] to produce a gray fill,
    and [~colors] by default to a red and a blue.  Other arguments are as 
    for [make_pdfs] except for arguments that are specific to 3D plots. *)
let make_setchain_bounds_pdfs ?(addl_2D_fn=fill_bounds)
                              ?(colors=Pl.[RGB (0,0,200); RGB (200,0,0)])
                              ?(leftright=true) ?(rows=1) ?(cols=1)
			      ?(sample_interval=1) ?plot_max ?fontsize 
                              basename tdistlists =
  make_pdfs ~pdfdim:TwoD ~addl_2D_fn ~colors
            ~leftright ~rows ~cols ~sample_interval ?plot_max ?fontsize
            basename tdistlists

