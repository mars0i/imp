(** setchainPDFs:

    Executable wrapper that uses Wright-Fisher diploid models with
    random mating, selection, and drift to create PDF plots of setchains
    of population-level allele frequencies. *)

module Command = Core.Command
module Spec = Core.Command.Spec
module SC = Models.Setchains
module WF = Models.Wrightfisher
module IO = Models.CredalsetIO
module T = Models.Tdists
module G = Utils.Genl
module Pl = Owl.Plot

(* TODO: 
 * Fix docstring(s).
 * Add option to control fill color.
 * Add option to control number of pmap forks.
 *)

let bottom_top_colors=Pl.[RGB (0,0,200); RGB (200,0,0)]

let sprintf = Printf.sprintf

let description = sprintf
"FIXME Creates plots of probabilities of population-level allele frequencies 
using \"credal sets\" of evolutionary transition probabilities, which are
generated by diploid Wright-Fisher models with random mating and selection.
Each plot will have its own PDF file.

BASENAME: PDF filenames consist of this string followed by generation number.

POPSIZE, INITFREQ: Integers specifying the 2N number of alleles in the
population and the initial frequency of the A allele.

STARTGEN LASTGEN: First and last generations for which to generate PDFs.
(Generation 0 in which INITFREQ has a  probability of 1 can't be plotted.)

[FITN ...]: Triples of fitness values in the order:
    w11 (AA homozygote), w12 (heterozygote), w22 (BB heterozygote)
At least one triple is required (despite brackets above).  The setchain
will be constructed from the the tight matrix interval containing the
the transition matrices for these fitness triples.

Example (FIXME): %s foo 500 250 2 6  1.0 0.95 0.8  0.8 0.95 1.0"
(Filename.basename(Sys.executable_name))

let rows_docstring = sprintf "integer number of rows for multi-plot pages (default %d)" 1
let cols_docstring = sprintf "integer number of columns for multi-plot pages (default %d)" 1
let plot_max_docstring = "float If present, sets max height for all plots."
let fontsize_docstring = "float If present, sets font size."
let sample_docstring = sprintf "sample data to plot only at every nth frequency (default %d)" 1
let skip_docstring = sprintf "skip to every nth generation (default %d)" 1
let updown_docstring = "If present arrange plots top bottom right; vs left right down."
let nofork_docstring = "If present, don't split computation across multiple cores/CPUs."

let commandline =
  Command.basic
    ~summary:"Make 3D pdfs for multiple generations with setchains."
    ~readme:(fun () -> description)
    Spec.(empty +> flag "-r" (optional_with_default 1 int) ~doc:rows_docstring
                +> flag "-c" (optional_with_default 1 int) ~doc:cols_docstring
                +> flag "-m" (optional float) ~doc:plot_max_docstring
                +> flag "-f" (optional float) ~doc:fontsize_docstring
                +> flag "-s" (optional_with_default 1 int) ~doc:sample_docstring
                +> flag "-k" (optional_with_default 1 int) ~doc:skip_docstring
                +> flag "-u" no_arg ~doc:updown_docstring
                +> flag "-1" no_arg ~doc:nofork_docstring
                +> anon ("basename" %: string)
                +> anon ("popsize" %: int)
                +> anon ("initfreq" %: int)
                +> anon ("startgen" %: int)
                +> anon ("lastgen" %: int)
                +> anon (sequence ("fitn" %: float)))
    (fun rows cols plot_max fontsize sample skip updown nofork basename popsize initfreq startgen lastgen fitn_floats () ->
      let fitn_recs = WF.group_fitns fitn_floats in
      Printf.printf "making matrix interval ... %!";
      let pmat, qmat = SC.make_wf_interval popsize fitn_recs in
      Printf.printf "making lazy bounds mats list ... %!";
      let bounds_mats =  SC.lazy_bounds_mats_list ~fork:(not nofork) pmat qmat in
      Printf.printf "making lazy prob intervals list ... %!";
      let selected_gens = G.lazy_ints ~skip:skip 1 in
      let tdistlists = T.add_gens (SC.lazy_prob_intervals_from_freq initfreq bounds_mats) in
      let selected_distlists = T.sublist startgen lastgen (T.select_by_gens selected_gens tdistlists) in
      Printf.printf "making pdfs ... \n%!";
      IO.make_setchain_bounds_pdfs ~colors:bottom_top_colors
                    ~rows ~cols ~sample_interval:sample ?plot_max ?fontsize ~leftright:(not updown)
                    basename selected_distlists) (* startgen lastgen FIXME *)

let () = Command.run ~version:"1.0" ~build_info:"setchainPDFS, (c) 2017 Marshall Abrams" commandline
