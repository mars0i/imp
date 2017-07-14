(** wrightfisherPDFs:

    Executable wrapper that uses a Wright-Fisher diploid model with
    random mating to create PDF plots of probabilities of population-level 
    allele frequencies using "credal sets" of evolutionary transition 
    probabilities involving selection and drift but no other forces.  *)

module WF = Wrightfisher
module Command = Core.Command
module Spec = Core.Command.Spec

let sprintf = Printf.sprintf

let description = sprintf
"Creates plots of probabilities of population-level allele frequencies 
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
At least one triple is required (despite brackets above).

Example: %s foo 500 250 2 6  1.0 0.95 0.8  0.8 0.95 1.0"
(Filename.basename(Sys.executable_name))


let default_alt = 20
let default_az = 300

let alt_docstring = sprintf "integer aLtitude of perspective: degrees in [0,90] (default: %d)"  default_alt
let az_docstring =  sprintf "integer aZimuth of perspective: degrees in in [0,360] (default: %d)" default_az
let rows_docstring = sprintf "integer number of rows for multi-plot pages (default %d)" 1
let cols_docstring = sprintf "integer number of columns for multi-plot pages (default %d)" 1
let plot_max_docstring = "float If present, sets max height for all plots."
let fontsize_docstring = "float If present, sets font size."
let every_docstring = sprintf "integer plot only at every nth frequency (default %d)" 1
let twoD_docstring = "make 2D plots; with -3 make 2D and 3D side-by-side (default 3D)"
let threeD_docstring = "make 3D plots; with -2 make 2D and 3D side-by-side (default 3D)"
let updown_docstring = "If present arrange plots top bottom right; vs left right down."

let commandline =
  Command.basic
    ~summary:"Make 3D pdfs for multiple generations with multiple probability distributions."
    ~readme:(fun () -> description)
    Spec.(empty +> flag "-l" (optional_with_default default_alt int) ~doc:alt_docstring
                +> flag "-z" (optional_with_default default_az int)  ~doc:az_docstring
                +> flag "-r" (optional_with_default 1 int) ~doc:rows_docstring
                +> flag "-c" (optional_with_default 1 int) ~doc:cols_docstring
                +> flag "-m" (optional float) ~doc:plot_max_docstring
                +> flag "-f" (optional float) ~doc:fontsize_docstring
                +> flag "-e" (optional_with_default 1 int) ~doc:every_docstring
                +> flag "-2" no_arg ~doc:twoD_docstring
                +> flag "-3" no_arg ~doc:threeD_docstring
                +> flag "-u" no_arg ~doc:updown_docstring
                +> anon ("basename" %: string)
                +> anon ("popsize" %: int)
                +> anon ("initfreq" %: int)
                +> anon ("startgen" %: int)
                +> anon ("lastgen" %: int)
                +> anon (sequence ("fitn" %: float)))
    (fun alt_int az_int rows cols plot_max fontsize every twoD threeD updown basename popsize initfreq startgen lastgen fitn_floats () ->
      let altitude = float alt_int in
      let azimuth = float az_int in
      let fitn_recs = WF.group_fitns fitn_floats in
      let distlists = WF.make_distlists popsize [initfreq] fitn_recs in
      let pdfdim = match twoD, threeD with
                   | true, true   -> WF.BothDs
                   | true, false  -> WF.TwoD
                   | false, true | false, false -> WF.ThreeD (* default *)
      in WF.make_pdfs ~rows ~cols ~altitude ~azimuth ~every ~pdfdim ?plot_max ?fontsize ~leftright:(not updown) basename startgen lastgen distlists)

let () = Command.run ~version:"1.1" ~build_info:"wrightfisherPDFS, (c) 2017 Marshall Abrams" commandline
