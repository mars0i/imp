(** wrightfisherPDFs:

    Executable wrapper that uses a Wright-Fisher diploid model with
    random mating to create PDF plots of probabilities of population-level 
    allele frequencies using "credal sets" of evolutionary transition 
    probabilities involving selection and drift but no other forces.  *)

module WF = Wrightfisher
module Command = Core.Command
module Spec = Core.Command.Spec

let description = 
"Uses a diploid Wright-Fisher model with random mating to create PDF
plots of probabilities of population-level  allele frequencies using
\"credal sets\" of evolutionary transition probabilities involving
selection and drift but no other forces.

BASENAME: a string with which all PDF names will begin.

POP_SIZE, INITIAL_FREQ: Integers specifying the 2N number of alleles in the
population and the initial frequency of the A allele.

START_GENERATION, LAST_GENERATION: The first and last generations for
which  to generate PDFs.  Generation 0, in which INITIAL_FREQ in effect
has a  probability of 1, cannot be plotted, so the first possible value
for START_GENERATION is 1.  Note that there will be 
<number of fitness triples>^t probability distributions over frequencies
at generation t, so each generation will take longer to process than the 
preceding one, and it will be impractical to generate PDFs for many 
generations beyond 10.

[FITNESS ...]: Triples of fitness values in the order:
    w11 (AA homozygote), w12 (heterozygote), w22 (BB heterozygote)
At least one triple is required (despite brackets above).  Each triple
will define a transition matrix, which will be applied to each probability 
distribution over frequencies at one generation in order to generate the 
probability distributions at the next generation.

Example:
<executablename> \"twodirs\"  1000 500   2 8   1.0 0.95 0.8   0.8 0.95 1.0
"

let commandline =
  Command.basic
    ~summary:"Make 3D pdfs for multiple generations with multiple probability distributions."
    ~readme:(fun () -> description)
    Spec.(empty +> anon ("basename" %: string)
                +> anon ("pop_size" %: int)
                +> anon ("initial_freq" %: int)
                +> anon ("start_generation" %: int)
                +> anon ("last_generation" %: int)
                +> anon (sequence ("fitness" %: float)))
    (fun basename pop_size init_freq start_gen last_gen fitn_floats () ->
      let fitn_recs = WF.group_fitns fitn_floats in
      let distlists = WF.make_distlists pop_size [init_freq] fitn_recs in
      WF.make_3D_pdfs basename distlists start_gen last_gen)

let () = Command.run ~version:"1.0" ~build_info:"wrightfisherPDFS, (c) 2017 Marshall Abrams" commandline
