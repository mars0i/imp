# imp
Miscellaneous experiments with imprecise probability in OCaml, with
particular interested in modeling evolutionary processes in which
probabilities are imprecise.

The following is not all self-explanatory.  Please feel free to write
to me for clarification.

*src/wrightfisher.ml* implements a diploid, single-locus Wright-Fisher
model with random mating and natural selection, using (discrete, not
continuous) sets of fitnesses, transition matrices, and initial
frequencies.  Plots of probability distributions at time t are easily
generated.  Perhaps it would be clearer to say that it implements
superimposed Wright-Fisher models.  Here is one simple inspiration for
the idea: Suppose that fitnesses differ depending on multiple 
environments, and one of these environments is experienced by the
population at (discrete) generation *t*.  Suppose further that which
environment is "chosen" for *t* is not *random* in the sense of there
being objective probabities for each environment, but rather
*erratic*, meaning that at the relevant level of description, there
are no relevant probabiities for which environment occurs.
Then this model displays the alternative probability distributions
that *may* characterize the frequencies in the population at each *t*.
(Please feel free to write to me for references or for draft papers or
presentations on this topic.)

The module `Wrightfisher` will end up in the probutils.cma package and
can be loaded into utop with `#load "probutils.cma"`.  This is compiled
by

    oasis setup
    make clean
    make

You will need to have installed opam and OCaml, and will need to use
opam to install the libraries Batteries, (Jane Street's) Core, and
Owl.

*src/wrightfisherPDFs.ml*: The above make invocation will also create
an executable that will run the Wright-Fisher simulation(s) and
generate PDF files.  (Restricted to a single intitial frequency.)  A
link to the executable named (probably) "wrightfisherPDFs.native" will
be generated as well.  Run it to see possible arguments and command
line options.

*src/setchains.ml* contains code specifically inspired by Hartiel's
*Markov Set-Chains*, Springer 1998.
