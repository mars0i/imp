It's clear from looking at owl_plot.ml that `create` doesn't cause an
Plplot functions to be called, so it can't possibly open a file through
Plplot.  The only thing that `create` does is to store information in
data structures.  Some of this information consists of functions, but
they're not executed afaics.  (`create` doesn't explicitly open a file,
either.)  It must be `output` that opens the actual file for writing,
probably in `plinit`, which is called in `_initialize`--one of the
internal functions that `output` calls.

Given how my code is written, this means that when I create new
processes using `Parmap`, the PDF file is not open.  The code run via
`Parmap` just helps to create data that is then incorporated into thunks
and data structures that get used when `Owl.Plot.output` is called
later.  So I don't have to worry about duplicated file handles.  I don't
know whether this sort of situation was one of the reasons for the
decision to delay opening the output file in `Owl.Plot`, but it's a good
design.
