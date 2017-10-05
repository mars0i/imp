(* 
   Example illustrating difference between setchain interval results
   where we start from three fitness specifications and use the
   resulthing three transition matrices to construct the minimal
   interval containing them all, vs. simply iterating those three tran
   matrices to get 3^t new dists at every time t.

   See
     SetchainsIntervalsVsIndivDists_3dists01to06.pdf
     SetchainsIntervalsVsIndivDists_setchain01to06.pdf

   Things to notice:

     At time 1, you can see exactly how the initial three dists map give
     you the setchain interval.

     At subsequent times, the setchain interval goes higher in many
     places.  That's because, presumably, the interval contains many
     more transition operators, so the opportunities for combining
     values to produce a wider interval are greater.
*)

#use "src/loadit.ml"

let fitnesses = [W.{w11=1.0; w12=0.7; w22=0.5}; W.{w11=0.5; w12=1.0; w22=0.5}; W.{w11=0.2; w12=0.8; w22=1.0}];;

let pmat, qmat = S.make_wf_interval 100 fitnesses;;
let distlists = S.lazy_prob_intervals_from_freq 50 (S.lazy_bounds_mats_list pmat qmat);;
Models.CredalsetPDF.(make_setchain_bounds_pdfs "SetchainsIntervalsVsIndivDists_setchain" ~rows:2 ~cols:3 ~plot_max:0.12 1 6 distlists);;

let distlists3 = W.make_distlists 100 [50] fitnesses;;
Models.CredalsetPDF.(make_pdfs "SetchainsIntervalsVsIndivDists_setchain_3dists" ~rows:2 ~cols:3 ~pdfdim:TwoD ~colors:[RGB (255,0,0)] ~plot_max:0.12 1 6 distlists3);;

(******************************)

(* loadit.ml contains: *)

(*
#require "parmap";;
#load "utils.cma";;
#load "models.cma";;
module W = Models.Wrightfisher;;
module S = Models.Setchains;;
module LL = Batteries.LazyList;;
*)
