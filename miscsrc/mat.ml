module M  = Owl.Mat
(* module I = Core.Interval;; *)

(* Something like this is what one wants, I think:
module Mat_interval = Interval.Make(struct 
  type t = Owl.Mat.mat
  let compare = Comparemat.compare
end);;
*)

(*
# module MC = Core.Comparable.Make(struct
  type t = Owl.Mat.mat
  let compare = Mat.compare end);;
Error: Signature mismatch:
       ...
       The value `sexp_of_t' is required but not provided
       The value `t_of_sexp' is required but not provided
# module MC = Core.Comparable.Make(Mat);;
Error: Signature mismatch:
       ...
       The value `sexp_of_t' is required but not provided
       The value `t_of_sexp' is required but not provided
*)

type t = Owl.Mat.mat

let compare = Utils.biased_compare


(* debug:
  let yo = ref 0 in
    yo := !yo + 1; Printf.printf "%d %d\n" acc !yo;
*)
