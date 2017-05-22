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

(** A compare function for matrices that returns zero if all elements of
    both matrices are equal, and if not returns -1 only if all elements 
    of m1 are less than or equal to corresponding elements of m2; otherwise
    returns 1, indicating that at least one element in m1 is greater than 
    the corresponding element in m2. *)
let compare m1 m2 =
  let f acc e1 e2 =
    if acc = 1 || e1 > e2 then 1  (* at least one pair has had e1 > e2 *)
    else match acc with           (* at this point we know that e1 <= e2 *)
         | -1 -> -1               (* all previous pairs were <= *)
         |  0 -> if e1 = e2 then 0 else -1
         |  _ -> raise (Failure "bug in compare: acc is not -1, 0, or 1")
  in Utils.fold2 f 0 m1 m2
