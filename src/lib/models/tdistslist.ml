(** Module to wrap lazy list functions so that they can be altered.
    Began with Batteries.LazyList, but will change to Core.Sequence. *)

module LL = Batteries.LazyList

type t = Tdists.t LL

let hd = LL.hd
let tl = LL.tl
let last = LL.last
let at = LL.at
let fold_left = LL.fold_left
let map = LL.map
let seq = LL.seq
let cons = LL.cons
let fold_right = LL.lazy_fold_right
let take = LL.take
let take_while = LL.take_while
let to_list = LL.to_list
let rev = LL.rev
