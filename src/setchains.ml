(** Functions inspired by Hartfiel's _Markov Set-Chains_, Springer 1998 *)

module L = Batteries.List

let vertices p q =
  let idxs = L.range 0 `To (L.length p) in
  let verts free =
    let unfree = L.remove idxs free in
    L.concat (L.map (fun i -> [L.at p i; L.at q i]) unfree

  L.concat 
  (L.init



