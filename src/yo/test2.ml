open Hdf5_caml

let () =
  let m = Owl.Mat.uniform 10 10 in
  let f = H5.create_trunc "foo.h5" in
  H5.write_float_genarray f "m" m;
  H5.close f;

  let f = H5.open_rdonly "foo.h5" in
  let m' = H5.read_float_genarray f "m" Bigarray.c_layout in
  H5.close f;
  Printf.printf "original data = data from file: %b\n" (m = m')
    
