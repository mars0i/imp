module H = Hdf5_caml.H5

module M = Owl.Mat

let () = 
  print_string "creating simple array\n";
  let m = [| 0.; 1.; 2.; 3. |] in
    let f = H.create_trunc "yow.h5" in
    print_string "writing data to file yow.h5\n";
      H.write_float_array f "m" m;
      H.close f;
    print_string "reading from file yow.h5\n";
    let f = H.open_rdonly "yow.h5" in
    let m' = H.read_float_array f "m" in
      H.close f;
      Printf.printf "original data is the same as data in file: %b\n\n" (m = m')

let () = 

  print_string "creating 10x10 random matrix\n";
  let m = M.uniform 10 10 in
  print_string "opening file foo.h5 for writing\n";
  let f = H.create_trunc "foo.h5" in
  print_string "writing data to file foo.h5\n";
  H.write_float_genarray f "m" m;
  H.close f;
  print_string "opening file foo.h5 for reading\n";
  let f = H.open_rdonly "foo.h5" in
  print_string "reading data from file foo.h5\n";
  let m' = H.read_float_genarray f "m" Bigarray.c_layout in
  print_string "closing file\n";
  H.close f;
  Printf.printf "original data is the same as data in file: %b\n" (m = m')
    

(*
open Hdf5_caml

let () =
  let h5 = H5.create_trunc "test.h5" in
  H5.write_float_array h5 "a" [| 0.; 1.; 2.; 3. |];
  let a = H5.open_dataset h5 "a" in
  H5.write_attribute_float a "b" 4.;
  H5.close h5
*)
