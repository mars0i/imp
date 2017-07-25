
let f m sd x = Stats.Pdf.gaussian (x +. m) sd in
let h = Plot.create "yo.pdf" in
for i=1 to 50 do 
  let i' = float i in
  let m, sd = i' /. 1000., i' /. 1000. in
  Plot.plot_fun ~h (f m sd) (-0.2) 0.1;
done;
Plot.output h;;

let gauss m sd x = Stats.Pdf.gaussian (x +. m) sd;;
let xs = Mat.(((sequential 1 1000) -$ 499.) /$ 100.);;
let ys = Mat.(map (gauss 0. 1.) xs);;
let h = Plot.create "yo.pdf" in Plot.plot ~h xs ys ; Plot.output h;;
