reset

unset title
unset key
#unset surface

#set terminal postscript eps enhanced font "Times-Roman" 24
#set output "permeabilityField.eps"

#set terminal pngcairo transparent enhanced font "arial,10" fontscale 1.0 size 600, 400
#set output "permeabilityField.png"
set view 60, 75, 1, 1
#set view map scale 1

#set dgrid3d 30,27,5
#set dgrid3d qnorm 2
#set dgrid3d 30,27 splines
#set dgrid3d cauchy
#set dgrid3d gauss 0.350000,0.500000
set dgrid3d 30,27 gauss 0.01,0.01

#set pm3d implicit at b
#set samples 20, 20
#set isosamples 21, 21
set contour base
set cntrlabel format '%8.3g' font 'arial,14' start 0 interval 50
set cntrparam order 8
set cntrparam bspline
#set cntrparam levels 10
#set cntrparam levels discrete 4e-16, 2e-15, 6e-15, 2e-14, 4e-13
set cntrparam levels incremental 0,1e-14,1e-13

set hidden3d back offset 1 trianglepattern 3 undefined 1 altdiagonal bentover
set title "{/Symbol k}, [m^2]"

set xlabel "{/Symbol q\260}" offset 0, 0, 0
set xtics offset 0, 0.5 rotate by 90 right
set xrange [*:*]

set ylabel "{/Symbol d}, [m]" offset 0, 0, 0
set ytics offset -0.5, 0 rotate by 0 right
set yrange [0.01:*]

#set logscale z

splot "permeabilitiesFirstDirection" using 1:2:3 with lines lw 5, "permeabilitiesFirstDirection" using 1:2:3 with labels
pause -1
