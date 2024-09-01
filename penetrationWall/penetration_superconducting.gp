# 设置输出格式和文件
set terminal pngcairo size 1600,400 enhanced 
set output 'penetration_superconducting.png'

set title "Penetration on Superconductor wall" font 'Verdana,18'
set xlabel "x" font 'Verdana,18'
set ylabel "Bx" font 'Verdana,18'

set key font 'Verdana,12'

set samples 1000

set xrange[-0.004:0.004]
set yrange[-0.2:1.2]

set arrow from 0,-0.2 to 0,1.2 nohead lc rgb "black" lw 0.1

plot (x < 0 ? 1.0 : exp(-x/4.6e-4)) with lines lw 1.5 title "lambda_{0}", (x < 0 ? 1.0 : exp(-x/8.0e-8)) with lines lw 1.5 title "lambda_{Nb3Sn}" , (x < 0 ? 1.0 : exp(-x/1.0e-7)) with lines lw 1.5 title "lambda_{NbTi}"


set output
