# 设置输出格式和文件
set terminal pngcairo size 1600,400 enhanced 
set output 'penetration_wall.png'

set title "Penetration on Resistive Wall" font 'Verdana,16'
set xlabel "x" font 'Verdana,16'
set ylabel "Bx" font 'Verdana,16'

set key font 'Verdana,10'

set samples 1000

set xrange[-0.001:0.007]
set yrange[-0.2:1.2]

set arrow from 0,-0.2 to 0,1.2 nohead lc rgb "black" lw 0.1

plot (x < 0 ? 1.0 : exp(-x/4.6e-4)) with lines lw 1.5 title "t=0", (x < 0 ? 1.0 : exp(-x*exp(-1.0/2)/4.6e-4)) with lines lw 1.5 title "t=1.0", (x < 0 ? 1.0 : exp(-x*exp(-2.0/2)/4.6e-4)) with lines lw 1.5 title "t=2.0", (x < 0 ? 1.0 : exp(-x*exp(-3.0/2)/4.6e-4)) with lines lw 1.5 title "t=3.0", (x < 0 ? 1.0 : exp(-x*exp(-4.0/2)/4.6e-4)) with lines lw 1.5 title "t=4.0" , (x < 0 ? 1.0 : exp(-x*exp(-5.0/2)/4.6e-4)) with lines lw 1.5 title "t=5.0", (x < 0 ? 1.0 : exp(-x*exp(-6.0/2)/4.6e-4)) with lines lw 1.5 title "t=6.0"


set output
