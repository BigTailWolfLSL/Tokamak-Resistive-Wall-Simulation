# 设置输出格式和文件
set terminal pngcairo size 1600,1000 enhanced 
set output 'TotalEnergy.png'

# 总标题
set title "Total Energy" font 'Verdana,26'

set xrange[0:0.2]

# 定义自定义虚线
set style line 1 lc rgb "magenta" lt 1 lw 2 dashtype (10,10)

set xtics font "Verdana,18"
set ytics font "Verdana,18"

plot 'eta00000FirstVersionDC/integrated_totalEnergy.dat' using 1:2 with lines title "eta=0.0" linestyle 1, \
     'eta00010ThirdVersionDC_low/integrated_totalEnergy.dat' using 1:2 with lines title "eta=0.001", \
     'eta00020ThirdVersionDC_low/integrated_totalEnergy.dat' using 1:2 with lines title "eta=0.002", \
     'eta00030ThirdVersionDC_low/integrated_totalEnergy.dat' using 1:2 with lines title "eta=0.003", \
     'eta00040ThirdVersionDC_low/integrated_totalEnergy.dat' using 1:2 with lines title "eta=0.004", \
     'eta00050ThirdVersionDC_low/integrated_totalEnergy.dat' using 1:2 with lines title "eta=0.005", \
     'eta00060ThirdVersionDC_low/integrated_totalEnergy.dat' using 1:2 with lines title "eta=0.006", \
     'eta00070ThirdVersionDC_low/integrated_totalEnergy.dat' using 1:2 with lines title "eta=0.007", \
     'eta00080ThirdVersionDC_low/integrated_totalEnergy.dat' using 1:2 with lines title "eta=0.008" lc rgb "#8B4513", \
     'eta00090ThirdVersionDC_low/integrated_totalEnergy.dat' using 1:2 with lines title "eta=0.009" lc rgb "#4682B4"

unset multiplot
set output
