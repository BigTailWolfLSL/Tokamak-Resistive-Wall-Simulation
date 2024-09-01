# 设置输出格式和文件
set terminal pngcairo size 2400,1600 enhanced 
set output 'InitialCylindricalEquilibrium.png'

# 设置坐标轴和 color box 刻度的字体和大小
set xtics font "Verdana,22"
set ytics font "Verdana,22"
set ztics font "Verdana,22"
set cbtics font "Verdana,22"  # 调整 color box 刻度的字体大小

# 启用pm3d，设置视图和表面选项
set pm3d
set view map
unset surface
set size ratio -1
set cbtics scale 0

# 总标题
set multiplot layout 2,3 title "Initial States" font 'Verdana,40' 

# 图1
set title "Density" font 'Verdana,40' offset 0,1.0
set object 1 circle at 0,0 size 2 front lw 0.01 fc rgb "black"
set zrange[0:1.0]
splot 'tokamak_t.dat' using 1:2:3 with dots title "" 
unset zrange

# 图2
set title "Bz" font 'Verdana,40' offset 0,1.0
set object 1 circle at 0,0 size 2 front lw 0.01 fc rgb "black"
splot 'tokamak_t.dat' using 1:2:11 with dots title ""

# 图3
set title "Vx" font 'Verdana,40' offset 0,1.0
set object 1 circle at 0,0 size 2 front lw 0.01 fc rgb "black"
set cbrange[-2:10]
splot 'tokamak_t.dat' using 1:2:4 with dots title ""
unset cbrange

# 图4
set title "Bx" font 'Verdana,40' offset 0,1.0
set object 1 circle at 0,0 size 2 front lw 0.01 fc rgb "black"
set cbrange[-0.2:0.1]
splot 'tokamak_t.dat' using 1:2:9 with dots title ""
unset cbrange

# 图5
set title "By" font 'Verdana,40' offset 0,1.0
set object 1 circle at 0,0 size 2 front lw 0.01 fc rgb "black"
set cbrange[-0.2:0.1]
splot 'tokamak_t.dat' using 1:2:10 with dots title ""
unset cbrange

# 图6
set title "Vy" font 'Verdana,40' offset 0,1.0
set object 1 circle at 0,0 size 2 front lw 0.01 fc rgb "black"
set cbrange[-2:10]
splot 'tokamak_t.dat' using 1:2:5 with dots title ""
unset cbrange

unset multiplot
set output
