# 设置输出格式和文件
set terminal pngcairo size 1700,1600 enhanced 
set output 'cylindricalEquilibrium.png'

# 定义彩虹色调色板
#set palette defined ( 0 "blue", 1 "cyan" , 2 "green", 3 "yellow", 4 "red")

# 启用pm3d，设置视图和表面选项
set pm3d
set view map
unset surface
set size ratio -1
set cbtics scale 0

# 总标题
set multiplot layout 2,2 title "" font 'Verdana,26'

set zrange[0:2.5]
set cbrange[0:2]
# 图1
set title "Perfect Conducting" font 'Verdana,36' offset 0,2.2
set ylabel "Density" font 'Verdana,36' offset -4.0,0
set object 1 circle at 0,0 size 2 front lw 0.01 fc rgb "black"
splot 'eta00000ThirdVersionDC/tokamak_t200.dat' using 1:2:3 with dots title "" 
unset ylabel


# 图2
set title "Resistive eta=0.05" font 'Verdana,36' offset 0,2.2
set object 1 circle at 0,0 size 2 front lw 0.01 fc rgb "black"
splot 'eta00500ThirdVersionDC/tokamak_t200.dat' using 1:2:3 with dots title ""
unset cbrange
unset zrange


set cbrange[0:5.3]

# 图1
set title "" font 'Verdana,26'
set ylabel "Bz" font 'Verdana,36' offset -4.0,0
set object 1 circle at 0,0 size 2 front lw 0.01 fc rgb "black"
splot 'eta00000ThirdVersionDC/tokamak_t200.dat' using 1:2:11 with dots title "" 
unset ylabel


# 图2
set title "" font 'Verdana,26'
set object 1 circle at 0,0 size 2 front lw 0.01 fc rgb "black"
splot 'eta00500ThirdVersionDC/tokamak_t200.dat' using 1:2:11 with dots title ""
unset cbrange




unset multiplot
set output
