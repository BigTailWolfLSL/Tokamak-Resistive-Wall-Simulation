set term pngcairo size 700,900

totalEnergyFile = 'integrated_totalEnergy.dat'

set pm3d 
set view map
unset surface
set size ratio -1
set xtics font "Verdana,24"
set ytics font "Verdana,24"
set cbtics font "Verdana,24"
set object 1 polygon from 0.500000,0.000000 to 0.500000,0.135849 to 0.500000,0.294340 to 0.500000,0.475472 to 0.500000,0.633962 to 0.500000,0.792453 to 0.500000,0.973585 to 0.510000,1.086792 to 0.522727,1.200000 to 0.545455,1.336000 to 0.568182,1.471698 to 0.613636,1.640189 to 0.681818,1.831321 to 0.772727,2.023000 to 0.810000,2.085000 to 0.863636,2.150943 to 0.930000,2.210000 to 1.000000,2.264151 to 1.136364,2.332075 to 1.318182,2.380000 to 1.477273,2.395000 to 1.636364,2.377358 to 1.795455,2.332075 to 1.954545,2.263000 to 2.136364,2.150000 to 2.318182,1.992453 to 2.500000,1.800000 to 2.659091,1.590000 to 2.772727,1.381132 to 2.863636,1.154717 to 2.931818,0.905660 to 2.977273,0.679245 to 3.000000,0.452830 to 3.000000,0.226415 to 3.000000,0.000000 to 3.000000,-0.226415 to 3.000000,-0.452830 to 2.977273,-0.679245 to 2.931818,-0.905660 to 2.863636,-1.154717 to 2.772727,-1.381132 to 2.659091,-1.590000 to 2.500000,-1.800000 to 2.318182,-1.992453 to 2.136364,-2.150000 to 1.954545,-2.263000 to 1.795455,-2.332075 to 1.636364,-2.377358 to 1.477273,-2.395000 to 1.318182,-2.380000 to 1.136364,-2.332075 to 1.000000,-2.264151 to 0.930000,-2.210000 to 0.863636,-2.150943 to 0.810000,-2.085000 to 0.772727,-2.023000 to 0.681818,-1.831321 to 0.613636,-1.640189 to 0.568182,-1.471698 to 0.545455,-1.336000 to 0.522727,-1.200000 to 0.510000,-1.086792 to 0.500000,-0.973585 to 0.500000,-0.792453 to 0.500000,-0.633962 to 0.500000,-0.475472 to 0.500000,-0.294340 to 0.500000,-0.135849 to 0.500000,0.000000 fc rgb "black" lw 0.1 front
#set palette defined ( 0 "blue", 1 "cyan" , 2 "green", 3 "yellow", 4 "red")
set zrange[0:5]
set cbrange[0:3]

do for [i=1:400] {
    #integratedTotalEnergy = system(sprintf("awk 'NR==%d {print $2}' %s", i, totalEnergyFile))

    #energy = real(integratedTotalEnergy)
    
    set output sprintf('tokamak_Density_%03d.png', i)
    
   #set title sprintf("Time=%.3f           Total Energy=%.3f", i*0.002, energy)
    set title sprintf("Time=%.3f", i*0.002) font "Verdana,32"

    splot sprintf('tokamak_t%d.dat', i*2) using 1:2:3 notitle

}


unset output
