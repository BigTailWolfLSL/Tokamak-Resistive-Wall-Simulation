set term pngcairo size 700,900

totalEnergyFile = 'integrated_totalEnergy.dat'

set pm3d 
set view map
unset surface
set size ratio -1
set object 1 circle at 0,0 size 2 front lw 0.01 fc rgb "black"
set zrange[0:2.5]
set cbrange[0:2.2]

do for [i=1:100] {
    integratedTotalEnergy = system(sprintf("awk 'NR==%d {print $2}' %s", i, totalEnergyFile))

    energy = real(integratedTotalEnergy)

    set output sprintf('tokamak_Density_%03d.png', i)
    
    set title sprintf("Time=%.3f           Total Energy=%.3f", i*0.002, energy)

    splot sprintf('tokamak_t%d.dat', i*2) using 1:2:3 notitle

}


unset output