set term pngcairo size 700,900

totalEnergyFile = 'integrated_totalEnergy.dat'

set pm3d 
set view map
unset surface
set size ratio -1
set object 1 circle at 0,0 size 2 front lw 0.01 fc rgb "black"
set cbrange[0:3]

do for [i=1:150] {
    integratedTotalEnergy = system(sprintf("awk 'NR==%d {print $2}' %s", i, totalEnergyFile))

    energy = real(integratedTotalEnergy)
    
    set output sprintf('tokamak_Bz_%03d.png', i)
    
    set title sprintf("Time=%.3f           Total Energy=%.3f", i*0.002, energy)

    splot sprintf('tokamak_t%d.dat', i*2) using 1:2:11 notitle

}


unset output
