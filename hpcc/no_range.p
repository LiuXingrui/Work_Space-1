set format x "10^{%L}"
set format y "10^{%L}"
set key left
set term png
set output "./forcopy/no_range_BP.png"
                    
    set log x
   set log y
# remove any log-scaling
      unset label                            # remove any previous labels
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically

      set xlabel "Bit-error rate"
      set ylabel "Logical-error rate"
#     set key 0.01,100
#  set label "Yield Point" at 0.003,260
#     set arrow from 0.0028,250 to 0.003,280
#    set xr [0.0:0.022]
#     set yr [0:325]
    set title "Quantum BP"

plot  "./no_range_data/d5_parallel.data" using 4:3 title 'd=5 toric parallel' with linespoints, \
   "./no_range_data/d5_serial.data" using 4:3 title 'd=5 toric serial' with linespoints, \
      "./no_range_data/d7_parallel.data" using 4:3 title 'd=7 toric parallel' with linespoints, \
	 "./no_range_data/d7_serial.data" using 4:3 title 'd=7 toric serial' with linespoints, \
	    "./no_range_data/d9_parallel.data" using 4:3 title 'd=9 toric parallel' with linespoints, \
	       "./no_range_data/d9_serial.data" using 4:3 title 'd=9 toric serial' with linespoints,\
		  2*x title "y=2x" 



