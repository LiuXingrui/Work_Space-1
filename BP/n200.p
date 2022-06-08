      set multiplot
      set size 1,0.5 
      set   autoscale                        # scale axes automatically
      unset log                              # remove any log-scaling
      unset label                            # remove any previous labels
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set title "n=200 parallel"
      set xlabel "p"
      set ylabel "fail rate"
#     set key 0.01,100
#  set label "Yield Point" at 0.003,260
#     set arrow from 0.0028,250 to 0.003,280
#    set xr [0.0:0.022]
#     set yr [0:325]
    set origin 0.0,0.5
   plot   "n200_range0.0_parallel2.data" using 4:3 title 'range=0' with lines, \
	      "n200_range0.1_parallel2.data" using 4:3 title 'range=0.1' with lines,\
	      "n200_range0.2_parallel2.data" using 4:3 title 'range=0.2' with lines,\
	         "n200_range0.3_parallel2.data" using 4:3 title 'range=0.3' with lines,\
		   "n200_range0.6_parallel2.data" using 4:3 title 'range=0.6' with lines,\
		     "n200_range0.9_parallel2.data" using 4:3 title 'range=0.9' with lines

 set title "n=200 serial"
		        set origin 0.0,0.0
		         plot   "n200_range0.0_serial2.data" using 4:3 title 'range=0' with lines, \
	      "n200_range0.1_serial2.data" using 4:3 title 'range=0.1' with lines,\
	      "n200_range0.2_serial2.data" using 4:3 title 'range=0.2' with lines,\
	         "n200_range0.3_serial2.data" using 4:3 title 'range=0.3' with lines,\
		   "n200_range0.6_serial2.data" using 4:3 title 'range=0.6' with lines,\
		     "n200_range0.9_serial2.data" using 4:3 title 'range=0.9' with lines



		       
 unset multiplot 
