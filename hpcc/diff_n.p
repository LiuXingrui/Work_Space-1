     set multiplot
      set size 1,0.3
      set   autoscale                        # scale axes automatically
#     set log x
#     set log y
# remove any log-scaling
      unset label                            # remove any previous labels
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
	 set title "p=0.1 range=0.9"
      set xlabel "n"
      set ylabel "fail rate"
#     set key 0.01,100
#  set label "Yield Point" at 0.003,260
#     set arrow from 0.0028,250 to 0.003,280
#    set xr [0.0:0.022]
#     set yr [0:325]
    set origin 0.0,0.75
   plot  "n_p0.100f_serial_2.data" using 2:3 title 'serial' with lines, \
	      "n_p0.100f_parallel_2.data" using 2:3 title 'parallel' with lines
	 

			set title "p=0.01 range=0.9"
		        set origin 0.0,0.5
           plot  "n_p0.010f_serial_2.data" using 2:3 title 'serial' with lines, \
	      "n_p0.010f_parallel_2.data" using 2:3 title 'parallel' with lines

		 	   set title "p=0.001 range=0.9"
		        set origin 0.0,0.25
           plot  "n_p0.001f_serial_2.data" using 2:3 title 'serial' with lines, \
	      "n_p0.001f_parallel_2.data" using 2:3 title 'parallel' with lines

		       
 unset multiplot 
