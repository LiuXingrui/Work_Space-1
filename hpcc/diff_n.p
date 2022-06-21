
set term png
 set output "test.png"
set multiplot
	
      set size 1,0.33
      set   autoscale                        # scale axes automatically
#     set log x
#     set log y
# remove any log-scaling
      unset label                            # remove any previous labels
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically

      set xlabel "n"
      set ylabel "fail rate"
#     set key 0.01,100
#  set label "Yield Point" at 0.003,260
#     set arrow from 0.0028,250 to 0.003,280
#    set xr [0.0:0.022]
#     set yr [0:325]
   	set title "p=0.1 range=0.9"
    set origin 0.0,0.66
   plot  "n_p0.100_serial_2.data" using 1:3 title 'serial range=0.9' with lines, \
	      "n_p0.100_parallel_2.data" using 1:3 title 'parallel  range=0.9' with lines
#	  "n_p0.100_range0_serial_2.data" using 1:3 title 'serial range=0' with lines, \
#     "n_p0.100_range0_parallel_2.data" using 1:3 title 'parallel range=0' with lines
	 

			set title "p=0.01 range=0.9"
		        set origin 0.0,0.33
           plot  "n_p0.010_serial_2.data" using 1:3 title 'serial  range=0.9' with lines, \
	      "n_p0.010_parallel_2.data" using 1:3 title 'parallel  range=0.9' with lines
#		  "n_p0.010_range0_serial_2.data" using 1:3 title 'serial range=0' with lines, \
#	      "n_p0.010_range0_parallel_2.data" using 1:3 title 'parallel range=0' with lines


		 	   set title "p=0.001 range=0.9"
		        set origin 0.0,0.0
		 plot  "n_p0.001_serial_2.data" using 1:3 title 'serial range=0.9' with lines, \
	      "n_p0.001_parallel_2.data" using 1:3 title 'parallel  range=0.9' with lines
#		   "n_p0.001_range0_serial_2.data" using 1:3 title 'serial range=0' with lines, \
#	      "n_p0.001_range0_parallel_2.data" using 1:3 title 'parallel range=0' with lines

	       
 unset multiplot 


