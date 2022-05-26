      set multiplot
      set size 1,0.5 
      set   autoscale                        # scale axes automatically
      unset log                              # remove any log-scaling
      unset label                            # remove any previous labels
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set title "n=200 parallel, try_again if reach max iter"
      set xlabel "p"
      set ylabel "fail rate"
#     set key 0.01,100
#  set label "Yield Point" at 0.003,260
#     set arrow from 0.0028,250 to 0.003,280
#    set xr [0.0:0.022]
#     set yr [0:325]
    set origin 0.0,0.5
      plot   "200_0.0_p2.data" using 4:3 title 'range=0' with lines, \
	     "200_0.1_p2.data" using 4:3 title 'range=0.1' with lines,\
	      "200_0.2_p2.data" using 4:3 title 'range=0.2' with lines,\
		"200_0.3_p2.data" using 4:3 title 'range=0.3' with lines,\
		  "200_0.6_p2.data" using 4:3 title 'range=0.6' with lines,\
		    "200_0.9_p2.data" using 4:3 title 'range=0.9' with lines

 set title "n=200 serial try_again, if reach max iter"
		        set origin 0.0,0.0
		          plot   "200_0.0_s.data" using 4:3 title 'range=0' with lines, \
	     "200_0.1_s2.data" using 4:3 title 'range=0.1' with lines,\
	      "200_0.2_s2.data" using 4:3 title 'range=0.2' with lines,\
		"200_0.3_s2.data" using 4:3 title 'range=0.3' with lines,\
		  "200_0.6_s2.data" using 4:3 title 'range=0.6' with lines,\
		    "200_0.9_s2.data" using 4:3 title 'range=0.9' with lines




		       
 unset multiplot 
