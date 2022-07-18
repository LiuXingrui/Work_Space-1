
#set term png
#set output "./forcopy/cla_BP.png"
                    
    set log x
   set log y
# remove any log-scaling
      unset label                            # remove any previous labels
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically

      set xlabel "bit error rate"
      set ylabel "logical rate"
#     set key 0.01,100
#  set label "Yield Point" at 0.003,260
#     set arrow from 0.0028,250 to 0.003,280
#    set xr [0.0:0.022]
#     set yr [0:325]
    set title "Classical BP"

plot  "./cla_data/M1_p" using 0:1 title 'n=96 parallel' with lines, \
	     "./cla_data/M1_s" using 0:1 title 'n=96 serial' with lines, \
           "./cla_data/M2_p" using 0:1 title 'n=204 parallel' with lines, \
	     "./cla_data/M2_s" using 0:1 title 'n=204 serial' with lines, \
	      "./cla_data/M3_p" using 0:1 title 'n=408 parallel' with lines, \
	     "./cla_data/M3_s" using 0:1 title 'n=408 serial' with lines, \
