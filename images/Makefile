# dt000.pov: dual_topology.py
# 	python dual_topology.py

# st000.pov: single_topology.py
# 	python single_topology.py

all: out.gif

%000.png: %000.pov
	for i in $$(seq -f %03.0f 0 100); do \
		povray +I$*$${i}.pov +O$*$${i}.png +A0.5; \
	done

comb000.png: st000.png dt000.png
	for i in $$(seq -f %03.0f 0 100); do \
		convert +append st$${i}.png dt$${i}.png tmp.png; \
		convert -pointsize 50 label:'Single topology             vs             Dual topology'  tmp.png -pointsize 50 -gravity Center -append comb$${i}.png; \
	done

out.gif: comb000.png
	convert -delay 10 -loop 0 comb*.png comb.gif

clean:
	rm -f *.png *.pov *.gif
