
output/%.png: output/%.pdf
	convert -density 300 $< $@
