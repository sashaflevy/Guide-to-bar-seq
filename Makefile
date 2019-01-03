




# Here's a rule to just convert PDFs to PNGs for good sharing
output/%.png: output/%.pdf
	convert -density 300 $< $@

# Here's a rule to just convert PDFs to JPGs for easy sharing
output/%.jpg: output/%.pdf
	convert $< $@
