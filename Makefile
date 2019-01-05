
.PHONY: all simulations clean
all: simulations

testsim: ./scripts/nextflow scripts/run_simulations.nf scripts/run_simulations.nfconfig
	$(word 1,$^) $(word 2,$^) -c $(word 3,$^) 

simulations: ./scripts/nextflow scripts/run_simulations.nf scripts/run_simulations.nfconfig
	$(word 1,$^) $(word 2,$^) -c $(word 3,$^)

# Here's a rule to just convert PDFs to PNGs for good sharing
output/%.png: output/%.pdf
	convert -density 300 $< $@

# Here's a rule to just convert PDFs to JPGs for easy sharing
output/%.jpg: output/%.pdf
	convert $< $@

clean:
	rm -r work .nextflow* reports tmp
