PAPERNAME = manual


# FIGURES = figures/likelihood.pdf \
#           figures/5taxon_tree_rooted.pdf \
#           figures/5taxon_tree_unrooted.pdf \
#           figures/quadcore_shared_l3.pdf \
#           figures/quadcore_separate_l2.pdf \
#           figures/caneland.pdf \
#           figures/periscope.pdf \
#           figures/mpi_wrapper.pdf \
#           figures/openmp_scheme.pdf \
#           figures/mpi_scheme_simple.pdf \
#           figures/mpi_scheme_partitioned.pdf \
#           figures/mpi_scheme_hybrid.pdf

# PLOTS = plots/growth.pdf \
#         plots/spec_pin_nopin.pdf \
# 	plots/runtimes.SPEC.2.pdf \
#         plots/runtimes.SPEC.1-4.min_max.pdf \
# 	plots/runtimes.SPEC.4-16.min_max.pdf \
# 	plots/top500.NUM.pdf \
# 	plots/top500.PERF.pdf \
#         plots/runtimes.OMP.d50_50000.pdf \
#         plots/runtimes.OMP.d500_5000.pdf \
#         plots/speedups.OMP.d50_50000.pdf \
#         plots/speedups.OMP.d500_5000.pdf \
#         plots/runtimes.PT.d50_50000.pdf \
#         plots/runtimes.PT.d500_5000.pdf \
#         plots/speedups.PT.d50_50000.pdf \
#         plots/speedups.PT.d500_5000.pdf \
#         plots/runtimes.MPI.d50_50000.pdf \
#         plots/runtimes.MPI.d50_500000.pdf \
#         plots/runtimes.MPI.d250_500000.pdf \
#         plots/speedups.MPI.d50_50000.pdf \
#         plots/speedups.MPI.d50_500000.pdf \
#         plots/speedups.MPI.d250_500000.pdf \
#         plots/runtimes.HYBRID.multi.pdf \
#         plots/runtimes.HYBRID.32trees.pdf

#TEXFILES = \
#           introduction.tex 
#           shared_memory.tex \
#           distributed_memory.tex \
#           parallelization.tex \

#AUXFILES = \
#           phylogenetic_trees.aux
#           shared_memory.aux \
#           distributed_memory.aux \
#           parallelization.aux \

PDFCMD = pdflatex -interaction=nonstopmode

all: $(PAPERNAME).pdf 

clean:
	rm -rf *~ *.aux *.log *.bbl *.blg *.toc

$(PAPERNAME).pdf: $(PAPERNAME).aux $(AUXFILES) $(PAPERNAME).bbl $(FIGURES) $(PLOTS) 
	$(PDFCMD) $(PAPERNAME)

$(PAPERNAME).aux: $(PAPERNAME).tex $(TEXFILES) $(FIGURES) $(PLOTS)
	$(PDFCMD) $(PAPERNAME)

$(PAPERNAME).bbl: references.bib $(AUXFILES)
	$(PDFCMD) $(PAPERNAME)
	bibtex $(PAPERNAME)
	$(PDFCMD) $(PAPERNAME)

# figures/%.eps: figures/%.dia
#	dia $< -t eps -e $@

figures/%.pdf: figures/%.eps
	epstopdf $<

figures/%.pdf: figures/%.svg
	inkscape --export-pdf=$@ $<

plots/%.pdf: plots/%.gp plots/settings.gp
	gnuplot $< |epstopdf --filter >$@

plots/%.png: plots/%.gp plots/settings.gp
	gnuplot $< >$@

