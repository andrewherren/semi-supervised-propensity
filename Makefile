# Runs project code and builds paper / presentation
# Ideas taken from: https://robjhyndman.com/hyndsight/makefiles/

# User inputs
LOGDIR=log
RDIR=R
PAPERTEXFILE=reports/paper/semi-supervised-propensity
PRESTEXFILE=reports/presentation/Presentation-Slides

# Latex paper
$(PAPERTEXFILE).pdf:	$(PAPERTEXFILE).tex
	latexmk -pdf -cd -quiet $(PAPERTEXFILE)
	
# Presentation
$(PRESTEXFILE).pdf:	$(PRESTEXFILE).tex
	latexmk -pdf -cd -quiet $(PRESTEXFILE)

# R output
$(LOGDIR)/batch_R_script.Rout:	$(RDIR)/run_cached_files.R
	R CMD BATCH $(RDIR)/run_cached_files.R $(LOGDIR)/batch_R_script.Rout

# To build just the paper, run `make paper` in project directory
paper:	$(PAPERTEXFILE).pdf
# To build just the presentation, run `make presentation` in project directory
presentation:	$(PRESTEXFILE).pdf
# To run only the R simulations, run `make r` in project directory
r:	$(LOGDIR)/batch_R_script.Rout

# Build everything
.PHONY:	r	presentation	paper