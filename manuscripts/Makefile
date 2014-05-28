# Fill these out
MANUSCRIPT=manuscript
BIBLIO=components/references.bib
CSL=components/ecology.csl
TEMPLATE=components/elsarticle.latex

all:
	make render
#	make install
#	make restore-cache
#	make $(MANUSCRIPT).pdf
#	make $(MANUSCRIPT).tex

render:
	Rscript -e 'rmarkdown::render("manuscript.Rmd")'

install: components/install.R
	Rscript components/install.R

$(MANUSCRIPT).md: $(MANUSCRIPT).Rmd
	Rscript components/knit.R $(MANUSCRIPT).Rmd	

## Customize template, latex engine, etc here
$(MANUSCRIPT).pdf: $(MANUSCRIPT).md
	pandoc -o $(MANUSCRIPT).pdf --template=$(TEMPLATE) --biblio $(BIBLIO) --csl $(CSL) --latex-engine xelatex $(MANUSCRIPT).md

$(MANUSCRIPT).tex: $(MANUSCRIPT).md
	pandoc -o $(MANUSCRIPT).tex --template=$(TEMPLATE) --biblio $(BIBLIO) --csl $(CSL) $(MANUSCRIPT).md



# Only runs if cache has been cleared
components/cache: components/restore-cache.R
	Rscript components/restore-cache.R

# just an alias 
restore-cache: 
	make components/cache
	
clear-cache:
	rm -rf components/cache

clean:
	rm -f $(MANUSCRIPT).md 


# Non portable, just for me to update the source cache
# Must be on same line to use cd. must use cd to get paths right
# when unpacking.  
deploy-cache:
	cd components; tar -pczf pdg_controlcache.tar.gz cache; scp pdg_controlcache.tar.gz ten:~/public_html/data/



