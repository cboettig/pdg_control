## Making supplement makes both
all:
	make supplement.pdf 


# combines knitr & pandoc calls
manuscript.pdf: manuscript.Rmd
	Rscript -e 'library("rmarkdown"); render("manuscript.Rmd")'

# supplement needs manuscript in cache
supplement.pdf: supplement.Rmd
	Rscript -e 'library("rmarkdown"); render("manuscript.Rmd"); render("supplement.Rmd")'


clean:
	rm -f $(MANUSCRIPT).md *.log *.aux *.pdf  
	make clear-cache


## Stuff for dealing with caching ##  

components/cache: components/restore-cache.R
	Rscript components/restore-cache.R

clear-cache:
	rm -rf components/cache

# just an alias 
restore-cache: 
	make components/cache

# Non portable, just for me to update the source cache
# Must be on same line to use cd. must use cd to get paths right
# when unpacking.  
deploy-cache:
	cd components; tar -pczf pdg_controlcache.tar.gz cache; scp pdg_controlcache.tar.gz ten:~/public_html/data/

deploy-presentation:
	cp presentation.html ../../gh-pages/pdg_control/



	

