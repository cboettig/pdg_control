all:
	make install
	make render


install: components/install.R
	Rscript components/install.R

components/cache: components/restore-cache.R
	Rscript components/restore-cache.R

# combines knitr & pandoc calls
render:
	Rscript -e 'rmarkdown::render("manuscript.Rmd")'


clean:
	rm -f $(MANUSCRIPT).md *.log *.aux *.pdf  
	make clear-cache
	make clear-packrat


clear-packrat: 
	rm -rf .Rprofile packrat/lib packrat/src packrat/bootstrap.R

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



