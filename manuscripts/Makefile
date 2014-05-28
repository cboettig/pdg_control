all:
	make install
	make restore-cache
	make render


install: components/install.R
	Rscript components/install.R

# just an alias 
restore-cache: 
	make components/cache

# combines knitr & pandoc calls
render:
	Rscript -e 'rmarkdown::render("manuscript.Rmd")'

# Set target as components/cache, so it only runs if cache has been cleared
components/cache: components/restore-cache.R
	Rscript components/restore-cache.R

clear-cache:
	rm -rf components/cache

clean:
	rm -f $(MANUSCRIPT).md *.log *.aux  



# Non portable, just for me to update the source cache
# Must be on same line to use cd. must use cd to get paths right
# when unpacking.  
deploy-cache:
	cd components; tar -pczf pdg_controlcache.tar.gz cache; scp pdg_controlcache.tar.gz ten:~/public_html/data/



