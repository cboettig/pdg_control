## Making supplement makes both, since it needs manuscript.Rmd env vars
all:
	make supplement 

# combines knitr & pandoc calls
manuscript: manuscript.Rmd
	Rscript -e 'library("rmarkdown"); render("manuscript.Rmd")'

# supplement needs manuscript in cache
supplement: supplement.Rmd
	Rscript -e 'library("rmarkdown"); render("manuscript.Rmd"); render("supplement.Rmd")'

clean:
	rm -f $(MANUSCRIPT).md *.log *.aux *.pdf  
	make clear-cache


