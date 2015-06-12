## Making supplement makes both, since it needs manuscript.Rmd env vars
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


