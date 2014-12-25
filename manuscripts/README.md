Instructions for compiling manuscripts
======================================

**Quickstart**:

Edit the `manuscript.Rmd` file. It contains the text, code, and various
formatting options.  Then to build the manuscript pdf, use at the (shell)
command line:

```bash
make
```

Or from R:

```r
rmarkdown::render("manuscript.Rmd")
```



Installation requirements
-------------------------


System requirements

- R (>= 3.0)
- pandoc (>= 1.12), pandoc-citeproc
- LaTeX environment

Installing [RStudio >= 0.98b] should provide everything necessary on most platforms.

Manuscripts in this package are written in [knitr]'s R Markdown format
(`.Rmd` files), and compiled into PDFs using [pandoc] with [LaTeX]. A
`.tex` file is also created (e.g. for journal submission systems). This
workflow can easily produce alternative formats, such as Microsoft
Office's `.docx`, HTML, etc.

Because markdown is platform independent plain text and easier to
learn than LaTeX, the hope is that this workflow is relatively portable
across users.  Markdown also does a better job than most alternatives
(tex included) at separating content from formatting, freeing the writer
to just write.

[RStudio >= 0.98b]: http://www.rstudio.com/ide/download/preview
[knitr]: http://yihui.name/knitr
[pandoc]: http://johnmacfarlane.net/pandoc/
[LaTeX]: http://www.latex-project.org/
[components/config_pandoc.txt]: http://github.com/cboettig/template/tree/master/manuscripts/components/config_pandoc.txt

Caching
-------

I've also enabled caching.  It can be annoying to have to have to rerun
all the R code just to make a textual change to the manuscript.

The first time you run `make` the cache will be downloaded from a
remote archive (since it is ignored by `git`).  To have the code
re-run locally, just delete the cache (`make clear-cache`) and
then recreate the document from scratch with `make render`.

The caching is modestly intelligent, in that if you edit a chunk it will be rerun by
default, (as will chunks with declared dependencies on it). See `knitr`'s
[caching documentation] for details.

**Restoring the default cache**

You can obtain my current cache by using `make restore-cache`.

**Interactive R**

For the sake of modular reproducibility, it can also be desirable
to pick up from somewhere in the middle of the manuscript and not
want to have to run all the previous code.  The easiest way to load
the pre-computed results is to launch an R session and call

```r
rmarkdown::render("manuscript.Rmd")
```

The manuscript will be rebuilt and all variables will be available
in the R environment.  Caching is also modular by chunk, allowing you
to restore the results only up to a particular chunk to investigate.
This can be very useful if a given variable is changed by later chunks.
Again see `knitr`'s [caching documentation] for details.


[caching documentation]: http://yihui.name/knitr/demo/cache/
