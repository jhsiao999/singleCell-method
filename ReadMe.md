## Single cell data methods development
na

View our results [here](https://jhsiao999.github.io/singleCell-method).

### Creating a new anaysis

* Open Rstudio project `singleCell-method.Rproj`.
* Creat a copy of [template.Rmd][].
* Chnage the author, title, and date at the top of the file.
* Add the analysis code.
* Use the RStudio button "Preview HTML" to view the results.
* Add the analysis to the table of content in [index.md][].
* Add, commit, and push the new analysis.


### Builing the site

```bash
git checkout gh-pages
git merge master
git add index.md *.html
git commit -m "Build site"
git push origin gh-pages
```


### Acknowledgements

* Thanks to Karl Broman, who provided easy-to-follow instructions on buildling a gh-pages for beginners.
* Thanks to John Blischak for the contributing guidelines. 

