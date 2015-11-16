## Single cell data methods development

View our website [here](https://jhsiao999.github.io/singleCell-method).

## Outline

*  [Site organization](#site-organization)
*  [Creating a new analysis](#creating-a-new-analysis)
*  [Collaborating on the project](#collaborating-on-the-project)
*  [Building the site](#building-the-site)
*  [Style guide](#style-guide)
*  [Acknowledgements](#acknowledgements)


## Site organization

For consistency, we structure the content as follows:

*  `project/analysis`: New analysis (Rmd and html).
*  `project/analysis/rdas`: Interim or end-of-analysis results that are 
	saved for later analysis.
*  `project/analysis/figures`: Image files. 
*  `project/data`: Input data, such as count matrix and phenotypes.
*  `project/code`: Data processing or analysis codes.
*  `project/slides`: Talk slides.


## Creating a new anaysis

*  Go to the local repository. 
*  Switch to `master` branch.
*  Open Rstudio project `project/analysis/singleCell-method.Rproj`.
*  Creat a copy of [template.Rmd][].
*  Chnage the author, title, and date at the top of the file.
*  Add the analysis code.
*  Use the RStudio button "Preview HTML" to view the results.
*  Add the analysis to the table of content in [index.md][].
*  Add, commit, and push the new analysis.



## Collaborating on the project

Person A makes changes in the local branch and attemps to merge these changes into the remote master and update the remote master.

```bash
git checkout work-branch  ## move the pointer to local work branch
git add new_edits
git commit -m "new_edits" 
git push origin work-branch  ## push edits to remove work branch

### At this point, you can make a merge and a pull request to review the edits 

git checkout master  ## move the pointer to local master
git pull origin master  ## fetch and merge remote master to local master
git merge work-branch  ## merge local work-branch into master and update master 
git checkout work-branch  ## move to local work-branch
git merge master  ## merge remote master to local work-branch
git push origin work-branch  ## move the pointer back to local work-branch

### At this point, the commit numbers of your work-branch and master should be the same!!!!!!!!!!!!!!!!!!!
```


## Builing the site

```bash
git checkout gh-pages
git merge master
git push origin gh-pages
```



## Style guide

We will follow the useful Style Guide created by John Blischak for a
related [single cell project]()

*  Name variables using snakecase, e.g. `gene_exp_mat`.
*  Name files with dashes, e.g. `this-is-a-long-filename.txt`.
*  Name directories with camelCase, e.g `data`, `rawData`.
*  Use `<-` for assignment.
*  Surround binary operators with spaces, e.g. `1 + 1`, not `1+1`.
*  Use two spaces for indentation.
*  Lines should not be longer than 80 characters.

When in doubt, use the style indicated either in [Google's R Style Guide][google-style] or [Ha$

When writing text, aim to write one sentence per line.
This makes it easier to understand edits when reviewing the version control log.
The limit of 80 characters for code described above does not need to be applied to text.

[google-style]: https://google-styleguide.googlecode.com/svn/trunk/Rguide.xml
[hadley-style]: http://r-pkgs.had.co.nz/style.html




## Acknowledgements

*  Thanks to Karl Broman for maintaining publically accessible 
   instructions of building a github website.
*  Thanks to John Blischak for sharing his contributing guidelines.
*  Thanks to Raman Shah for helpful tips on Git collaboration.
