# ashlar

Workflow template for statistical computing projects at [Stephens Lab](http://stephenslab.uchicago.edu/). 

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Make your own ashlar](#make-your-own-ashlar)
  - [1. Making a copy of ashlar](#1-making-a-copy-of-ashlar)
  - [2. Configuring settings](#2-configuring-settings)
  - [3. Adding and committing local files](#3-adding-and-committing-local-files)
  - [4. Publishing the content](#4-publishing-the-content)
  - [5. Adding new analysis](#5-adding-new-analysis)
- [Resources](#resources)
- [Collaborate with other using git](#collaborate-with-other-using-git)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->



## Make your own ashlar

Below we give a step-by-step guide for how to set up an ashlar repository, 
and how to publish the repository content to a website.

### 1. Making a copy of ashlar

Move to a directory which is going to keep a local copy of the ashlar repository. 
*git clone* to copy the remote ashlar folder to the local directory, and rename to
*ashlar-trial*.

```
cd GitHub-local-dir
git clone https://github.com/jhsiao999/ashlar.git ashlar-trial
```

### 2. Configuring settings

Create a repository for *ashlar-trial* on github.com. Then, link this remote 
directory with the local folder created in the first step.

```
git remote rm origin
git remote add origin https://github.com/jhsiao999/ashlar-trial.git
```

### 3. Personalizing ashlar

* analysis/index.Rmd: this is your table of content. Follow markdown syntax to make bullet lists.


### 4. Adding and committing local files

*git add -f --all* forces add all local files ot the remote directory. *-f* force option
overrides *.gitignore*. Following the steps below, you can now see the local files on github.com.

```
git add -f --all
git commit -m "first commit"
git push -u origin master
```

### 5. Publishing the content

Switch to a new branch *gh-pages*. By convention, GitHub publishes the content 
of this branch to a website.

```
git checkout -b gh-pages 
```

Add and commit changes to *gh-pages* branch.

```
git add -f --all
git commit -m "Build site"
```

Push to the remote *gh-pages* branch. This step also generates your website.

```
git push origin gh-pages
```

Note: The website is under the analysis directory:
*https://jhsiao999.github.io/ashlar-trial/analysis*


### 6. Adding new analysis

This workflow is set up to separate source codes from their output files. 
The source codes are kept in the master branch, and their output is stored and published
in the gh-pages branch.

```
git checkout master
cd analysis
git add new-analysis.Rmd index.Rmd
git commit -m "add new analysis"
git push origin master

git checkout gh-pages
git merge master
make
git add *Rmd *html figure/*
git commit -m "add new analysis"
git push origin gh-pages
```


## Resources 

* [singleCellSeq][singleCellSeq] inspires the creation of *ashlar*. For advanced users,
check out [singleCellSeq][singleCellSeq] for tips on writing a paper with markdown.


[singleCellSeq]: https://github.com/jdblischak/singleCellSeq
[site]: http://jhsiao999.github.io/ashlar/analysis
[contrib]: https://github.com/jdblischak/singleCellSeq/blob/master/CONTRIBUTING.md


## Collaborate with other using git

We suggest that each collaborator creates his/her own master copy of the repository - 
that is, forking the repository. Then, all changes to the main copy need to go through
a pull request, which can be reviewed by all collaborators.

How to fork a repo: https://help.github.com/articles/fork-a-repo/

How to sync a fork with the remote repository: https://help.github.com/articles/syncing-a-fork/









