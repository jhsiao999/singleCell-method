# ashlar: A workflow template for statistical computing projects

[*ashlar*](http://github.com/jhsiao999/ashlar) is our attempt to streamline workflow and to do reproducible research here at the University of Chicago [Stephens Lab](http://stephenslab.uchicago.edu/). 

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Making your own ashlar](#making-your-own-ashlar)
  - [Cloning [*ashlar*](http://github.com/jhsiao999/ashlar)](#cloning-ashlarhttpgithubcomjhsiao999ashlar)
  - [Reset git remote directory](#reset-git-remote-directory)
  - [Producing and publishing the website](#producing-and-publishing-the-website)
    - [Option 1: All contents for my eyes only](#option-1-all-contents-for-my-eyes-only)
    - [Option 2: Publish it! Keep a two-branch workflow.](#option-2-publish-it-keep-a-two-branch-workflow)
  - [A typical git workflow](#a-typical-git-workflow)
- [Resources](#resources)
- [Some tips on how to collaborate with others using GitHub](#tip-collaborate)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->




## Making your own ashlar

### Cloning [*ashlar*](http://github.com/jhsiao999/ashlar) 

[*ashlar*](http://github.com/jhsiao999/ashlar) is inspired by [*singleCellSeq*](https://github.com/jdblischak/singleCellSeq) - a collaborative project between biologists, bioinformaticians and statisticians that aim at exploring and understanding batch effects in single-cell RNA sequencing data. Both projects adopt the popular [*rmarkdown*](http://rmarkdown.rstudio.com/) website layout.

I suggest cloning into a new folder to distinguish your work from the example repository.

```
git clone https://github.com/jhsiao999/ashlar.git ashlar-trial
```



### Reset git remote directory 

At this point, your remote directory of the clone is still [*ashlar*](http://github.com/jhsiao999/ashlar). Make sure you change the name of the remote repository to match your local directory.

```
git remote rm origin
git remote add origin https://github.com/jhsiao999/ashlar-trial.git
```

Create a repository at [*github.com*](http://github.com). Then, push contents of the entire directory to the *master* branch. We use *git add -f* option to force add *html* files to the master branch, such as *index.html* for table of content. The default *.gitignore* in [*ashlar*](http://github.com/jhsiao999/ashlar) ignores *htmls*. 

```
git add -f --all
git commit -m "first commit"
git push origin master
```



### Personalize your repository

Now edit these files 



### Producing and publishing the website 

#### Option 1: All contents for my eyes only

Open index.html. This is the homepage of your unpubished website. You are DONE!

If you choose this option, you only have the master branch. The *gitignore* is set up to not to push *htmls, pngs, pdfs, etc* to the remote *master* brach, so edit the *.gitignore* to add these files if you want to add them to the remote directory. 


#### Option 2: Publish it! Keep a two-branch workflow.

Copy the master branch into a new branch named gh-pages. 

```
git checkout -b gh-pages 
```

Use knitr to turn Rmds into htmls. I mostly use RStudio to render Rmds. Sometimes *make*
file is useful for reproducing a large number of htmls or for rendering Rmds under
customized options.

```
cd analysis
make
```

Add files to git index and record changes in the local repository.

```
git add -f --all
git commit -m "Build site"
```

Update remote repository with the local changes.

```
git push origin gh-pages
```


The site address is under the analysis directory since the contents are under the analysis directory.

*https://jhsiao999.github.io/ashlar-trial/analysis*


### A typical git workflow

This two-branch workflow is set up to separate source codes from their output files.
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

* [John Blischak's tips on worflow management.][contrib]


[site]: http://jhsiao999.github.io/ashlar/analysis
[contrib]: https://github.com/jdblischak/singleCellSeq/blob/master/CONTRIBUTING.md


## Some tips on how to collaborate with others using GitHub

Many of us collaobrate with scientists from different discipline on code, analysis, or bioinformatic pipelines. There are two common options

1. Each collaborator creats his/her own work branch in the repository.

2. Each collaborator creates his/her own fork of the repository.

Both of these options allow the other collaborators to review changes in a pull request, comment on the changes, and merge the changes. Version control process is also similar in these two options. From my experience, the major difference is your work habits. When I work with a branch inside a repository, I merge the changes to the master or the gh-pages branch when done with the changes. This often results in me forgetting to change the branch back to my work branch, and the next time I work in the repository, I'll end up making changes to the master branch...

I suggest using the forking option. With this option, you are always working under your fork and don't EVER have to worry about working under the wrong branch. The steps are simple. 

How to fork a repo:

https://help.github.com/articles/fork-a-repo/


How to sync a fork with the remote repository"

https://help.github.com/articles/syncing-a-fork/









