# ashlar

Workflow template for statistical computing projects at [Stephens Lab](http://stephenslab.uchicago.edu/). 

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Make your own ashlar](#make-your-own-ashlar)
  - [1. Making a copy of ashlar](#1-making-a-copy-of-ashlar)
  - [2. Configuring git settings](#2-configuring-git-settings)
  - [3. Editing basic information and site layout](#3-editing-basic-information-and-site-layout)
  - [4. Adding and committing local files](#4-adding-and-committing-local-files)
  - [5. Publishing the content](#5-publishing-the-content)
- [Adding new analysis](#adding-new-analysis)
- [What if I don't want to publish the content?](#what-if-i-dont-want-to-publish-the-content)
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

### 2. Configuring git settings

Create a repository for *ashlar-trial* on github.com. Then, link this remote 
directory with the local folder created in the first step.

```
cd GitHub-local-dir/ashlar-trial
git remote rm origin
git remote add origin https://github.com/jhsiao999/ashlar-trial.git
```

### 3. Editing basic information and site layout

Now that you have a copy of ashlar. Edit the following content and add your project information.

* analysis/index.Rmd: Homepage of the website, which is typically the table of content. Follow markdown syntax to make bullet lists.
* README.md: readme of the GitHub repository. 
* analysis/include/before_body.html: Edit the name of the repository in line 10. 
Edit the hyperlink of the reposotyr in line 17.
* analysis/license.Rmd: Change the software license if necessary.
* analysis/about.Rmd: Edit basic information about yourself or about the project.
 
### 4. Adding and committing local files

*git add --all* add all local files to the remote directory, except for files listed in *.gitignore*.

```
git add --all
git commit -m "first commit"
git push -u origin master
```

### 5. Publishing the content

Switch to a new branch *gh-pages*. By convention, GitHub publishes the content 
of this branch to a website. All content of the repository becomes public, even 
if the repository was private.

```
git checkout -b gh-pages 
```

Push to the remote *gh-pages* branch. This step also generates your website.

```
git push origin gh-pages
```

Note: The website is under the analysis directory:
*https://jhsiao999.github.io/ashlar-trial/analysis*


## Adding new analysis

The idea is to work from the master branch and then add to both master branch and gh-pages branch.
Using this workflow, gh-pages branch mirrors master branch. Advanced users can modify *.gitignore* 
and choose to not to publish selected files.

```
git checkout master
cd analysis
git add *Rmd *html figure/* (or git add --all)
git commit -m "add new analysis"
git push origin master

git checkout gh-pages
git merge master
git push origin gh-pages
```


## What if I don't want to publish the content?

Do not make a gh-pages branch. Just work with the master branch. This means you will
skip the step "Publishing the content". 

You can still view the site locally, which can be accessed via its
homepage "index.html".


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









