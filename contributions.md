# How to contribute to this repository
## Documentation editing

To contribute to the documentation (and especially the **specifications.md**), you can edit directly the files in your browser. After any modifications, you have to commit your changes by giving a short description of the modifications made. 


## Code
If you want to contribute significantly to the code, please follow the next steps in order to keep things clean and ordered:

### If it's your first contribution:

* Clone this repository to your local machine
* Create your own git branch by running the command:

>
git checkout -b *nameOfYourNewBranch*

* Do your modifications (try to do one type of modification per commit) following the good practice rules defined below.
* Edit the .gitignore file to ignore some useless files in your local folder
* Commit the changes:

>
git commit --message="Changed temperature of LN2"

* Push your changes to the new branch:

>
git push origin *nameOfYourNewBranch*

* When you think your changes are ready to be included in the main code, open a pull-request by using the *New Pull-Request* button on the web interface and select `to:master` and `from:nameOfYourBranch`.

### If it's not your first contribution:

* Pull the current state of the repository:

>
git pull origin *nameOfYourNewBranch*

* Do your modification and commit the changes:

>
git commit --message="Changed temperature of LN2"

* Push your changes to your branch:

>
git push origin *nameOfYourNewBranch*

* When you think your changes are ready to be included in the main code, open a pull-request by using the *New Pull-Request* button on the web interface and select `to:master` and `from:nameOfYourBranch`.

### Send the code by email
If this is just some small changes, you can also just send me the patches by email and I can integrate it directly. 

## Good practice for code editing in Matlab
In order to have some coherence for this new retrieval code, here are some basic rules I would suggest to use when collaborating to the project.

### Modules
Every module has to make 1 (and only 1) action at a time and be named accordingly. Suggestions for naming the modules:

*module_doing_this_action(var1,var2,...)*

### Variables 
No hardcoded variables, except in the main module. Suggestions for naming the variables:

*nameOfTheVariable*

## Good practice for code editing in Python
 TODO
