# Contributing to SPECFEM3D

You want to contribute to the SPECFEM code? Great, let's check out this set of simple guidelines to follow for contributions.

## Contributing as a user

Software needs constant maintaining and updating to keep up with newest hardware and compilers.
You found an issue when running the code? Please consider creating a "*New issue*" on the [SPECFEM3D Issues github page](https://github.com/SPECFEM/specfem3d/issues).

Describe the problem as detailed as possible to help us reproduce the error. You can also attach console outputs from the executables to make it easier for debugging. Feel free to add an appropriate label to highlight the issue by checking out "*Labels*" on the right side of the page. Once done, click the button "*Submit new issue*". Catching bugs is always welcome, many thanks to you for improving the code!

## Contributing as a developer

You have a new feature, bug fix, or other modification you want to contribute to the code? In this case, consider submitting a "*Pull request*" to the **devel branch** of our github code repository.

This will require a few steps to setup your own github fork and be able to create a PR to the official devel version of the SPECFEM code (**note that only Pull requests towards devel are accepted**). The most basic setup looks the following:

#### 1. Create your fork of the repository:
Go to the main [SPECFEM3D github repository](https://github.com/SPECFEM/specfem3d) and click the "*Fork*" button at the top of the page. This will create a copy of the SPECFEM3D repository in your personal GitHub account.

#### 2. Clone your fork to your local workstation/laptop:
```
git clone --recursive --branch devel https://github.com/<your-github-account-name>/specfem3d.git
```
Once you change into your local folder `cd specfem3d/`, all git commands will be recognized.  
Now, add the remote address of the SPECFEM3D repository:
```
git remote add upstream https://github.com/SPECFEM/specfem3d.git
```
With `git remote -v`, you should see now an output like:
```
origin  https://github.com/<your-github-account-name>/specfem3d.git (fetch)
origin  https://github.com/<your-github-account-name>/specfem3d.git (push)
upstream  https://github.com/SPECFEM/specfem3d.git (fetch)
upstream  https://github.com/SPECFEM/specfem3d.git (push)
```
Please make sure that your local version is up-to-date with the upstream version. To so, you can either call
```
git pull upstream devel
```
or
```
git fetch upstream
git rebase --interactive upstream/devel
```

#### 3. Apply your modifications to your local code and commit your changes:
Check first the status of the files, and/or differences applied before committing, by typing
```
git status
```
and
```
git diff
```
Then, create the commit with a short descriptive message by
```
git commit -m "<my-commit-description-message>"
```
Remember that it is good practice to have small commits with a single purpose.

#### 4. Push your modifications to your github repository:
```
git push
```
or more explicit
```
git push origin devel
```

#### 5. Open a new "*Pull request*" on the github website:
Go to your personal github repository website, and possibly click on the branch button to change to the **devel** branch.
Github might already recognize that you have new changes uploaded, sees that your branch is "*ahead of SPECFEM:devel*" and provides you with the option to "*Contribute*" by clicking on the "*Open pull request*" button.<br>

In the "*Open a pull request*" page, double-check that the base repository is: **SPECFEM/specfem3d, base: devel**. The head repository should be your personal fork with "*compare: devel*". Write an appropriate title and add more descriptions to your pull request in the "*Leave a comment*" section. Finally, click the "*Create pull request*" button.

#### 6. Final merge:
We'll do the rest by reviewing your code changes, checking if the Github Actions, Travis and Azure checks all look okay. We might follow up with you by commenting on the PR, as you can still fix smaller issues in the PR by committing them to your github fork version.<br>

Finally, if the are no merge conflicts, the new version still compiles and tests pass, we'll merge your PR into the SPECFEM devel version - **with big thanks to you from the maintainers and the whole community!**


## Further informations

More detailed developer informations and instructions for how to contribute to SPECFEM codes are available at:<br>
[https://github.com/SPECFEM/specfem3d/wiki](https://github.com/SPECFEM/specfem3d/wiki).
