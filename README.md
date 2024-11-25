<p align="center">
  <img src="figures/img/DNB-BINF-Core-logo.png" alt="ScATACSeq repository logo" width="560px" />
</p>
<p align="center">
  <a href="https://www.repostatus.org/#active"><img src="https://www.repostatus.org/badges/latest/active.svg?style=for-the-badge" alt="The project has reached a stable, usable state and is being actively developed." /></a>

</p>

# Single cell/Single-nucleus ATAC Seq workflow (ScATACSeq)

This repository contains a collection of analysis modules designed to process and analyze single cell ATAC (scATAC) data from 10X sequencing technology.  


## To access the code in this repository:

### Clone option
1. Clone the repository
```
git clone https://github.com/stjude-dnb-binfcore/sc-atac-seq.git
```

### Fork option
1. Fork the repository on your own account from the main page of the `stjude-dnb-binfcore/sc-atac-seq` by clicking the “Fork” button

<img width="650" alt="how-to-fork-repo-1" src="https://github.com/user-attachments/assets/1fc0a459-2c8c-4d2e-ab6b-6abaafae963e">


2. Change the name if you like, but probably not; click “Create fork”

<img width="650" alt="how-to-fork-repo-2" src="https://github.com/user-attachments/assets/914a3db5-6e87-41fb-baf2-a50ffdb2a7c0">


3. Enjoy your new project repo!

<img width="650" alt="how-to-fork-repo-3" src="https://github.com/user-attachments/assets/073abb78-3993-4527-a574-859fd3046d39">


## To run the code in this repository:

1. Replace the `project_parameters.Config.yaml` with your file paths and parameters.

2. Navigate to an analysis module and run the shell script of interest:
```
cd ./sc-atac-seq/analyses/<module_of_interest>
```

**If you have forked the repo, you will need to do the following steps before running the script of interest. If you have cloned the repo, you can skip this.**

You need to do `sync fork` of your project repo at GitHub before running a module, if your branch is behind the main branch of the `stjude-dnb-binfcore/sc-atac-seq:main`. This will update the main branch of your project repo with the new code and modules (if any). This will add code and not break any analyses already run in your project repo. 

Then navigate to your `./sc-atac-seq` project repo and ensure you are at the `main` branch (in case you have featured branches). If not, you will need to `git checkout` to the main branch.
```
git branch
git checkout main
```

Finally, `git pull` to get the most updated changes and code in your project repo. Please be mindful of any changes in files in your project repo that you have done, e.g., `project_parameters.Config.yaml`. You will need to commit or stash (or restore) the changes to the yaml before completing the pull.
```
git pull
```


### Below is the main directory structure listing the analyses and data files used in this repository

```
├── analyses
|  ├── ...
|  └── README.md
├── figures
├── LICENSE
├── project_parameters.Config.yaml
├── README.md
└── SECURITY.md
```

## Contact

Contributions, issues, and feature requests are welcome! Please feel free to check [issues](https://github.com/stjude-dnb-binfcore/sc-atac-seq/issues).

---

*These tools and pipelines have been developed by the Bioinformatic core team at the [St. Jude Children's Research Hospital](https://www.stjude.org/). These are open access materials distributed under the terms of the [BSD 2-Clause License](https://opensource.org/license/bsd-2-clause), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
