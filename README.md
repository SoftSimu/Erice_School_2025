# Material for the Erice Summer School 

CLONE THIS SITE: [https://github.com/SoftSimu/Erice_School_2025](https://github.com/SoftSimu/Erice_School_2025)

<hr>

- Material for the Erice [7th Workshop and School on Frontiers in Water Biophysics (FWB)](https://www.waterbiophysics.eu/Main/HomePage).
- EuMINe Training School on Machine Learning in Hydrated Biosystems.
- COST Action [Eumine](https://www.eumine-cost.eu)

- Sunday, July 6, 2025, Erice, Sicily, Italy.
  - Lecture: 8:45-9:55
  - Hands-on session: 14:45-15:45
- Material by [Mikko Karttunen](https://www.softsimu.net/mikko/) and [Matt Davies](https://www.researchgate.net/profile/Matthew-Davies-48). If using the codes here, please cite the two articles below:
  - The $g_3$ three-point correlation function was originally introduced in: Sukhomlinov, S. V.; Müser, M. H. A Mixed Radial, Angular, Three-Body Distribution Function as a Tool for Local Structure Characterization: Application to Single-Component Structures. [J. Chem. Phys. 2020, 152, 194502](https://doi.org/10.1063/5.0007964).
  - It was further worked on, and integrated to an ML workflow in:  Davies, M.; Reyes-Figueroa, A. D.; Gurtovenko, A. A.; Frankel, D.; Karttunen, M. Elucidating Lipid Conformations in the Ripple Phase: Machine Learning Reveals Four Lipid Populations. [Biophys. J. 2023, 122, 442–450](https://doi.org/10.1016/j.bpj.2022.11.024).

- New to python? Here is a gentle introduction with links to more resources: [Brief introduction to Python](https://mejk.github.io/SciComp/class/PythonBasics/python-intro.html#plotting-in-python) 

<hr>

## 1. Installation instructions

The instrcutions below use the command line, and should work on Linux, MacOS terminal, and WSL. 
They have been tested on Ubuntu 24.04 LTS with python 3.12.3.

### 1.1 Check which version of python you have

```
python3 --version
```

Anything above python 3.10 should be ok.

### 1.2 Choose your preferred method

Here, we are talking command line options. There is also the Anaconda GUI, but that is not covered here. For command line, there are two main options:

1. A `conda`-based installation and
1. `pip` (Python Package Installer). This is what we will use.

Both are good choices but there are some differences. Conda installation gives immediately a very large number of packages. `pip` usually gets new packages faster and all packages must be separately installed. This also means that Jupytheer must be installed separately. 

##### 1.2.1 Check if you have `pip`

To check if you have `pip`, run the command

```
python3 -m pip --version
```
If it doesn't show anything, then you don't have `pip` installed. Try
```
python3 -m pip3 --version
```

If you don't have it, the command below is how it can be installed. The below assumes that you have superuser rights to your computer:
```
sudo apt install python3-pip
```
Now, check that it is there, and let's move on.

### 1.3 IMPORTANT: Create a virtual environment

Independent if you use `pip` or `conda`,  it is *very* important to use a virtual environment.
If you already have python virtual environment installed and are familiar with them, create a virtual environment for this exercise, activate it, and skip to Section 1.5. If you are unfamiliar with virtual environments, read on and install. 

**Why use a virtual environment?** There is a vast number of python packages both distributed with the basic installation, and contributed by the users and various projects. When packages are updated or/and changed, compatibility may (and often does) break.  The solution to this is the so-called *virtual environment*. Virtual environment provides a sandbox - an isolated environment - for a given project. This allows the project to have its own dependencies *without interfering* with other projects and, importantly, the main installation. Using virtual environments is highly recommended, and not doing it may lead to wrecking your current python installation and the package dependencies. Using a virtual environment eliminates lots of potential errors and problems with package updates.


#### 1.3.1 Using `pip`

Check if you have virtual environment installed:
```
virtualenv --version
```
If not, then, depending on the version of python, the following may work:
Install the package `virtualenv` using
```
python3 -m pip install --user virtualenv
```
If that gives an error something like this
```
error: externally-managed-environment

× This environment is externally managed
╰─> To install Python packages system-wide, try apt install
    python3-xyz, where xyz is the package you are trying to
    install.
```
Then try the following (as per the instructions above -> you need to know the version of python you have). In my case, I had this come with a new installation, and my python version is 3.12.3. With that, I did
```
sudo apt install python3.12-venv
```
and voilà!

Create the environment, the syntax is 
```
python3 -m venv my_new_env
```
replace `my_new_env` by whatever you like your environment to be called. Here, let's call the environment `g3_erice`. If you want to use that, just create the environment:
```
python3 -m venv g3_erice
```
and activate it:
```
source g3_erice/bin/activate
```
And to deactivate
```
deactivate
```

#### 1.3.2 Using `conda`

Create a new virtual environment. To use your current version of Python:
```
conda create -n my_new_environment
```
Alternatively, if you want to use some specific version of Python, instead use
```
conda create -n my_new_environment python=x.y
```
and `replace x.y` by the version number accordingly

Activate the new virtual environment
```
conda activate  my_new_environment
```
After this, the name of your new environment will appear in front of the prompt on your terminal.

Install the packages that you want for the environment called `my_new_environment`
```
conda install -n my_new_environment <package name>
```
Note: the option `-n` specifies the name of the environment. That allows for easy installation of packages into any environment. Leaving `-n` out installs in the current environment.

Deactivate the virtual environment
```
conda deactivate
```
After this, the text base will appear in front of the command prompt on your terminal.

If you want to remove your virtual environment:
```
conda env remove -n  my_new_environment
```

## 1.4 Activate your virtual environment


### 1.4.1 Install the necessary packages and dependencies

**First activate your virtual environment.** After that,  install

```
pip install jupyterlab
```
The following also installs `seaborn` for pretty plots, `scikit-learn` that is used in ML in the notebooks, `numpy`, a compatible version of `matplotlib` and many other dependecies:
```
pip install --upgrade MDAnalysis[analysis]
```

Install also the test cases for `MDAnalysis`

```
pip install --upgrade MDAnalysisTests
```
**Optional but recommended:** Run the MDAnalysis tests (more info [here](https://userguide.mdanalysis.org/stable/installation.html#testing)). This took about 7 minutes on AMD Ryzen 9 (5800 Series) based laptop running Ubuntu 24.04 LTS and aboiut 9 minutes on the same laptop using WSL. The following uses serial operations:

```
pytest --disable-pytest-warnings --pyargs MDAnalysisTests
```
There is also an option to run in parallel. Replace the number after `--numprocesses` by the number of cores you have. This took about 2 minutes om AMD Ryzen 9 (5800 Series) based laptop with 8 cores (16 CPUs).
```
pip install pytest-xdist
pytest --disable-pytest-warnings --pyargs MDAnalysisTests --numprocesses 4
```
Get `ipywidgets` for widgets and such
```
pip install ipywidgets
```
Get the `scikit-image` package
```
pip install scikit-image
```

### 1.4.2 Check if you have a given package

Using `pip`, this is how you can check if you have a given package installed

```
pip show package_name
```

## 1.5 Clone this repository

With the above, we can download the notebooks, codes, and data for the exercises. Use your preferred method under the tab `code`.

**Note 1:** If you haven't cloned GitHub repositories before, you may need to install `gh`
- [Installing the `gh` client](https://github.com/cli/cli#installation). This link has instructions for all major operating systems.

**NOTE 2:** You can also just simply download the zip file

**Make sure you have activated your virtual environment!**

## Files and directrories
Once the download is available, the main directory 

- `notebooks`: Jupyter notebooks for the exercises
- `notebooks/data`: the data files for analysis
- `notebooks/src`: pythin source file
- `notebooks/pics`: pictures used for demonstrations in the notebooks

# 2. Start Jupyter lab

Move to the notebook directory and your're ready to start.

## 2.1 How to proceed

Use the files in this order (directory `notebooks`):

1. `notebook_intro.ipynb`
2. `notebook_fcc.ipynb`
3. `notebook_spc_water.ipynb`
4. `notebook_g3_pca.ipynb`


