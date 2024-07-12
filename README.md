# Wildfire CO<sub>2</sub> Flux Project - Summer 2024

This directory includes scripts and data for the 2024 Summer Wildfire CO<sub>2</sub> Flux Project. Please see the repository structure section below.

# Repository Structure

|Folder|Description|
|----|----|
|DATA|Collection of datasets for this project. See `info.txt` for subfolder organization information.|
|ELT CALIB|Scripts for calibrating ELT S300-3V CO<sub>2</sub> sensors. Ambient colocations are used for ambient sensors. Ambient and closed loop colocations are used for ambient sensors.|
|FLUXES|Scripts for calculating fluxes from datasets. Both static chamber method and dynamic flux chamber methods are included.|
|LICOR CALIB|Calibration files for LICOR LI-7810.|
|MFC CALIB|Calibration scripts for Mass Flow Controller (MFC).|
|MODELING|Theoretical applications of box model.|
|PCB|Hardware and firmware for data acquisition system PCB.|
|UTILS|Collection of utility functions for all scripts.|

# Contributing

To ensure smooth collaboration and avoid unnecessary merge requests, please follow these simple steps for contribution:

1) Fork the Repository
    - Click the Fork button at the top right of this repository page to create a copy of the repository under your GitHub account.
2) Clone Your Fork
    - Clone your forked repository to your local machine:
    - `git clone https://github.com/BikeLinc/co2_flux_summer_24.git`
    - `cd your-repo`
3) Create a New Branch
    - Create a new branch for your work. Name your branch based on the feature or issue you are working on:
    - `git checkout -b feature/your-feature-name`
4) Make Your Changes
    - Make your changes to the code.
5) Commit Your Changes
    - Commit your changes with a clear and concise commit message:
    - `git add .`
    - `git commit -m "Add description of changes"`
6) Push Your Branch
    - Push your branch to your forked repository:
    - `git push origin feature/your-feature-name`
7) Open a Pull Request
    - Go to the original repository on GitHub. You should see a prompt to open a pull request. Click on it and provide a clear description of your changes.
8) Wait for Review
    - Your pull request will be reviewed. Make any necessary changes based on feedback.

For an extensive guide to git and GitHub, see [GitHub's Git Guide](https://github.com/git-guides).


# Questions

Please email Lincoln (lascheer@fortlewis.edu) for any questions about this repository.
