# Executable Environment for OSF Project [xjm6z](https://osf.io/xjm6z/)

This repository was automatically generated as part of a project to test the reproducibility of open science projects hosted on the Open Science Framework (OSF).

**Project Title:** Bayesian analysis of cross-sectional network psychometrics: A tutorial in R and JASP

**Project Description:**
> Network psychometrics is a new direction in psychological research that conceptualizes psychological constructs as systems of interacting variables. In network analysis, variables are represented as nodes and their interactions yield (partial) associations. Current estimation methods mostly use a frequentist approach, which does not allow for proper uncertainty quantification of the model and its parameters. Here, we outline a Bayesian approach to network analysis that offers three main benefits. In particular, applied researchers can use Bayesian methods to (1) determine structure uncertainty, (2) obtain evidence for edge inclusion and exclusion (i.e., distinguish conditional (in)dependence between variables), and (3) quantify parameter precision. The paper provides a conceptual introduction to Bayesian inference, describes how researchers can facilitate the three benefits for networks, and reviews the available R packages. In addition, we present two user-friendly software solutions: a new R package easybgm for fitting, extracting, and visualizing the Bayesian analysis of networks, and a graphical user interface implementation in JASP. The methodology is illustrated with a worked-out example of a network of personality traits and mental health. 

**Original OSF Page:** [https://osf.io/xjm6z/](https://osf.io/xjm6z/)

---

**Important Note:** The contents of the `xjm6z_src` folder were cloned from the OSF project on **12-03-2025**. Any changes made to the original OSF project after this date will not be reflected in this repository.

The `DESCRIPTION` file was automatically added to make this project Binder-ready. For more information on how R-based OSF projects are containerized, please refer to the `osf-to-binder` GitHub repository: [https://github.com/Code-Inspect/osf-to-binder](https://github.com/Code-Inspect/osf-to-binder)

## flowR Integration

This version of the repository has the **[flowR Addin](https://github.com/flowr-analysis/rstudio-addin-flowr)** preinstalled. flowR allows visual design and execution of data analysis workflows within RStudio, supporting better reproducibility and modular analysis pipelines.

To use flowR, open the project in RStudio and go to `Addins` > `flowR`.

## How to Launch:

**Launch in your Browser:**

🚀 **MyBinder:** [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/code-inspect-binder/osf_xjm6z-f/HEAD?urlpath=rstudio)

   * This will launch the project in an interactive RStudio environment in your web browser.
   * Please note that Binder may take a few minutes to build the environment.

🚀 **NFDI JupyterHub:** [![NFDI](https://nfdi-jupyter.de/images/nfdi_badge.svg)](https://hub.nfdi-jupyter.de/r2d/gh/code-inspect-binder/osf_xjm6z-f/HEAD?urlpath=rstudio)

   * This will launch the project in an interactive RStudio environment on the NFDI JupyterHub platform.

**Access Downloaded Data:**
The downloaded data from the OSF project is located in the `xjm6z_src` folder.

## Run via Docker for Long-Term Reproducibility

In addition to launching this project using Binder or NFDI JupyterHub, you can reproduce the environment locally using Docker. This is especially useful for long-term access, offline use, or high-performance computing environments.

### Pull the Docker Image

```bash
docker pull meet261/repo2docker-xjm6z-f:latest
```

### Launch RStudio Server

Run the container (with a name, e.g. `rstudio-dev`):
```bash
docker run -it --name rstudio-dev --platform linux/amd64 -p 8888:8787 --user root meet261/repo2docker-xjm6z-f bash
```

Inside the container, start RStudio Server with no authentication:
```bash
/usr/lib/rstudio-server/bin/rserver --www-port 8787 --auth-none=1
```

Then, open your browser and go to: [http://localhost:8888](http://localhost:8888)

> **Note:** If you're running the container on a remote server (e.g., via SSH), replace `localhost` with your server's IP address.
> For example: `http://<your-server-ip>:8888`

## Looking for the Base Version?

For the original Binder-ready repository **without flowR**, visit:
[osf_xjm6z](https://github.com/code-inspect-binder/osf_xjm6z)

