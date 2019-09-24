---
title: 'Setup and installation'
prev_page:
  url: /intro/about.html
  title: 'About the course'
next_page:
  url: /intro/overview.html
  title: 'Workflow overview'
comment: "***PROGRAMMATICALLY GENERATED, DO NOT EDIT. SEE ORIGINAL FILES IN /content***"
---
# Setup

## Format

The workshop consists of explanatory discussions interspersed with hands-on exercises. We strongly encourage you to bring a laptop with all required packages installed in order to fully participate. We will also have office hours with CZI computational biology staff available in the morning to help troubleshoot any installation issues before class begins.

## Installation guide (before class)

### 1. Install Python via Anaconda

Even if you have previously installed Python, please install [Anaconda](https://docs.anaconda.com/anaconda/install/) for **Python version 3** (available on OSX, Linux, and Windows).

Anaconda is _package manager_, which means that it helps coordinate your Python installation and related _packages_ (useful code written by other people for performing specific tasks) for you so that you have a consistent _environment_ (the version of Python and the version of the code in each package that your computer looks at when doing your analysis).

You should now see an icon for "Anaconda Navigator" in your Applications folder (mac) or the Start menu (Windows). Please contact the instructors if you do not see this icon!

### 2. Install Bioconda

Not everyone who uses Anaconda is a biologist. As a result, some biology-specific packages are only available in the Bioconda _channel_ (collection of packages).

1. Open Anaconda Navigator.
2. Click on `Environments` in the left sidebar
3. Click on `Channels` in the top middle of the screen
4. Click on `Add...`
5. In the bottom text box, type `bioconda` and press `Enter`
6. Press `Update channels`

<img src="https://github.com/chanzuckerberg/scRNA-python-workshop/raw/master/content/figures/anaconda-channel.png" width=700>

### 3. Download the course config file

To handle installation of all the Python packages required for the workshop (both days, all tracks), we have prepared a configuration file that tells Anaconda how to configure your environment.

1. Download the configuration file called [`sfn-workshop.yaml` here](https://github.com/chanzuckerberg/scRNA-python-workshop/raw/master/sfn-workshop.yaml). Save it an a spot you'll remember.
2. Open Anaconda Navigator
3. Click on `Environments` in the left sidebar
4. Click on `Import` in the bottom left
5. Enter `sfn-workshop` as the Name, and browse for the `sfn-workshop.yaml` file you just downloaded
6. Press `Import`

<img src="https://github.com/chanzuckerberg/scRNA-python-workshop/raw/master/content/figures/anaconda-env.png" width=700>

### 4. Check your installation

1. Open Anaconda Navigator
2. Click on `Environments` in the left sidebar
3. Select the `sfn-workshop` environment from the list. This may take several minutes.
4. Click on the green `play` button that appears
5. Select `Open with Jupyter Notebook` from the list

<img src="https://github.com/chanzuckerberg/scRNA-python-workshop/raw/master/content/figures/anaconda-jupyter.png" width=700>

You should get a browser tab that says "Jupyter" at the top and lists all the files on your computer. This might not seem like much, but is all you need to get started! :)

## **If this does not work, please come to office hours before class so the CZI computational biology team can help you**

### 5. Download course files

We'll update this once the curriculum has been finalized; please check back soon.

### Bonus: install fun visualization tools (optional)

1. Open Anaconda Navigator
2. Click on `Environments` in the left sidebar
3. Select the `sfn-workshop` environment from the list
4. Click on the green `play` button that appears
5. Select `Open Terminal` _(it's not scary, we promise! :)_

#### For the imaging track:

6. In the terminal window that opens, copy and paste the following and press `Enter`:  
   `python3 -m pip install napari`

To test your installation of napari, type `napari` in the next terminal line and press `Enter`. You should see an application window open up.

#### For the single-cell analysis track:

6. In the terminal window that opens, copy and paste the following and press `Enter`:  
   `python3 -m pip install cellxgene`

To test your installation of cellxgene, type `cellxgene` in the next terminal line and press `Enter`. You should see something like this:  
<code>Usage: cellxgene [OPTIONS] COMMAND [ARGS]...
Options:
--version Show the version and exit.
--help Show this message and exit.
Commands:
launch Launch the cellxgene data viewer.
prepare Preprocesses data for use with cellxgene.</code>

#### If this does not work, and you would like to try out the visualization tools in class, please come to office hours before class so we can help you.