# Exploratory analysis with cellxgene

## Introduction to cellxgene

Cellxgene (pronounced "cell-by-gene") is an interactive data explorer for single-cell transcriptomics datasets. You can see a working example [here](https://www.kidneycellatlas.org/fetal-kidney-immune.)

As input, it takes an `h5ad` file (AnnData) that contains a pre-computed low-dimensional embedding (e.g., tSNE or UMAP) and, optionally, additional metadata (e.g., annotated tissue types or precomputed cluster labels).

Let's load the brain data we've been working with into cellxgene, and then spend some time exploring the dataset in a bit more detail.

## Launching cellxgene

0. If you haven't done so already, follow the instructions at the bottom of [this page](https://chanzuckerberg.github.io/scRNA-python-workshop/intro/setup.html) to install cellxgene and verify your installation
1. Open a terminal window (Mac: `Applications > Utilities > Console`; Windows: `Start > All Programs > Accessories > Command Prompt`
1. In the terminal, type `cellxgene launch path/to/course/folder/content/data/brain_clusters.h5ad --open`  
   You need to replace `path/to/course/folder` with the correct file path to where your course folder lives on your computer.  
   For example, for a mac user with the course folder on the desktop, I would enter `Desktop/scRNA-python-workshop/content/data/brain_clusters.h5ad`

## Exercises: exploring your data

Cellxgene can be used for many different exploratory tasks; here are a few examples to get you started.

### Looking for confounders

Try coloring by the `mouse.id`; does it look like the structure in the embedded data is determined by which mouse a cell came from, or are they intermingled? What does that tell you?  
Try this for each other possible confounder in the metadata. What do you observe?

### QC cluster assignments

Color by `mouse.sex` and expand `louvain` (cluster labels). Are the clusters we identified made up of all cells from one mouse, or are they intermingled? What does that tell you?

### Visualizing gene expression

Add your favorite gene and color by its expression. Because most genes are not expressed in most cells, you may need to use the `clip` function (menu in the top right) to see the distribution.
Now look at its distribution across cell types (expand the `cell_ontology_class` metadata). Are there any interesting differences?

### Compare groups of cells

Select an interesting group of cells and assign it to `Set1`. Repeat for another group of cells and save them in `Set2`. Click `Compute differential expression`. What do you make of the returned differentially expressed genes? Are the p-values meaningful? What are some of the caveats for interpreting these p-values? How about the log-fold change?

### Assign cell types

Practice identifying a population of cells based on their marker genes, and save your observations by creating a `New label`.

## Resources

For more information about cellxgene, [read the docs here](https://chanzuckerberg.github.io/cellxgene/).

To request a feature or report a bug, [check out the github repo](https://github.com/chanzuckerberg/cellxgene/issues).

To chat with other users or contact the developers with questions, comments, or suggestions, [join the `#cellxgene-users` slack channel](https://join-cellxgene-users.herokuapp.com/).
