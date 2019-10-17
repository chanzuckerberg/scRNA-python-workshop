# Analysis of single cell RNA-seq data (Python)

## Pre- and post-surveys

Before the workshop begins, please fill in [this pre-survey](https://forms.gle/tMrSMyKMZQ6XKJGo6).  
At the conclusion of the workshop, please fill in [this post-survey](https://forms.gle/g5oMj8WzkwRVFyzo8).

Thanks! This will help us improve future workshops.

## Course goals

Today it is possible to obtain genome-wide transcriptome data from single cells using high-throughput sequencing (scRNA-seq). The goal of this one-day workshop is to help scientists from all backgrounds (computational or otherwise) feel empowered to explore their scRNA-seq data. Specifically, we hope students leave with the ability to:

- Understand a general workflow for dealing with scRNA-seq data
- Anticipate and avoid some of the most common pitfalls in scRNA-seq analysis
- Build intuition around the tradeoffs inherent in analytical choices
- Feel comfortable and confident working with current Python-based tools for single cell analysis
- Have better conversations with their collaborators
- Know where to find additional information and assistance

### Topics covered

We'll cover the basics of:

- Quality control
- Normalization
- Dimensionality reduction
- Clustering
- Differential expression
- Exploratory analysis

## Course non-goals

We have tried to incorporate current best practices and build intuition around methodological choices throughout this workshop. However, single cell data is complex and the field is evolving rapidly. There are also many aspects of this analysis where the field has not yet reached consensus on best practices. A one-day course simply cannot cover all of the relevant considerations and tools. Here, we have prioritized topics and tools that build foundational intuition, and are available in Python with reasonable runtime. Given the one-day timetable, we opt to start from the expression matrix and do not cover processing raw data. While not a comprehensive guide, we hope this serves as a stepping stone to making single cell analysis approachable.

## Prerequisites & resources

The workshop consists of explanatory discussions interspersed with hands-on exercises. **We strongly encourage you to bring a laptop with all required packages installed in order to fully participate.** Please follow the instructions [here](https://chanzuckerberg.github.io/scRNA-python-workshop/intro/setup.html)

The course is intended for those who have basic familiarity with Python (e.g., at the level covered in a software carpentry workshop). Basic familiarity with the Jupyter notebooks and the command line is helpful but not required.

We recommend the following introductory materials:

- **Python**: Software Carpentry workshop on ["Plotting and Programming in Python"](http://swcarpentry.github.io/python-novice-gapminder/)
- **Python**: [Codecademy Python3 course](https://www.codecademy.com/learn/learn-python-3) (free with trial).
- **Command line**: [Codecademy command line course](https://www.codecademy.com/learn/learn-the-command-line) (free with trial).
- **Jupyter notebooks**: [Coderefinery Jupyter workshop](https://coderefinery.github.io/jupyter/).

## Recommended reading

## [Analysis of single cell RNA seq data](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html)

The original source of much of the course content; as a 2-day course, it includes a broader overview of analytical methods and a guide to generating the expression matrix from raw data. While its examples are implemented in R, the conceptual underpinnings are broadly applicable.

## [Current best practices in single‐cell RNA‐seq analysis: a tutorial](https://www.embopress.org/doi/full/10.15252/msb.20188746)

Clear explanations of many of the methodological tradeoffs and our current understanding of best practices. Includes an in-depth tutorial, implemented primarily in Python.

## Development, reuse and contributing

### Content

This course was [originally developed by the Hemberg Lab](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html), and has been abbreviated and adapted for Python (with permission) by [Sidney Bell](https://twitter.com/sidneymbell) and the Computational Biology team at the Chan Zuckerberg Initiative.  
We gratefully acknowledge the work of the original authors of the course material: Vladimir Kiselev, Tallulah Andrews, Jennifer Westoby, Davis McCarthy, Maren Büttner, Jimmy Lee, Krzysztof Polanski, Sebastian Y. Müller, Elo Madissoon, Stephane Ballereau, Maria Do Nascimento Lopes Primo, Rocio Martinez Nunez and Martin Hemberg.

We have also incorporated practices and framing put forth in [Luecken and Theis, 2019](https://www.embopress.org/doi/full/10.15252/msb.20188746), and thank them for their work.

This curriculum was originally taught during one day of a CZI-sponsored workshop in Chicago, IL on October 18, 2019.

### Contributing

We warmly welcome and encourage members of the scientific community to submit updates and improvements through [github](https://github.com/chanzuckerberg/scRNA-python-workshop).

We adhere to the license of the original materials:

> All of the course material is licensed under GPL-3. Anyone is welcome to go through the material in order to learn about analysis of scRNA-seq data. If you plan to use the material for your own teaching, we would appreciate if you tell the Hemberg lab about it in addition to providing a suitable citation to the original materials.

### Contact

For conceptual questions about the original source material, please contact [Vladimir Kisilev](vladimir.yu.kiselev@gmail.com).  
For questions about the Python adaptation, please contact [Sidney Bell](twitter.com/sidneymbell).
