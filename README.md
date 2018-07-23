# HUMICHIP PROJECT

## Project Overview

## Directory Structure

```
.
├── README.md		# This README
├── code		    # Source code (R, Python, Bash, etc)
├── makefile		# Source code (R, Python, Bash, etc)
├── data
│   ├── processed	# Data that has been cleaned and is ready to be analyzed
│   └── raw		    # Raw unmanipulated data. *DO NOT MODIFY*
├── exploratory		
│   ├── notebook	# Subdirectories for exploratory analyses
│   └── scratch		# Temporary files that can be safely deleted
├── docs		    # Collection of files (powerpoint, word, excel, etc) used for
│			          reference
└── results
    ├── figures	    # Figures, likely generated from R scripts 
    │                 (all figs must have make command)
    ├── pictures	# Pictures, likely to be incorporated into Rmarkdown reports
    ├── report.Rmd	# Report of findings, to be sent to collaborators as status updates
    ├── table		# Text tables to be rendered with kable in R 
    │                 (all figs must have make command)
```
