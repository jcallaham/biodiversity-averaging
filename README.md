## Modeling Biological Response to Climate Change

This code was the beginning of a project that didn't really bear fruit. The goal was to expand on the findings of the 2017 [paper](https://www.biorxiv.org/content/early/2017/03/08/113753) by Kristina Riemer, Rob Guralnick, and Ethan White and account for uncertainty in temperature measurement and specimen collection location.

My code integrates the NOAA climate data with the VertNet biological collection records with a weighted average described in 'weighted-averaging.pdf'. I adapted the data cleaning code to Python from the [original](https://github.com/KristinaRiemer/MassResponseToTemp), and 'analysis.py' actually carries out the weighted averaging.

In the end, as demonstrated by 'shift-hist.png', the averaging procedure doesn't make a difference in most cases, and doesn't significantly alter the original conclusions.
