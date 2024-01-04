[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7736212.svg)](https://doi.org/10.5281/zenodo.7736212)

# Discrepancy measures for global sensitivity analysis

Arnald Puy, Pamphile Roy and Andrea Saltelli

This is the ``R`` code of the report, whose abstract reads as follows:

*While sensitivity analysis improves the transparency and reliability of mathematical models, its uptake by modelers is still scarce. This is partially explained by its technical requirements, which may be hard to understand and implement by the non-specialist. Here we propose a sensitivity analysis approach based on the concept of discrepancy that is as easy to understand as the visual inspection of input-output scatterplots. Firstly, we show that some discrepancy measures are able to rank the most influential parameters of a model almost as accurately as the variance-based total sensitivity index. We then introduce an ersatz-discrepancy whose performance as a sensitivity measure matches that of the best-performing discrepancy algorithms, is orders of magnitude faster, simple to implement and much easier to interpret.*

## The code

As indicated in the R code, you first need to install the ``sensobol`` package, which includes the function ``load_packages``. This function allows to install and load all the required packages for the analysis in one go. Also note that the code makes use of several functions in ``C++`` (``.cpp``), which were developed for the ``R`` package ``sensitivity`` (Iooss et al. 2020). In order to run them you need to have installed a C++ compiler such as Clang, version 15.0.0.

Some sections of the code involve computationally intensive operations and hence the simulation may take a while to complete.

Hereby we inform on which code snippet reproduces which figure in the manuscript:
* Figure 1: "plot_list_functions"
* Figure 2: "plot_sampling"
* Figure 3: No code snippet as it is a tikzdevice figure.
* Figure 4: "scatter_mood"
* Figure 5: "plot_scatter"
* Figure 6: "plot_time_complexity"

## References

Iooss, B., Janon, A., Pujol, G., with contributions from Baptiste Broto, Boumhaout, K., Veiga, S. D., Delage, T., Amri, R. E., Fruth, J., Gilquin, L., Guillaume, J., Le Gratiet, L., Lemaitre, P., Marrel, A., Meynaoui, A., Nelson, B. L., Monari, F., Oomen, R., Rakovec, O., Ramos, B., Roustant, O., Song, E., Staum, J., Sueur, R., Touati, T. & Weber, F. 2020. *sensitivity: Global Sensitivity Analysis of Model Outputs*. R package version 1.27.0.

Puy, A., Lo Piano, S., Saltelli, A., Levin, S.A., 2022. sensobol: an R package to compute variance-based sensitivity indices. *Journal of Statistical Software* 102, 1–37. https://doi.org/10.18637/jss.v102.i05


