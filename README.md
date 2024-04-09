# RobSparseMVA
A package for robust and sparse multivariate methods for multivariate high-dimensional analysis.

The package implements an algorithm for robust and sparse maximum association, ccaMM(), based on a combination of the method-of-multipliers algorithm and adaptive gradient descent, a more general variant of
canonical correlation analysis, that is described in detail in our paper:

Pfeiffer, P., Alfons, A. & Filzmoser, P. (2023+). Efficient Computation of Sparse and Robust Maximum Association Estimators. https://arxiv.org/pdf/2311.17563.pdf

In addition, the package implements an algorithm for cellwise robust and sparse PCA, pcaSCRAMBLE(), implemented via Riemannian gradient descent. If you find it helpful, please cite
Pfeiffer, P. & Filzmoser, P. (2024). Cellwise robust and sparse principal component analysis. Technical Report. 

For information on how to install the package from github, see https://cran.r-project.org/web/packages/githubinstall/vignettes/githubinstall.html.
