# CARP: The Clustering Algorithms Referee Package

### Authors: Volodymyr Melnykov and Ranjan Maitra

### Package: CARP (version 1.0)

#### Description

The C-package CARP is a convenient and easy tool for evaluating performance of clustering algorithms. The underlying methodology is based on Gaussian mixture model generation according to prespecified levels of average and maximum pairwise overlaps. The concept of overlap is defined as the sum of two misclassification probabilities Maitra and Melnykov (2010). Upon simulating Gaussian mixtures, datasets are simulated from them. This concludes the first step of the procedure. At the second step, a clustering algorithm which is being tested has to be run on simulated datasets. As an example, the hierarchical agglomerative algorithm (“hierclust”) is included. At the final step, original true classification vectors are compared with estimated classifications by means of the Adjusted Rand index. A user can also specify some other program comparing two partitions. CARP is released under the GNU GPL license.

If you use this package, please cite:

Melnykov, V. & R. Maitra. (2011) "CARP: Software for fishing out good clustering algorithms." *Journal of Machine Learning Research*, 12:31-35. URL: <https://jmlr.csail.mit.edu/papers/v12/melnykov11a.html>.

Maitra, R. & V. Melnykov. (2010) "Simulating data to study performance of finite mixture modeling and model-based clustering algorithms." *Journal of Computational and Graphical Statistics*, 19(2):354-376. DOI: 10.1198/jcgs.2009.08054. <https://www.tandfonline.com/doi/abs/10.1198/jcgs.2009.08054>

#### Installation

After download, compile using at the prompt:

`make`

#### Usage

'./CARP <set of parameters>'

''''
Parameters:
-b: average overlap (no default value)
-m: maximum overlap (no default value)
-p: number of dimensions (2 by default)
-K : number of mixing components (2 by default)
-n: number of observations generated from every mixture (0 by default)
-#: number of simulated mixtures (1 by default)
-0 : name of clustering program (no default name)
-1 : name of partition analyzing program (“AdjRand” by default)
-s: spherical covariance matrix structure (non-spherical by default if option is unspecified)
-e: maximum eccentricity (0.90 by default)
-z : smallest mixing proportion (equal mixing proportions 1/K by default)
-u: upper bound for Uniform(0, upper-bound ) distribution from which mean vectors are generated
-r : maximum number of resimulations (100 by default)
-a: accuracy of estimation (1e-06 by default)
-l : maximum number of integration terms (1e06 by default)
-P : name of the file containing mixing proportions (“Pi.dat” by default)
-M : name of the file containing mean vectors (“Mu.dat” by default)
-S : name of the file containing covariance matrices in triangular form (“S.dat” by default)
-D: name of working directory (“DATA” by default)
-I : name of the file containing numbers of observations generated from every cluster (“Nk.dat” by default)
-i : name of the file containing estimated classifications (“idEst.dat” by default)
-X : name of the file containing simulated datasets (“x.dat” by default)
-W : name of the file containing maps of pairwise overlaps (“overMap.dat” by default)
-C : name of the file containing characteristics of simulated mixtures (“overBarMax.dat” by default)
-R: name of the file containing index values (Adjusted Rand index and “AR.dat” by default)
''''

#### Details

Upon launching CARP, three stages of the program have to be fulfilled. The first stage is responsible for the simulation of Gaussian mixtures with pre-specified level of complexity expressed in terms of pairwise overlap and generating datasets from these mixtures. For the second stage, a user-specified clustering method has to be run for the simulated datasets. For the illustration purposes, the hierarchical clustering algorithm `hierclust` is employed. At the third stage, the original classification is compared with the estimated one to assess the performance of the clustering algorithm run at the second stage. The Adjusted Rand index is incorporated as a default measure of similarity between the true and estimated classifications. A user, however, can specify other programs for comparing the obtained and true partitions.

If both options `-b` and `-m` are specified, CARP produces a mixture satisfying both characteristics. If one option `-b` or `-m` is specified, a mixture satisfying the prespecified value is generated. The obtained parameters will be saved to the files specified by options `-P -M -S`. The working directory is specified by the option `-D`. In addition to these, the map of misclassification probabilities will be saved into the file specified by the option `-W` while average and maximum overlaps as well as the row and column numbers of the components that produced maximum overlap are stored in the file given by the option `-C`. The element with the index `(i, j)` in the misclassification map represents the probability that `X` simulated from the `i`th component is classified to the `j`th component. If both options `-b` and `-m` are not specified, CARP reads parameters from files specified by options `-P -M -S` and computes misclassification probabilities for the components of the given mixture. Note that options `-p` and `-K` have to be appropriately specified. For existing datasets, this feature can be used for assessing the level of clustering complexity based on the set of estimated parameters.

Options `-n` and `-#` specify the sample size of a dataset simulated from every generated mixture and the number of such mixtures correspondingly. If it is desired to use datasets stored in a file specified by the option `-X`, the sample size should be given as a negative number, for example: `-n-100` instead of the more rational `-n100`. Then, new datasets will not be simulated. Note that this capability is available only in the mode when mixture parameters are not simulated but also read from files.

#### Examples

% simulate a 3-dimensional 4-component mixture with spherical covariance matrices, equal mixing proportions, maximum overlap 0.1 and average overlap 0.05; generate a sample of size 100 and analyze using "hierclust"

> ./CARP -m0.1 -b0.05 -p3 -K4 -s -n100 -0hierclust

% simulate a 2-dimensional 4-component mixture with maximum overlap 0.1 and mixing proportions greater or equal than 0.10; generate a sample of size 100
> ./CARP -m0.1 -p2 -K4 -z0.10 -n100

% compute the overlap for the 2-dimensional 4-component mixture with parameters specified in "DATA/Pi.dat", "DATA/Mu.dat" and "DATA/LTSigma.dat"; generate a sample of size 100
> ./CARP -p2 -K4 -n100

% get help
> ./CARP

##### Note

There also exists a R package on CRAN titled `MixSim`

##### References

Melnykov, V. & R. Maitra. (2011) "CARP: Software for fishing out good clustering algorithms." *Journal of Machine Learning Research*, 12:31-35. URL: <https://jmlr.csail.mit.edu/papers/v12/melnykov11a.html>.

Maitra, R. & V. Melnykov. (2010) "Simulating data to study performance of finite mixture modeling and model-based clustering algorithms." *Journal of Computational and Graphical Statistics*, 19(2):354-376. DOI: 10.1198/jcgs.2009.08054. <https://www.tandfonline.com/doi/abs/10.1198/jcgs.2009.08054>
