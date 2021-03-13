.. image:: https://travis-ci.org/hkraemer/PECUZAL_Matlab.svg?branch=main
    :target: https://travis-ci.org/hkraemer/PECUZAL_Matlab

.. image:: https://img.shields.io/badge/docs-dev-blue.svg
    :target: https://hkraemer.github.io/PECUZAL_Matlab/

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4451879.svg
   :target: https://doi.org/10.5281/zenodo.4451879

.. image:: https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg
  :target: https://de.mathworks.com/matlabcentral/fileexchange/86004-pecuzal-embedding-algorithm-for-matlab


PECUZAL Matlab
==============

We introduce the PECUZAL automatic embedding of time series method for Matlab. It is solely based
on the paper [kraemer2021]_ `(Open Source) <http://iopscience.iop.org/article/10.1088/1367-2630/abe336>`_, where the functionality is explained in detail. Here we
give an introduction to its easy usage in three examples. Enjoy Embedding!

.. image:: icon.png


Getting started
===============

There are two ways of using the proposed PECUZAL method:

- Install as Toolbox. This is the easiest way and allows the usage of the the function `pecuzal_embedding.m` independently from your current working directory. It gets treated as a built-in Matlab-function and you do not have to copy any files etc. For this simply download the 'pecuzal-embedding.mltbx' from this repository or from `Matlab-Central <https://de.mathworks.com/matlabcentral/fileexchange/86004-pecuzal-embedding-algorithm-for-matlab>`_ and double-click `pecuzal-embedding.mltbx` for installation. That's it.
- You can also download this repository and copy the folder into the MATLAB user's directory. This is usually the user's "Documents" folder appended with "MATLAB" (you can find out using the function `userpath`). Add the toolbox by the `addpath` command, e.g., `addpath ~/Documents/MATLAB/PECUZAL_Matlab` on a Linux system. For everyday use, copy this command to a `startup.m` file in the MATLAB user's directory.


NOTE
====

For performance reasons we recommend to use the implementation
in the `Julia language <https://juliadynamics.github.io/DynamicalSystems.jl/latest/embedding/unified/>`_,
in order to get fast results, especially in the multivariate case. Moreover,
it is well documented and embedded in the
`DynamicalSystems.jl <https://juliadynamics.github.io/DynamicalSystems.jl/dev/>`_ ecosystem.
The computation times can be magnitudes higher than in the Julia implementation.


Documentation and basic usage
=============================

There is a `documentation <https://hkraemer.github.io/PECUZAL_Matlab/>`_ and a
`Matlab Live-Script <https://github.com/hkraemer/PECUZAL_Matlab/blob/main/html/pecuzal_examples.mlx>`_ available including some basic usage examples.


Citing and reference
====================
If you enjoy this tool and find it valuable for your research please cite

.. [kraemer2021] Kraemer et al., "A unified and automated approach to attractor reconstruction", New Journal of Physics. `doi:10.1088/1367-2630/abe336 <https://doi.org/10.1088/1367-2630/abe336>`_ (2021).

or as BiBTeX-entry:

::

﻿   @article{Kraemer2021,
	author={Kai Hauke Kraemer and George Datseris and Jürgen Kurths and Istvan Z. Kiss and Jorge L. Ocampo-Espindola and Norbert Marwan},
	title={A unified and automated approach to attractor reconstruction},
	journal={New Journal of Physics},
	url={http://iopscience.iop.org/article/10.1088/1367-2630/abe336},
	year={2021},
	abstract={We present a fully automated method for the optimal state space reconstruction from univariate and multivariate time series. The proposed methodology generalizes the time delay embedding procedure by unifying two promising ideas in a symbiotic fashion. Using non-uniform delays allows the successful reconstruction of systems inheriting different time scales. In contrast to the established methods, the minimization of an appropriate cost function determines the embedding dimension without using a threshold parameter. Moreover, the method is capable of detecting stochastic time series and, thus, can handle noise contaminated input without adjusting parameters. The superiority of the proposed method is shown on some paradigmatic models and experimental data from chaotic chemical oscillators.}
    }


Licence
=======
This is program is free software and runs under `MIT Licence <https://opensource.org/licenses/MIT>`_.
