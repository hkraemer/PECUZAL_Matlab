.. image:: https://travis-ci.org/hkraemer/PECUZAL_Matlab.svg?branch=main
    :target: https://travis-ci.org/hkraemer/PECUZAL_Matlab

.. image:: https://img.shields.io/badge/docs-dev-blue.svg
    :target: https://hkraemer.github.io/PECUZAL_Matlab/


PECUZAL Matlab
==============

We introduce the PECUZAL automatic embedding of time series method for Matlab. It is solely based
on the paper [kraemer2020]_ `(Open Source) <https://arxiv.org/abs/2011.07040>`_, where the functionality is explained in detail. Here we
give an introduction to its easy usage in three examples. Enjoy Embedding! 


Getting started
===============

There are two ways of using the proposed PECUZAL method:
- Install as Toolbox. This is the easiest way and allows the usage of the the function `pecuzal_embedding.m` independently from your current working directory.
It gets treated as a built-in Matlab-function and you do not have to copy any files etc. For this simply download the 'pecuzal-embedding.mltbx' from this repository 
or from Matlab-Central (hyperref HERE) and double-click for installation. That's it.
PIC-HERE
- You can also download this repository and copy all functions contained in the `/src`-folder into the working directory, in which you'd like to use the function `pecuzal_embedding.m`.


TODO NOTES
==========
- Proper citation, when accepted
- Finish installation guide
- correct hyperrefs to DynamicalSystems.jl

NOTE
====

For performance reasons we recommend to use the implementation
in the `Julia language <https://juliadynamics.github.io/DynamicalSystems.jl/dev/>`_,
in order to get fast results, especially in the multivariate case. Moreover,
it is well documented and embedded in the 
`DynamicalSystems.jl <https://juliadynamics.github.io/DynamicalSystems.jl/dev/>`_ ecosystem.
The computation times can be magnitudes higher than in the Julia implementation.


Documentation and basic usage
=============================

There is a `documentation <https://hkraemer.github.io/PECUZAL_Matlab/>`_ and a 
`Matlab Live-Script <https://github.com/hkraemer/PECUZAL_Matlab/blob/main/docs/pecuzal_examples.mlx>`_ available including some basic usage examples.


Citing and reference
====================
If you enjoy this tool and find it valuable for your research please cite

.. [kraemer2020] Kraemer et al., "A unified and automated approach to attractor reconstruction",  `arXiv:2011.07040 [physics.data-an] <https://arxiv.org/abs/2011.07040>`_, 2020.

or as BiBTeX-entry:

::

    @misc{kraemer2020,
    title={A unified and automated approach to attractor reconstruction}, 
    author={K. H. Kraemer and G. Datseris and J. Kurths and I. Z. Kiss and J. L. Ocampo-Espindola and N. Marwan},
    year={2020},
    eprint={2011.07040},
    archivePrefix={arXiv},
    primaryClass={physics.data-an}
    url={https://arxiv.org/abs/2011.07040}
    }


Licence
=======
This is program is free software and runs under `MIT Licence <https://opensource.org/licenses/MIT>`_.
