
PECUZAL Matlab
==============

We introduce the PECUZAL automatic embedding of time series method for Matlab. It is solely based
on the paper [kraemer2020]_ `(Open Source) <https://arxiv.org/abs/2011.07040>`_, where the functionality is explained in detail. Here we
give an introduction to its easy usage in three examples. Enjoy Embedding! 


Getting started
===============

The method is available as a Toolbox. Simply double-click the ... [...]. The Toolbox is also available on Matlab Central [...]

TODO NOTES
==========
- Proper citation, when accepted
- installation guide
- correct hyperrefs to DynamicalSystems.jl
- link to doc/usage html-file

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

There is a documentation available including some basic usage examples.


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
