# DrSwine

> Minimal LePIGS implementation in Racket.

DrSwine is PIGS, MD (Path Integral Ground State Molecular Dynamics).
Specifically, it's an implementation of [LePIGS](http://pubs.acs.org/doi/abs/10.1021/jp4015178) (Langevin equation Path Integral Ground State) in Racket with not much else.
It is especially lacking in the speed department.
(If you're looking for a serious MD package for path integrals, consider [MMTK](https://bitbucket.org/khinsen/mmtk/branch/path_integrals).)

The aim of this package is to provide the user with the ability to quickly set up simulations with different path configurations, without the restrictions that may normally be applied when simulating physical systems.
For example, the paths do not correspond to any "atoms" or "particles", and the simulation temperature may be regulated independently of the path lengths.


## Dependencies

* The Fourier transforms are done using the [Science Collection](http://planet.racket-lang.org/package-source/williams/science.plt/4/8/planet-docs/science/index.html), which is available under the LGPLv3.


## License

Provided under the terms of the MIT license.
See LICENSE.txt for more information.
