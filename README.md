[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
[![Testing suite](https://github.com/irukoa/WannInt/actions/workflows/CI.yml/badge.svg)](https://github.com/irukoa/WannInt/actions/workflows/CI.yml)
# WannInt
### Wannier Interpolation

This is a tiny modern Fortran library to provide utilities for [Wannier interpolation](https://link.aps.org/doi/10.1103/PhysRevB.75.195121). The library is meant to serve as a building block for codes that compute the resolution of quantum mechanical operators in the Brillouin zone of a crystal.

### Working principle

The [Wannierization](https://link.aps.org/doi/10.1103/RevModPhys.84.1419) procedure consists on post-processing the results of a [DFT](https://en.wikipedia.org/wiki/Density_functional_theory) calculation into Wannier states. These states are defined on real space, and provide an obvious basis to express the reciprocal basis Bloch states. Many quantities in solid state physics (Berry curvature, optical responses, transport properties...) are expressed in reciprocal space in terms of the Hamiltonian and Berry connection operators and its derivatives. This library computes, given a matrix representation of the Hamiltonian and Berry connection in real space, $H_{nm}(\textbf{R})$ and $A^j_{nm}(\textbf{R})$ respectively, the matrix representation $H_{nm}(\textbf{k})$ and $A^j_{nm}(\textbf{k})$ and its derivatives in reciprocal space,

$$
H_{nm}^{lp\cdots }(\textbf{k}) = \sum_{R}\left[\left(iR^l\right)\left(iR^p\right)\cdots\right] e^{i\textbf{k}\cdot \textbf{R}}H_{nm}(\textbf{R}),
$$

$$
A_{nm}^{j \ lp\cdots }(\textbf{k}) = \sum_{R}\left[\left(iR^l\right)\left(iR^p\right)\cdots\right] e^{i\textbf{k}\cdot \textbf{R}}A^j_{nm}(\textbf{R}).
$$

The actual quantities $H_{nm}(\textbf{R})$ and $A^j_{nm}(\textbf{R})$ are given as inputs to the API's routines. For tight-binding models these can be given by hand. For real materials, these can be read from a file generated by the [WANNIER90](https://wannier.org/) code, by setting

```
write_tb = .true.
```
in the WANNIER90 file input card.

# Dependencies

This library uses [MAC](https://github.com/irukoa/MAC)'s containers when storing the result of the Hamiltonian and Berry connection when derivatives are requested.

# API

The derived type
``` fortran
type, public :: crystal
```
is defined.

# Build

An automated build is available for [Fortran Package Manager](https://fpm.fortran-lang.org/) users. This is the recommended way to build and use WannInt in your projects. You can add WannInt to your project dependencies by including

```
[dependencies]
MAC = { git="https://github.com/irukoa/WannInt.git" }
```
to the `fpm.toml` file.

[MAC](https://github.com/irukoa/MAC)'s objects
``` fortran
type, public :: container_specifier
type, extends(container_specifier), public :: container
```
are made public already by WannInt.
