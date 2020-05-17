# ICTS 2020, Numerical Hydrodynamics

Material for the (remote) [ICTS workshop on Gravitation Wave Astrophysics](https://www.icts.res.in/program/gws2020).

This repository contains handouts, reading lists, exercises and code examples for the *Numerical Hydrodynamics* section of the workshop.

The slides for the lectures can be [found online](http://ianhawke.github.io/slides).

## Use of coding exercises

Full solutions are provided in Python, as well as stub codes (to remove some of the tedium of setting up boilerplate code and simplifying plotting). To use this code you'll need to install the libraries

* `numpy`
* `scipy`
* `matplotlib`
* `numba` (useful, not essential)

I recommend the [Anaconda distribution](https://www.anaconda.com/products/individual) if you don't already have this installed.

### Stub codes

In the `coding_exercises` directory there are the `solutions` (mostly complete, but could be extended) and the `stubs`. The bulk of the code is there: all class and function definitions are complete, as is all the plotting code. Look for `#! To be completed` to find the required lines. Much that needs completing is in the `systems` directory, containing the classes for each separate system, but some lines are missing from the methods code at the top level as well.
