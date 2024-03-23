[![PyPI](https://img.shields.io/pypi/v/balm?style=flat-square)](https://pypi.org/project/balm/) 
[![Api Docs](https://img.shields.io/badge/docs-api-yellowgreen?style=flat-square)](https://todo.docs.url)
[![Continuous integration](https://img.shields.io/github/actions/workflow/status/jcrozum/balm/test.yml?branch=master&style=flat-square)](https://github.com/jcrozum/balm/actions?query=workflow%3Atest)
[![Coverage](https://img.shields.io/codecov/c/github/jcrozum/balm?style=flat-square)](https://codecov.io/gh/jcrozum/balm) 
[![GitHub issues](https://img.shields.io/github/issues/jcrozum/balm?style=flat-square)](https://github.com/jcrozum/balm/issues) 
[![License](https://img.shields.io/pypi/l/balm?style=flat-square)](https://github.com/jcrozum/balm/blob/main/LICENSE)

# Boolean Attractor Landscape Mapper (BALM)

BALM is a Python library for exploring the attractor landscape of large-scale Boolean networks with hundreds or thousands of variables. It combines symbolic (BDD) and automated (ASP) reasoning to efficiently construct a *succession diagram* of a Boolean network: an inclusion-based acyclic graph of the network's trap spaces. BALM can then use this succession diagram to accelerate attractor search and infer control strategies for target trap spaces.

### Installation

BALM is on PyPI: **TODO: PyPI release coming soon. Use git method (below).**

```
pip install balm
```

The base installation allows you to generate succession diagrams and control strategies, plus some easier-to-find attractors. However, to enable the full attractor detection functionality, you need to also install `pint` and `mole`:

 - Native binaries of `pint` can be obtained [here](https://loicpauleve.name/pint/doc/#Binaries).
 - Download `mole` [here](http://www.lsv.fr/~schwoon/tools/mole/), compile it (simply run `make`), and make sure the result is in your `$PATH`.

You can also install the latest version of BALM directly from github:

```
pip install git+https://github.com/jcrozum/balm.git@main
```

### Referencing BALM

**TODO: A publication describing BALM in detail will be available soon. Until then, please link this github repository instead.**

### Using BALM

To learn more about how BALM functions, you can explore the example notebooks listed below. Alternatively, BALM's API documentation is also available [online](https://todo.docs.url).

**TODO: Usage examples coming soon.**

