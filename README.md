[![PyPI](https://img.shields.io/pypi/v/biobalm?style=flat-square)](https://pypi.org/project/biobalm/) 
[![Api Docs](https://img.shields.io/badge/docs-api-yellowgreen?style=flat-square)](https://jcrozum.github.io/biobalm/)
[![Continuous integration](https://img.shields.io/github/actions/workflow/status/jcrozum/biobalm/test.yml?branch=main&style=flat-square)](https://github.com/jcrozum/biobalm/actions?query=workflow%3Atest)
[![Coverage](https://img.shields.io/codecov/c/github/jcrozum/biobalm?style=flat-square)](https://codecov.io/gh/jcrozum/biobalm) 
[![GitHub issues](https://img.shields.io/github/issues/jcrozum/biobalm?style=flat-square)](https://github.com/jcrozum/biobalm/issues) 
[![License](https://img.shields.io/pypi/l/biobalm?style=flat-square)](https://github.com/jcrozum/biobalm/blob/main/LICENSE)

# Boolean Attractor Landscape Mapper (biobalm)

biobalm is a Python library for exploring the attractor landscape of large-scale Boolean networks with hundreds or thousands of variables. It combines symbolic (BDD) and automated (ASP) reasoning to efficiently construct a *succession diagram* of a Boolean network: an inclusion-based acyclic graph of the network's trap spaces. biobalm can then use this succession diagram to accelerate attractor search and infer control strategies for target trap spaces.

### Installation

biobalm is on PyPI: **TODO: PyPI release coming soon. Use git method (below).**

```
pip install biobalm
```

The base installation allows you to generate succession diagrams and control strategies, plus some easier-to-find attractors. However, to enable the full attractor detection functionality, you need to also install `pint` and `mole`:

 - Native binaries of `pint` can be obtained [here](https://loicpauleve.name/pint/doc/#Binaries).
 - Download `mole` [here](http://www.lsv.fr/~schwoon/tools/mole/), compile it (simply run `make`), and make sure the result is in your `$PATH`.

You can also install the latest version of biobalm directly from github:

```
pip install git+https://github.com/jcrozum/biobalm.git@main
```

### Referencing biobalm

**TODO: A publication describing biobalm in detail will be available soon. Until then, please link this github repository instead.**

### Using biobalm

To learn more about how biobalm functions, you can explore the example notebooks listed below. Alternatively, biobalm's API documentation is also available [online](https://jcrozum.github.io/biobalm/).

**TODO: Usage examples coming soon.**

