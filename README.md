Lagrange (Python)
=================

Lagrange is a Python package implementing likelihood models for
geographic range evolution on phylogenetic trees, with methods for
inferring rates of dispersal and local extinction and ancestral
ranges.

This software implements methods described in: Ree, R H and S A
Smith. 2008. Maximum likelihood inference of geographic range
evolution by dispersal, local extinction, and cladogenesis. Systematic
Biology 57(1):4-14.

See also: Ree, R H, B R Moore, C O Webb, and M J Donoghue. 2005. A
likelihood framework for inferring the evolution of geographic range
on phylogenetic trees. Evolution 59(11):2299-2311.

Analysis scripts can be created using the web-based configuration tool
at http://www.reelab.net/lagrange. This tool lets you configure and
download your analysis as a Python script which is then run on your
own computer.

### Quickstart
- requires `docker` and `make`
- `make build`
- `make test`
- `make run # gives shell in docker container with package installed`
