import os
VERSION = open(
    os.path.join(os.path.dirname(__file__), "VERSION")
    ).read().strip()

msg = """Lagrange: likelihood analysis of geographic range evolution
Version: %s
Author: Richard Ree <rree@fieldmuseum.org>
https://github.com/rhr/lagrange-python
""" % VERSION

try:
    from . import input, output, nchoosem, optimize, ascii, newick, phylo, decmodel_mp
except ImportError:
    import input, output, nchoosem, optimize, ascii, newick, phylo, decmodel_mp
from .decmodel_mp import DECModel
from .nexus import Nexus
from .tree import Tree

try:
    from . import graph
    #from ratemodel import RateModel, RateModelGE, RateModelGE2, Ancsplit
except ImportError:
    try:
        import graph
    except ImportError:
        pass
