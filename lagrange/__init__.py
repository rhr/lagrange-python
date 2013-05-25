import os
VERSION = file(
    os.path.join(os.path.dirname(__file__), "VERSION")
    ).read().strip()

msg = """Lagrange: likelihood analysis of geographic range evolution
Version: %s
Authors: Richard Ree <rree@fieldmuseum.org>
         Stephen Smith <sasmith@nescent.org>
http://lagrange.googlecode.com
""" % VERSION

import input, output, nchoosem, optimize, ascii, newick, phylo, decmodel_mp
try:
    import graph
    #from ratemodel import RateModel, RateModelGE, RateModelGE2, Ancsplit
except:
    pass
from decmodel_mp import DECModel
from nexus import Nexus
from tree import Tree
