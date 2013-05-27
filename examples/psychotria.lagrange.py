#!/usr/bin/env python
import os
import lagrange
data = """\
### begin data
{'area_adjacency': [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]],
 'area_dispersal': [[[1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0]]],
 'area_labels': ['K', 'O', 'M', 'H'],
 'base_rates': '__estimate__',
 'dispersal_durations': [10.0],
 'dm_symmetric_entry': True,
 'excluded_ranges': [],
 'lagrange_version': '20120508',
 'max_range_size': 2,
 'model_name': 'psychotria',
 'newick_trees': [{'included': '__all__',
                   'name': 'Tree0',
                   'newick': '((((((((P_hawaiiensis_WaikamoiL1:0.010853, P_mauiensis_Eke:0.010853)N2:0.007964, (P_fauriei2:0.013826, P_hathewayi_1:0.013826)N5:0.004991)N6:0.001986, (P_kaduana_PuuKukuiAS:0.020803, P_mauiensis_PepeAS:0.020803)N9:1e-05)N10:0.003762, P_kaduana_HawaiiLoa:0.024565)N12:0.003398, (P_greenwelliae07:0.012715, P_greenwelliae907:0.012715)N15:0.015248)N16:0.018984, ((((P_mariniana_MauiNui:0.02241, P_hawaiiensis_Makaopuhi:0.02241)N19:0.008236, P_mariniana_Oahu:0.030646)N21:0.002893, P_mariniana_Kokee2:0.033539)N23:0.005171, P_wawraeDL7428:0.03871)N25:0.008237)N26:0.008255, (P_grandiflora_Kal2:0.027864, P_hobdyi_Kuia:0.027864)N29:0.027338)N30:0.003229, ((P_hexandra_K1:0.026568, P_hexandra_M:0.026568)N33:0.005204, P_hexandra_Oahu:0.031771)N35:0.026659)N36;',
                   'root_age': 5.2}],
 'ranges': [(),
            (0,),
            (0, 1),
            (0, 2),
            (0, 3),
            (1,),
            (1, 2),
            (1, 3),
            (2,),
            (2, 3),
            (3,)],
 'taxa': ['P_mariniana_Kokee2',
          'P_mariniana_Oahu',
          'P_mariniana_MauiNui',
          'P_hawaiiensis_Makaopuhi',
          'P_wawraeDL7428',
          'P_kaduana_PuuKukuiAS',
          'P_mauiensis_PepeAS',
          'P_hawaiiensis_WaikamoiL1',
          'P_mauiensis_Eke',
          'P_fauriei2',
          'P_hathewayi_1',
          'P_kaduana_HawaiiLoa',
          'P_greenwelliae07',
          'P_greenwelliae907',
          'P_grandiflora_Kal2',
          'P_hobdyi_Kuia',
          'P_hexandra_K1',
          'P_hexandra_M',
          'P_hexandra_Oahu'],
 'taxon_range_data': {'P_fauriei2': (1,),
                      'P_grandiflora_Kal2': (0,),
                      'P_greenwelliae07': (0,),
                      'P_greenwelliae907': (0,),
                      'P_hathewayi_1': (1,),
                      'P_hawaiiensis_Makaopuhi': (3,),
                      'P_hawaiiensis_WaikamoiL1': (2,),
                      'P_hexandra_K1': (0,),
                      'P_hexandra_M': (0,),
                      'P_hexandra_Oahu': (1,),
                      'P_hobdyi_Kuia': (0,),
                      'P_kaduana_HawaiiLoa': (1,),
                      'P_kaduana_PuuKukuiAS': (2,),
                      'P_mariniana_Kokee2': (0,),
                      'P_mariniana_MauiNui': (2,),
                      'P_mariniana_Oahu': (1,),
                      'P_mauiensis_Eke': (2,),
                      'P_mauiensis_PepeAS': (2,),
                      'P_wawraeDL7428': (0,)}}
### end data
"""

i = 0
while 1:
    if not i:
        outfname = "psychotria.results.txt"
    else:
        outfname = "psychotria.results-"+str(i)+".txt"
    if not os.path.exists(outfname): break
    i += 1
outfile = open(outfname, "w")
lagrange.output.log(lagrange.msg, outfile, tee=True)
model, tree, data, nodelabels, base_rates = lagrange.input.eval_decmodel(data)
lagrange.output.ascii_tree(outfile, tree, model, data, tee=True)
if base_rates != "__estimate__":
    d, e = base_rates
else:
    d, e = lagrange.output.optimize_dispersal_extinction(outfile, tree, model, tee=True)
if nodelabels:
    if nodelabels == "__all__":
        nodelabels = None
    lagrange.output.ancsplits(outfile, tree, model, d, e, nodelabels=nodelabels, tee=True)
