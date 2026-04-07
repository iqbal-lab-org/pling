from ete3 import Tree, TreeStyle, NodeStyle
import argparse
import pandas as pd
import numpy as np
import math

parser = argparse.ArgumentParser()
parser.add_argument('input')
parser.add_argument('output')
args = parser.parse_args()

filepath = args.input
t = Tree(filepath, format=0, quoted_node_names=True)

ts = TreeStyle()
ts.show_leaf_name = True
ts.scale = 250
ts.branch_vertical_margin = 5

for n in t.traverse():
	if n.is_leaf():
		nstyle = NodeStyle()
		nstyle['shape'] = 'circle'
		nstyle['size'] = 4
		nstyle['fgcolor'] = '#ff0000'
		nstyle['hz_line_width'] = 2
		nstyle['vt_line_width'] = 2
		n.set_style(nstyle)
	else:
		nstyle=NodeStyle()
		nstyle['size']=0
		nstyle['hz_line_width'] = 2
		nstyle['vt_line_width'] = 2
		n.set_style(nstyle)

t.render(str(args.output), tree_style=ts)
