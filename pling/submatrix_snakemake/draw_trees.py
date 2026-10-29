import os
#render without a display (no cocoa/xcb needed), must be set before Qt is imported
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from ete4 import Tree
from ete4.treeview import TreeStyle, NodeStyle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('input')
parser.add_argument('output')
args = parser.parse_args()

filepath = args.input
#quicktree splits the newick over multiple lines, which ete4's parser doesn't accept
with open(filepath) as f:
	t = Tree(f.read().replace("\n", "").replace("\r", ""), parser=0)

ts = TreeStyle()
ts.show_leaf_name = True
ts.scale = 250
ts.branch_vertical_margin = 5

for n in t.traverse():
	if n.is_leaf:
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
