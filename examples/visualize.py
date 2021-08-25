import networkx as nx
import matplotlib.pyplot as plt
import sys
args = sys.argv

delimeters = [";", ":", "."]
colors = ["b", "r", "g", "y", "c"]

graph_file = str(args[1])
if len(args) > 2:
  group_file = str(args[2])
else:
  group_file = ""
if len(args) > 3:
  save_file = str(args[3])
else:
  save_file = ""


# Read graph

try:
    f = open(graph_file, "r")
except IOError:
    print("Could not open file " + str(graph_file))
    exit()
G = nx.DiGraph()
#G = nx.Graph()

for x in f:
  
  for i in delimeters:
    x = x.replace(i, "")
  src = 0
  
  s = "".join(x.split())
  src = s[0]
  if(src[0] == "!") :continue
  for i in range (1, len(s)):
    G.add_edge(src, s[i])

#Read groups
try :
    f = open(group_file, "r")
    groups = []
except IOError: 
    print("No group file found")
    groups = [G.nodes()]

for x in f:
  temp = []
  s = "".join(x.split())
  for i in s:
    temp.append(i)
  groups.append(temp)

pos = nx.spring_layout(G)
for i in range (len(groups)):
  nx.draw_networkx_nodes(G, pos, nodelist=groups[i], node_color=colors[i % len(colors)])
nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.5, arrows=False)

if save_file != "":
  plt.savefig(save_file)
else :
  plt.show()
