# Powered by Python 3.7

# To cancel the modifications performed by the script
# on the current graph, click on the undo button.

# Some useful keyboards shortcuts :
#   * Ctrl + D : comment selected lines.
#   * Ctrl + Shift + D  : uncomment selected lines.
#   * Ctrl + I : indent selected lines.
#   * Ctrl + Shift + I  : unindent selected lines.
#   * Ctrl + Return  : run script.
#   * Ctrl + F  : find selected text.
#   * Ctrl + R  : replace selected text.
#   * Ctrl + Space  : show auto-completion dialog.

from tulip import *
from collections import deque
import csv

def importData(data):
  interraction_dict ={}
  wd = '/autofs/unityaccount/cremi/ahucteau/Semestre_9/DEA/Tulip_projet/Projet_DEA/'
  file = 'interactions_chromosome6.csv'
  with open (wd+file, newline ='') as csvfile:
    data = list(csv.reader(csvfile, delimiter ='\t'))
  data.pop(0)
  lendata = len(data)
  for i in range(lendata):
    nom_interraction = data[i][1] + "_" + data[i][2]
    interraction_dict[nom_interraction] = {"locus1" : data[i][1], "locus2" : data[i][2], "interraction" : data[i][3], "distance": data[i][4]}
    i +=1
  return interraction_dict

def create_nodes_edges(gr,dico,nodes):
  for key in dico.keys():
    locus=dico[key]["locus1"]
    if locus not in nodes.keys():
      nodes[locus]=gr.addNode()
    locus2=dico[key]["locus2"]
    if locus2 not in nodes.keys():
      nodes[locus2]=gr.addNode()
    dico[key]['edge']=gr.addEdge(nodes[locus],nodes[locus2])


def construireGraph(gr,locus1,locus2,distances,edges,nodes):
  for n in range(len(locus1)):
    gr.addNode()


# The updateVisualization(centerViews = True) function can be called
# during script execution to update the opened views

# The pauseScript() function can be called to pause the script execution.
# To resume the script execution, you will have to click on the "Run script " button.

# The runGraphScript(scriptFile, graph) function can be called to launch
# another edited script on a tlp.Graph object.
# The scriptFile parameter defines the script name to call (in the form [a-zA-Z0-9_]+.py)

# The main(graph) function must be defined
# to run the script on the current graph

def main(graph):
  data = []
  locus1 = []
  locus2 = []
  interaction = []
  distances = []
  nodes = {}
  edges = {}
  interract_dict = {}

  interraction_dict=importData(data)
  create_nodes_edges(graph, interraction_dict, nodes)

  viewBorderColor = graph.getColorProperty("viewBorderColor")
  viewBorderWidth = graph.getDoubleProperty("viewBorderWidth")
  viewColor = graph.getColorProperty("viewColor")
  viewFont = graph.getStringProperty("viewFont")
  viewFontSize = graph.getIntegerProperty("viewFontSize")
  viewIcon = graph.getStringProperty("viewIcon")
  viewLabel = graph.getStringProperty("viewLabel")
  viewLabelBorderColor = graph.getColorProperty("viewLabelBorderColor")
  viewLabelBorderWidth = graph.getDoubleProperty("viewLabelBorderWidth")
  viewLabelColor = graph.getColorProperty("viewLabelColor")
  viewLabelPosition = graph.getIntegerProperty("viewLabelPosition")
  viewLayout = graph.getLayoutProperty("viewLayout")
  viewMetric = graph.getDoubleProperty("viewMetric")
  viewRotation = graph.getDoubleProperty("viewRotation")
  viewSelection = graph.getBooleanProperty("viewSelection")
  viewShape = graph.getIntegerProperty("viewShape")
  viewSize = graph.getSizeProperty("viewSize")
  viewSrcAnchorShape = graph.getIntegerProperty("viewSrcAnchorShape")
  viewSrcAnchorSize = graph.getSizeProperty("viewSrcAnchorSize")
  viewTexture = graph.getStringProperty("viewTexture")
  viewTgtAnchorShape = graph.getIntegerProperty("viewTgtAnchorShape")
  viewTgtAnchorSize = graph.getSizeProperty("viewTgtAnchorSize")

  for n in graph.getNodes():
    print(n)
