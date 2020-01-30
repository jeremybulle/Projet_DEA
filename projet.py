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

def importData():
  interraction_dict ={}
  wd = '/net/stockage/TulipDEA/'
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

def import_chromosome6_fragments_expressions(nodes, expression):
  wd = '/net/stockage/TulipDEA/'
  file = 'chromosome6_fragments_expressions.csv'
  with open(wd+file, newline='') as csvfile:
    read=list(csv.reader(csvfile, delimiter='\t'))
  read.pop(0)
  lenread=len(read)
  for ligne in range(lenread):
    if read[ligne][1] in nodes.keys():
      nodes[read[ligne][1]]["expression"]=read[ligne][2]

def import_metabos(nodes, edges):
  wd = '/net/stockage/TulipDEA/'
  file = 'KEGG.symbols.csv'
  with open(wd+file, newline='') as csvfile:
    read=list(csv.reader(csvfile, delimiter='\t'))
  for ligne in range(len(read)):
    for index_gene in range(2,len(read[ligne])):
      if read[ligne][index_gene] in nodes.keys():
        nodes[read[ligne][index_gene]]["metabo"].append(read[ligne][0])
  file = 'REACTOME.symbols.csv'
  with open(wd+file, newline='') as csvfile:
    read=list(csv.reader(csvfile, delimiter='\t'))
  for ligne in range(len(read)):
    for index_gene in range(2,len(read[ligne])):
      if read[ligne][index_gene] in nodes.keys():
        nodes[read[ligne][index_gene]]["reactome"].append(read[ligne][0])

def create_nodes_edges(gr,dico,nodes):
  for key in dico.keys():
    locus=dico[key]["locus1"]
    if locus not in nodes.keys():
      nodes[locus]={}
      nodes[locus]={"node":gr.addNode(), "metabo":[], "reactome":[]}
    locus2=dico[key]["locus2"]
    if locus2 not in nodes.keys():
      nodes[locus2]={}
      nodes[locus2]={"node":gr.addNode(), "metabo":[], "reactome":[]}
    dico[key]["edge"]=gr.addEdge(nodes[locus]["node"],nodes[locus2]["node"])

def ajouter_metrics(nodes, metrics, edges, color, label):
  for node in nodes.keys():
#    label[nodes[node]["node"]]=node
    exp=nodes[node]["expression"]
    metrics["expression"][nodes[node]["node"]]=exp
    if exp=="intergenic":
      color[nodes[node]["node"]]=tlp.Color(255,255,255)
    elif exp=="NA":
      color[nodes[node]["node"]]=tlp.Color(0,0,0)
    elif exp=="up":
      color[nodes[node]["node"]]=tlp.Color(0,255,0)
    elif exp=="down":
      color[nodes[node]["node"]]=tlp.Color(255,0,0)
    elif exp=="stable":
      color[nodes[node]["node"]]=tlp.Color(125,125,125)
    if len(nodes[node]["reactome"])!=0:
      metrics["reactome"][nodes[node]["node"]]=nodes[node]["reactome"]
    if len(nodes[node]["metabo"])!=0:
      metrics["metabo"][nodes[node]["node"]]=nodes[node]["metabo"]
  for edge in edges.keys():
    if edges[edge]["interraction"]=="gain":
      color[edges[edge]["edge"]]=tlp.Color(0,255,0)
    if edges[edge]["interraction"]=="loss":
      color[edges[edge]["edge"]]=tlp.Color(255,0,0)
    if edges[edge]["interraction"]=="stable":
      color[edges[edge]["edge"]]=tlp.Color(200,200,200)
    metrics["distance"][edges[edge]["edge"]]=int(edges[edge]["distance"])
    metrics["interraction"][edges[edge]["edge"]]=edges[edge]["interraction"]
  return 0

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
  '''
  viewBorderColor = graph.getColorProperty("viewBorderColor")
  viewBorderWidth = graph.getDoubleProperty("viewBorderWidth")
  '''
  viewColor = graph.getColorProperty("viewColor")
  '''
  viewFont = graph.getStringProperty("viewFont")
  viewFontSize = graph.getIntegerProperty("viewFontSize")
  viewIcon = graph.getStringProperty("viewIcon")
  '''
  viewLabel = graph.getStringProperty("viewLabel")
  '''
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
  '''
  
  Metrics={}
  Metrics["expression"]=graph.getStringProperty("expression")
  Metrics["distance"]=graph.getDoubleProperty("distance")
  Metrics["interraction"]=graph.getStringProperty("interaction")
  Metrics["reactome"]=graph.getStringVectorProperty("reactome")
  Metrics["metabo"]=graph.getStringVectorProperty("metabo")

  locus1 = []
  locus2 = []
  interaction = []
  distances = []
  nodes = {}
  edges = {}
  interract_dict = {}

  interraction_dict=importData()
  create_nodes_edges(graph, interraction_dict, nodes)
  import_chromosome6_fragments_expressions(nodes, Metrics["expression"])
  import_metabos(nodes, interraction_dict)
  ajouter_metrics(nodes, Metrics, interraction_dict, viewColor, viewLabel)  
  
