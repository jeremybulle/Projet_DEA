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

def ouvrir_fichier(fichier, header=True):
  #Ouvre un fichier .tsv et le renvoie sous forme d'une liste
  #fichier = nom du fichier à ouvrir
  #header = True si le fichier contient un header
  #header = False si le fichier ne contient pas de header
  #return données sous forme de liste.
  wd = '/net/stockage/TulipDEA/'
  with open (wd+fichier, newline ='') as csvfile:
    data = list(csv.reader(csvfile, delimiter ='\t'))
  if (header):
    data.pop(0)
  return data

def importData():
  #importe et structure les données du fichier interactions_chromosome6.csv dans un dictionnaire d'arêtes
  #chaque clé correspond à une arête avec comme valeurs les noms des deux gènes en relation, le type et la distance d'interraction qui lie ces deux gènes
  #renvoit le dictionnaire des arêtes annotées
  interraction_dict ={}
  data = ouvrir_fichier('interactions_chromosome6.csv')
  lendata = len(data)
  for i in range(lendata):
    nom_interraction = data[i][1] + "_" + data[i][2]
    interraction_dict[nom_interraction] = {"locus1" : data[i][1], "locus2" : data[i][2], "interraction" : data[i][3], "distance": data[i][4]}
    i +=1
  return interraction_dict

def import_chromosome6_fragments_expressions(nodes):
  #ajoute les informations des gènes préimplémentés contenues dans le fichier chromosome6_fragments_expressions.csv
  #nodes = dictionnaire d'accès aux noeuds du graphe et à leurs annotations
  read=ouvrir_fichier('chromosome6_fragments_expressions.csv')
  lenread=len(read)
  for ligne in range(lenread):
    if read[ligne][1] in nodes.keys():
      nodes[read[ligne][1]]["expression"]=read[ligne][2]

def import_metabos(nodes):
  #ajoute les informations des gènes préimplémentés contenues dans les fichiers KEGG.symbols.csv et REACTOME.symbols.csv
  #nodes = dictionnaire d'accès aux noeuds du graphe et à leurs annotations
  Kegg={}
  read=ouvrir_fichier('KEGG.symbols.csv', False)
  for ligne in range(len(read)):
    Kegg[read[ligne][0]]={}
    for index_gene in range(2,len(read[ligne])):
      if read[ligne][index_gene] in nodes.keys():
        Kegg[read[ligne][0]][read[ligne][index_gene]]=nodes[read[ligne][index_gene]]
        nodes[read[ligne][index_gene]]["metabo"][read[ligne][0]]=Kegg[read[ligne][0]]
  Reactome={}
  read=ouvrir_fichier('REACTOME.symbols.csv', False)
  for ligne in range(len(read)):
    Reactome[read[ligne][0]]={}
    for index_gene in range(2,len(read[ligne])):
      if read[ligne][index_gene] in nodes.keys():
        Reactome[read[ligne][0]][read[ligne][index_gene]]=nodes[read[ligne][index_gene]]
        nodes[read[ligne][index_gene]]["reactome"][read[ligne][0]]=Reactome[read[ligne][0]]
  return Kegg, Reactome


def create_nodes_edges(gr,edges,nodes):
  #Crée les noeuds et arêtes d'après le dictionnaire des arêtes annotées
  #et les rend accessible respectivement par le biais des disctionnaires des noeuds et des arêtes
  #gr = graphe sur lequel appliquer les opérations
  #edges = dictionnaire des arêtes
  #nodes = dictionnaire des noeuds
  for key in edges.keys():
    locus=edges[key]["locus1"]
    if locus not in nodes.keys():
      nodes[locus]={}
      nodes[locus]={"node":gr.addNode(), "metabo":{}, "reactome":{}}
    locus2=edges[key]["locus2"]
    if locus2 not in nodes.keys():
      nodes[locus2]={}
      nodes[locus2]={"node":gr.addNode(), "metabo":{}, "reactome":{}}
    edges[key]["edge"]=gr.addEdge(nodes[locus]["node"],nodes[locus2]["node"])

def ajouter_metrics(nodes, metrics, edges):
  #ajoute toutes les annotations des différents fichiers sur les noeuds/arêtes
  #nodes = dictionnaire des noeuds
  #edges = dictionnaire des arêtes
  #metrics = liste des différentes metrics
  for node in nodes.keys():
    list_reac=[]
    list_meta=[]
    metrics["ID"][nodes[node]["node"]]=node
    exp=nodes[node]["expression"]
    metrics["expression"][nodes[node]["node"]]=exp
    if exp=="intergenic":
      metrics["color"][nodes[node]["node"]]=tlp.Color(150,150,150)
    elif exp=="NA":
      metrics["color"][nodes[node]["node"]]=tlp.Color(0,0,0)
    elif exp=="up":
      metrics["color"][nodes[node]["node"]]=tlp.Color(0,255,0)
    elif exp=="down":
      metrics["color"][nodes[node]["node"]]=tlp.Color(255,0,0)
    elif exp=="stable":
      metrics["color"][nodes[node]["node"]]=tlp.Color(0,100,255)
    if len(nodes[node]["reactome"].keys())!=0:
      for reac in nodes[node]["reactome"].keys():
        list_reac.append(reac)
      metrics["reactome"][nodes[node]["node"]]=list_reac
    if len(nodes[node]["metabo"].keys())!=0:
      for meta in nodes[node]["metabo"].keys():
        list_meta.append(meta)
      metrics["metabo"][nodes[node]["node"]]=list_meta
  for edge in edges.keys():
    metrics["ID"][edges[edge]["edge"]]=edge
    if edges[edge]["interraction"]=="gain":
      metrics["color"][edges[edge]["edge"]]=tlp.Color(0,255,0)
    if edges[edge]["interraction"]=="loss":
      metrics["color"][edges[edge]["edge"]]=tlp.Color(255,0,0)
    if edges[edge]["interraction"]=="stable":
      metrics["color"][edges[edge]["edge"]]=tlp.Color(200,200,200)
    metrics["distance"][edges[edge]["edge"]]=int(edges[edge]["distance"])
    metrics["distanceGraph"][edges[edge]["edge"]]=int(edges[edge]["distance"])/100000
    metrics["interraction"][edges[edge]["edge"]]=edges[edge]["interraction"]
  return 0

def create_subgraph(gr, metrics, nodes):
  list_nodes=[]
  for node in nodes.keys():
    if (len(nodes[node]["metabo"].keys())!=0 or len(nodes[node]["reactome"].keys())!=0) and (nodes[node]["expression"]=="up" or nodes[node]["expression"]=="down"):
      list_nodes.append(nodes[node]["node"])
  gr.inducedSubGraph(list_nodes, gr, "upDown")

def create_subReac(gr, nodes, Keggs, Metabos, layout):
  for reaction in Keggs.keys():
    listNodes=[]
    if len(Keggs[reaction].keys())!=0:
      for gene in Keggs[reaction].keys():
        if nodes[gene]["expression"] in ["up","down"]:
          listNodes.append(nodes[gene]["node"])
          for node in gr.getInOutNodes(nodes[gene]["node"]):
            listNodes.append(node)
        else :
          continue
      if len(listNodes)!=0:
        subgraph=gr.inducedSubGraph(listNodes, gr, reaction)
        subgraph.applyLayoutAlgorithm("FM^3 (OGDF)", layout)
  for reaction in Metabos.keys():
    listNodes=[]
    if len(Metabos[reaction].keys())!=0:
      for gene in Metabos[reaction].keys():
        if nodes[gene]["expression"] in ["up","down"]:
          listNodes.append(nodes[gene]["node"])
          for node in gr.getInOutNodes(nodes[gene]["node"]):
            listNodes.append(node)
        else :
          continue
      if len(listNodes)!=0:
        subgraph=gr.inducedSubGraph(listNodes, gr, reaction)
        subgraph.applyLayoutAlgorithm("FM^3 (OGDF)", layout)

  return True

def creer_express_graph(gr, nodes, Metrics):
  liste_noeud=[]
  for noeud in gr.getNodes():
    liste_noeud.append(noeud)
  clone=gr.addCloneSubGraph("clone", False, True)
  noeuds_vu=[]
  for noeud in gr.getNodes():
    List_noeud=[]
    if Metrics["expression"][noeud] in ["up","down"] and (len(Metrics["reactome"][noeud])!=0 or len(Metrics["metabo"][noeud])!=0) and (noeud not in noeuds_vu):
      noeuds_vu.append(noeud)
      explorerUpDown(gr, nodes, noeud, List_noeud, Metrics, noeuds_vu)
    if len(List_noeud)>1:
      clone.inducedSubGraph(List_noeud, clone, Metrics["ID"][noeud])
  return clone

def explorerUpDown(gr, nodes, noeud, List_noeud, Metrics, noeuds_vu):
  List_noeud.append(noeud)
  for voisin in gr.getInOutNodes(noeud):
    noeuds_vu.append(voisin)
    if Metrics["expression"][voisin] in ["up","down"] and voisin not in List_noeud and (len(Metrics["reactome"][voisin])!=0 or len(Metrics["metabo"][voisin])!=0):
      explorerUpDown(gr, nodes, voisin, List_noeud, Metrics, noeuds_vu)
  
def etiquettes(gr, nodes, edges, Metrics):
  for subgraph in gr.subGraphs():
    for noeud in subgraph.getNodes():
      if len(Metrics["reactome"][noeud])!=0:
        for reac in Metrics["reactome"][noeud] :
          nodes[reac]={"node":subgraph.addNode()}
          Metrics["ID"][nodes[reac]["node"]]=reac
          Metrics["label"][nodes[reac]["node"]]=reac
          edges[reac+"_"+Metrics["ID"][noeud]]={}
          edges[reac+"_"+Metrics["ID"][noeud]]["edge"]=subgraph.addEdge(nodes[reac]["node"],noeud)
      if len(Metrics["metabo"][noeud])!=0:
        for reac in Metrics["metabo"][noeud] :
          nodes[reac]={"node":subgraph.addNode()}
          Metrics["ID"][nodes[reac]["node"]]=reac
          Metrics["label"][nodes[reac]["node"]]=reac
          edges[reac+"_"+Metrics["ID"][noeud]]={}
          edges[reac+"_"+Metrics["ID"][noeud]]["edge"]=subgraph.addEdge(nodes[reac]["node"],noeud)
    viewLayout=subgraph.getLayoutProperty("viewLayout")
    subgraph.applyLayoutAlgorithm("FM^3 (OGDF)", viewLayout)
  
def afficher_reseau(gr, Metrics, nodes, Keggs, Metabos, viewLayout):
  create_subgraph(gr, Metrics, nodes)
  create_subReac(gr, nodes, Keggs, Metabos, viewLayout)

def afficher_interactions_reseaux(gr, nodes, Metrics, edges):
  clone=creer_express_graph(gr, nodes, Metrics)
  etiquettes(clone, nodes, edges, Metrics)


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
  viewFont = graph.getStringProperty("viewFont")
  viewFontSize = graph.getIntegerProperty("viewFontSize")
  viewIcon = graph.getStringProperty("viewIcon")
  viewLabelBorderColor = graph.getColorProperty("viewLabelBorderColor")
  viewLabelBorderWidth = graph.getDoubleProperty("viewLabelBorderWidth")
  viewLabelColor = graph.getColorProperty("viewLabelColor")
  viewLabelPosition = graph.getIntegerProperty("viewLabelPosition")
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
  viewLayout = graph.getLayoutProperty("viewLayout")
  
  Metrics={}
  Metrics["color"] = graph.getColorProperty("viewColor")
  Metrics["label"] = graph.getStringProperty("viewLabel")
  Metrics["ID"]=graph.getStringProperty("ID")
  Metrics["expression"]=graph.getStringProperty("expression")
  Metrics["distance"]=graph.getDoubleProperty("distance")
  Metrics["distanceGraph"]=graph.getDoubleProperty("distanceGraph")
  Metrics["interraction"]=graph.getStringProperty("interaction")
  Metrics["reactome"]=graph.getStringVectorProperty("reactome")
  Metrics["metabo"]=graph.getStringVectorProperty("metabo")

  nodes = {}
  interract_dict = {}
  
  
  
  interraction_dict=importData()
  create_nodes_edges(graph, interraction_dict, nodes)
  import_chromosome6_fragments_expressions(nodes)
  Keggs, Metabos=import_metabos(nodes)
  ajouter_metrics(nodes, Metrics, interraction_dict)  
 
  params = tlp.getDefaultPluginParameters("FM^3 (OGDF)", graph)
  params["Edge Length Property"] = Metrics["distanceGraph"]
  success = graph.applyLayoutAlgorithm("FM^3 (OGDF)", viewLayout, params)
  
  viewLayout.perfectAspectRatio()
  updateVisualization(True)
  original=graph.addCloneSubGraph("original", False, False)
  sub = graph.addCloneSubGraph("sub", False, False)
  
  #afficher_reseau(sub, Metrics, nodes, Keggs, Metabos, viewLayout)

  #afficher_interactions_reseaux(sub, nodes, Metrics, interraction_dict)
