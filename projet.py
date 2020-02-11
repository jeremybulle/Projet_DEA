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

def ouvrir_fichier(fichier, enTete=True):
  #Ouvre un fichier .tsv et le renvoie sous forme d'une liste
  #fichier = nom du fichier à ouvrir
  #enTete = True si le fichier contient un en-tête
  #enTete = False si le fichier ne contient pas de en-tête
  #return données sous forme de liste.
  wd = '/net/stockage/TulipDEA/'
  with open (wd+fichier, newline ='') as csvfile:
    donnees = list(csv.reader(csvfile, delimiter ='\t'))
  if (enTete):
    donnees.pop(0)
  return donnees

def importDonnees():
  #importe et structure les données du fichier interactions_chromosome6.csv dans un dictionnaire d'arêtes
  #chaque clé correspond à une arête avec comme valeurs les noms des deux gènes en relation, le type et la distance d'interraction qui lie ces deux gènes
  #renvoit le dictionnaire des arêtes annotées
  interraction_dict ={}
  donnees = ouvrir_fichier('interactions_chromosome6.csv')
  longueur_donnees = len(donnees)
  for i in range(longueur_donnees):
    nom_interraction = donnees[i][1] + "_" + donnees[i][2]
    interraction_dict[nom_interraction] = {"locus1" : donnees[i][1], "locus2" : donnees[i][2], "interraction" : donnees[i][3], "distance": donnees[i][4]}
    i +=1
  return interraction_dict

def import_chromosome6_fragments_expressions(noeuds):
  #ajoute les informations des gènes préimplémentés contenues dans le fichier chromosome6_fragments_expressions.csv
  #noeuds = dictionnaire d'accès aux noeuds du graphe et à leurs annotations
  lecture=ouvrir_fichier('chromosome6_fragments_expressions.csv')
  longueur_lecture=len(lecture)
  for ligne in range(longueur_lecture):
    if lecture[ligne][1] in noeuds.keys():
      noeuds[lecture[ligne][1]]["expression"]=lecture[ligne][2]

def import_metabos(noeuds):
  #ajoute les informations des gènes préimplémentés contenues dans les fichiers KEGG.symbols.csv et REACTOME.symbols.csv
  #noeuds = dictionnaire d'accès aux noeuds du graphe et à leurs annotations
  Kegg={}
  lecture=ouvrir_fichier('KEGG.symbols.csv', False)
  for ligne in range(len(lecture)):
    Kegg[lecture[ligne][0]]={}
    for index_gene in range(2,len(lecture[ligne])):
      if lecture[ligne][index_gene] in noeuds.keys():
        Kegg[lecture[ligne][0]][lecture[ligne][index_gene]]=noeuds[lecture[ligne][index_gene]]
        noeuds[lecture[ligne][index_gene]]["metabo"][lecture[ligne][0]]=Kegg[lecture[ligne][0]]
  Reactome={}
  lecture=ouvrir_fichier('REACTOME.symbols.csv', False)
  for ligne in range(len(lecture)):
    Reactome[lecture[ligne][0]]={}
    for index_gene in range(2,len(lecture[ligne])):
      if lecture[ligne][index_gene] in noeuds.keys():
        Reactome[lecture[ligne][0]][lecture[ligne][index_gene]]=noeuds[lecture[ligne][index_gene]]
        noeuds[lecture[ligne][index_gene]]["reactome"][lecture[ligne][0]]=Reactome[lecture[ligne][0]]
  return Kegg, Reactome


def creer_noeuds_aretes(gr,aretes,noeuds):
  #Crée les noeuds et arêtes d'après le dictionnaire des arêtes annotées
  #et les rend accessible respectivement par le biais des disctionnaires des noeuds et des arêtes
  #gr = graphe sur lequel appliquer les opérations
  #aretes = dictionnaire des arêtes
  #noeuds = dictionnaire des noeuds
  for cle in aretes.keys():
    locus=aretes[cle]["locus1"]
    if locus not in noeuds.keys():
      noeuds[locus]={}
      noeuds[locus]={"noeud":gr.addNode(), "metabo":{}, "reactome":{}}
    locus2=aretes[cle]["locus2"]
    if locus2 not in noeuds.keys():
      noeuds[locus2]={}
      noeuds[locus2]={"noeud":gr.addNode(), "metabo":{}, "reactome":{}}
    aretes[cle]["arete"]=gr.addEdge(noeuds[locus]["noeud"],noeuds[locus2]["noeud"])

def ajouter_metriques(noeuds, metriques, aretes):
  #ajoute toutes les annotations des différents fichiers sur les noeuds/arêtes
  #noeuds = dictionnaire des noeuds
  #aretes = dictionnaire des arêtes
  #metriques = liste des différentes metriques
  for noeud in noeuds.keys():
    list_reac=[]
    list_meta=[]
    metriques["ID"][noeuds[noeud]["noeud"]]=noeud
    exp=noeuds[noeud]["expression"]
    metriques["expression"][noeuds[noeud]["noeud"]]=exp
    if exp=="intergenic":
      metriques["color"][noeuds[noeud]["noeud"]]=tlp.Color(150,150,150)
    elif exp=="NA":
      metriques["color"][noeuds[noeud]["noeud"]]=tlp.Color(0,0,0)
    elif exp=="up":
      metriques["color"][noeuds[noeud]["noeud"]]=tlp.Color(0,255,0)
    elif exp=="down":
      metriques["color"][noeuds[noeud]["noeud"]]=tlp.Color(255,0,0)
    elif exp=="stable":
      metriques["color"][noeuds[noeud]["noeud"]]=tlp.Color(0,100,255)
    if len(noeuds[noeud]["reactome"].keys())!=0:
      for reac in noeuds[noeud]["reactome"].keys():
        list_reac.append(reac)
      metriques["reactome"][noeuds[noeud]["noeud"]]=list_reac
    if len(noeuds[noeud]["metabo"].keys())!=0:
      for meta in noeuds[noeud]["metabo"].keys():
        list_meta.append(meta)
      metriques["metabo"][noeuds[noeud]["noeud"]]=list_meta
  for arete in aretes.keys():
    metriques["ID"][aretes[arete]["arete"]]=arete
    if aretes[arete]["interraction"]=="gain":
      metriques["color"][aretes[arete]["arete"]]=tlp.Color(0,255,0)
    if aretes[arete]["interraction"]=="loss":
      metriques["color"][aretes[arete]["arete"]]=tlp.Color(255,0,0)
    if aretes[arete]["interraction"]=="stable":
      metriques["color"][aretes[arete]["arete"]]=tlp.Color(200,200,200)
    metriques["distance"][aretes[arete]["arete"]]=int(aretes[arete]["distance"])
    metriques["distanceGraph"][aretes[arete]["arete"]]=int(aretes[arete]["distance"])/100000
    metriques["interraction"][aretes[arete]["arete"]]=aretes[arete]["interraction"]
  return 0

def creer_sousGraph(gr, metriques, noeuds):
  #Creer un sous graphe de gr contenant uniquement les gènes dont l'expression est impactée par l'activation d'un lymphocyte T.
  #gr = graphe sur lequel appliquer les opérations
  #metriques = liste des différentes metriques
  #noeuds = dictionnaire des noeuds
  list_noeuds=[]
  for noeud in noeuds.keys():
    if (len(noeuds[noeud]["metabo"].keys())!=0 or len(noeuds[noeud]["reactome"].keys())!=0) and (noeuds[noeud]["expression"]=="up" or noeuds[noeud]["expression"]=="down"):
      list_noeuds.append(noeuds[noeud]["noeud"])
  gr.inducedSubGraph(list_noeuds, gr, "upDown")

def creer_subReac(gr, noeuds, Keggs, Metabos, layout):
  #Creer des sous graphes pour chaque réaction dont au moins un gène qui le compose est impactée par l'activation d'un lymphocyte T
  #gr = graphe sur lequel appliquer les opérations
  #noeuds = dictionnaire des noeuds
  #Keggs = dictionnaire structurant les données du fichier Kegg.
  #Metabos = dictionnaire structurant les données du fichier Reactome
  #layout = Propriété modifiée pour la visualisation
  for reaction in Keggs.keys():
    listNoeuds=[]
    if len(Keggs[reaction].keys())!=0:
      for gene in Keggs[reaction].keys():
        if noeuds[gene]["expression"] in ["up","down"]:
          listNoeuds.append(noeuds[gene]["noeud"])
          for noeud in gr.getInOutNodes(noeuds[gene]["noeud"]):
            listNoeuds.append(noeud)
        else :
          continue
      if len(listNoeuds)!=0:
        sousGraph=gr.inducedSubGraph(listNoeuds, gr, reaction)
        sousGraph.applyLayoutAlgorithm("FM^3 (OGDF)", layout)
  for reaction in Metabos.keys():
    listNoeuds=[]
    if len(Metabos[reaction].keys())!=0:
      for gene in Metabos[reaction].keys():
        if noeuds[gene]["expression"] in ["up","down"]:
          listNoeuds.append(noeuds[gene]["noeud"])
          for noeud in gr.getInOutNodes(noeuds[gene]["noeud"]):
            listNoeuds.append(noeud)
        else :
          continue
      if len(listNoeuds)!=0:
        sousGraph=gr.inducedSubGraph(listNoeuds, gr, reaction)
        sousGraph.applyLayoutAlgorithm("FM^3 (OGDF)", layout)
  return True

def creer_express_graph(gr, noeuds, Metriques):
  #Crée un graph pour chaque gène dont l'expression est impactée par l'activation d'un lymphocyte T et tous les gènes en intéractions entre eux qui sont eux aussi impactées par l'activation
  #gr = graphe sur lequel appliquer les opérations
  #noeuds = dictionnaire des noeuds
  #Metriques = liste des différentes metriques
  #return clone = renvoie un sous graphe pour ne pas modifier le graphe d'origine (ne fonctionne pas)
  liste_noeud=[]
  for noeud in gr.getNodes():
    liste_noeud.append(noeud)
  clone=gr.addCloneSubGraph("clone", False, True)
  noeuds_vu=[]
  for noeud in gr.getNodes():
    List_noeud=[]
    if Metriques["expression"][noeud] in ["up","down"] and (len(Metriques["reactome"][noeud])!=0 or len(Metriques["metabo"][noeud])!=0) and (noeud not in noeuds_vu):
      noeuds_vu.append(noeud)
      explorerUpDown(gr, noeuds, noeud, List_noeud, Metriques, noeuds_vu)
    if len(List_noeud)>1:
      clone.inducedSubGraph(List_noeud, clone, Metriques["ID"][noeud])
  return clone

def explorerUpDown(gr, noeuds, noeud, List_noeud, Metriques, noeuds_vu):
  #Parcourt tous les voisins d'un noeud si celui-ci est impacté par l'activation d'un lymphocyte T et si il est impliqué dans une voie métabolique
  #gr = graphe sur lequel appliquer les opérations
  #noeuds = dictionnaire des noeuds
  #List_noeud = accumulateur de l'exploration
  #Metriques = liste des différentes metriques
  #noeuds_vu = liste de noeuds déjà parcouru par l'algorythme afin d'éviter de dupliquer les graphes  
  List_noeud.append(noeud)
  for voisin in gr.getInOutNodes(noeud):
    noeuds_vu.append(voisin)
    if Metriques["expression"][voisin] in ["up","down"] and voisin not in List_noeud and (len(Metriques["reactome"][voisin])!=0 or len(Metriques["metabo"][voisin])!=0):
      explorerUpDown(gr, noeuds, voisin, List_noeud, Metriques, noeuds_vu)
  
def etiquettes(gr, noeuds, aretes, Metriques):
  #Ajoute des étiquettes aux gènes correspondant aux voies métaboliques dans lesquels ils sont impliqués
  #gr = graphe sur lequel appliquer les opérations
  #noeuds = dictionnaire des noeuds
  #aretes = dictionnaire des arêtes
  #Metriques = liste des différentes metriques
  for sousGraph in gr.subGraphs():
    for noeud in sousGraph.getNodes():
      if len(Metriques["reactome"][noeud])!=0:
        for reac in Metriques["reactome"][noeud] :
          noeuds[reac]={"noeud":sousGraph.addNode()}
          Metriques["ID"][noeuds[reac]["noeud"]]=reac
          Metriques["label"][noeuds[reac]["noeud"]]=reac
          aretes[reac+"_"+Metriques["ID"][noeud]]={}
          aretes[reac+"_"+Metriques["ID"][noeud]]["arete"]=sousGraph.addEdge(noeuds[reac]["noeud"],noeud)
      if len(Metriques["metabo"][noeud])!=0:
        for reac in Metriques["metabo"][noeud] :
          noeuds[reac]={"noeud":sousGraph.addNode()}
          Metriques["ID"][noeuds[reac]["noeud"]]=reac
          Metriques["label"][noeuds[reac]["noeud"]]=reac
          aretes[reac+"_"+Metriques["ID"][noeud]]={}
          aretes[reac+"_"+Metriques["ID"][noeud]]["arete"]=sousGraph.addEdge(noeuds[reac]["noeud"],noeud)
    viewLayout=sousGraph.getLayoutProperty("viewLayout")
    sousGraph.applyLayoutAlgorithm("FM^3 (OGDF)", viewLayout)
  
def afficher_reseau(gr, Metriques, noeuds, Keggs, Metabos, viewLayout):
  #Afin de gérer les différentes modifications des graphes, afficher_reseau et afficher_interactions_reseaux ne peuvent être lancées ensembles et doivent être commenté afin de visualiser différents aspects du graph
  #Afficher réseau permet d'afficher le graphe des gènes d'intérêt
  #gr = graphe sur lequel appliquer les opérations
  #Metriques = liste des différentes metriques
  #noeuds = dictionnaire des noeuds
  #Keggs = dictionnaire structurant les données du fichier Kegg.
  #Metabos = dictionnaire structurant les données du fichier Reactome
  #viewLayout = Propriété modifiée pour la visualisation
  creer_sousGraph(gr, Metriques, noeuds)
  creer_subReac(gr, noeuds, Keggs, Metabos, viewLayout)

def afficher_interactions_reseaux(gr, noeuds, Metriques, aretes):
  #Afin de gérer les différentes modifications des graphes, afficher_reseau et afficher_interactions_reseaux ne peuvent être lancées ensembles et doivent être commenté afin de visualiser différents aspects du graph
  #Afficher réseau permet d'afficher le graphe des gènes d'intérêt
  #gr = graphe sur lequel appliquer les opérations
  #Metriques = liste des différentes metriques
  #noeuds = dictionnaire des noeuds  
  #aretes = dictionnaire des arêtes
  clone=creer_express_graph(gr, noeuds, Metriques)
  etiquettes(clone, noeuds, aretes, Metriques)


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
  viewLayout = graph.getLayoutProperty("viewLayout")
  Metriques={}
  Metriques["color"] = graph.getColorProperty("viewColor")
  Metriques["label"] = graph.getStringProperty("viewLabel")
  Metriques["ID"]=graph.getStringProperty("ID")
  Metriques["expression"]=graph.getStringProperty("expression")
  Metriques["distance"]=graph.getDoubleProperty("distance")
  Metriques["distanceGraph"]=graph.getDoubleProperty("distanceGraph")
  Metriques["interraction"]=graph.getStringProperty("interaction")
  Metriques["reactome"]=graph.getStringVectorProperty("reactome")
  Metriques["metabo"]=graph.getStringVectorProperty("metabo")

  noeuds = {}
  interract_dict = {}
  
  interraction_dict=importDonnees()
  creer_noeuds_aretes(graph, interraction_dict, noeuds)
  import_chromosome6_fragments_expressions(noeuds)
  Keggs, Metabos=import_metabos(noeuds)
  ajouter_metriques(noeuds, Metriques, interraction_dict)  
 
  params = tlp.getDefaultPluginParameters("FM^3 (OGDF)", graph)
  params["Edge Length Property"] = Metriques["distanceGraph"]
  success = graph.applyLayoutAlgorithm("FM^3 (OGDF)", viewLayout, params)
  
  viewLayout.perfectAspectRatio()
  updateVisualization(True)
  original=graph.addCloneSubGraph("original", False, False)
  sub = graph.addCloneSubGraph("sub", False, False)
  
  #afficher_reseau(sub, Metriques, noeuds, Keggs, Metabos, viewLayout)

  afficher_interactions_reseaux(sub, noeuds, Metriques, interraction_dict)
