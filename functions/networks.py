#!/usr/bin/env python
# coding: utf-8

# # PARSE HIPPIE INTERACTOME

# 1. Parse the database
# 2. checks if there's some genes that are deprecated
# 3. build a network
# 4. simplifies the network by removing the duplicated elements and edges
# 5. takes the biggest cluster in the network and get rid of all non connected element in the graph

import pandas as pd
import igraph as ig
import wrappers as wr
import os

def parse_hippie(save_out=False,path=os.getcwd()):
    
    hippie=pd.read_csv("http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/HIPPIE-current.mitab.txt",sep='\t')
    #Parse the columns headers
    hippie.columns=[x.replace(" ", "_") for x in hippie.columns]
    edges=hippie[hippie.Confidence_Value>0.65]
    edges['ID_Interactor_A']=edges['ID_Interactor_A'].apply(lambda x: x.replace('entrez gene:',''))
    edges['ID_Interactor_B']=edges['ID_Interactor_B'].apply(lambda x: x.replace('entrez gene:',''))

    genes=pd.read_csv("https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",sep='\t')

    def convert_to_str(integer):
        return str(integer)
    genes.GeneID=genes.GeneID.apply(lambda x: convert_to_str(x))

    #keep only the genes present in the human genes database
    edges = edges[(edges.ID_Interactor_A.isin(genes.GeneID.tolist()))&(edges.ID_Interactor_B.isin(genes.GeneID.tolist()))]

    #build a list of single genes to be used as nodes
    nodes=pd.DataFrame(set(edges.ID_Interactor_A.tolist()).union(edges.ID_Interactor_B.tolist()),columns=['Nodes'])

    graph_full=ig.Graph.DataFrame(edges, directed=False, vertices=nodes)
    graph_full.ecount(),graph_full.vcount()#(104249, 13631)

    simplified_graph=ig.Graph.simplify(graph_full)
    simplified_graph.ecount(),graph_full.vcount()#(101820, 13631)

    aux=ig.Graph.clusters(simplified_graph)#list the clusters of the graph 
    largest_cluster_index=list(aux.sizes()).index(max(list(aux.sizes()))) # gives the cluster number with the largest connected component
    nodes_in_largest_cluster=aux[largest_cluster_index]# returns the genes in the LCC

    gc1=ig.Graph.induced_subgraph(simplified_graph,nodes_in_largest_cluster) # creates a subgraph with the genes in the LCC
    gc1.vcount(), gc1.ecount() # 13443  101755

    #store the edges in the dataframe and replace the node id (number of the node from (0-N) with edge name (entrez of the given gene
    ve=gc1.get_vertex_dataframe()
    ed=gc1.get_edge_dataframe()
    ed['source'].replace(ve['name'], inplace=True)
    ed['target'].replace(ve['name'], inplace=True)

    ed.insert(1,'scores',1) #insert the score column between edges for guild
    if save_out:
        ed.to_csv(path+'/hippie_interactome.sif',index=False,header=False, sep=' ')#save the file 
    return ed, gc1

#build a list of single genes to be used as nodes

def parse_huri(save_out=False,path=os.getcwd()):
    huri=pd.read_csv('http://www.interactome-atlas.org/data/HuRI.tsv',sep='\t',names=['intA','intB'])

    nodes_gfull=pd.DataFrame(set(huri.intA.tolist()).union(huri.intB.tolist()),columns=['Nodes'])

    graph_full=ig.Graph.DataFrame(huri, directed=False, vertices=nodes_gfull)
    graph_full.ecount(),graph_full.vcount()#(52167, 8231)

    simplified_graph=ig.Graph.simplify(graph_full)
    simplified_graph.ecount(),graph_full.vcount()#(51671, 8231)

    aux=ig.Graph.clusters(simplified_graph)#list the clusters of the graph 

    largest_cluster_index=list(aux.sizes()).index(max(list(aux.sizes()))) # gives the cluster number with the largest connected component

    nodes_in_largest_cluster=aux[largest_cluster_index]# returns the ids of the nodes represented in the largest cluster

    gc1=ig.Graph.induced_subgraph(simplified_graph,nodes_in_largest_cluster) # creates a subgraph with the genes in the LCC
    gc1.vcount(), gc1.ecount() #(8108, 51619)

    #build the interaction dataframe
    ve=gc1.get_vertex_dataframe()
    ed=gc1.get_edge_dataframe()
    #replace source and target that initially are represented with nodes id with nodes names (genes_id)
    ed['source'].replace(ve['name'], inplace=True)
    ed['target'].replace(ve['name'], inplace=True)
    #insert the score column between edges for guild
    ed.insert(1,'scores',1)
    if save_out:
        ed.to_csv(path+'/huri_interactome.sif',index=False,header=False, sep=' ')#save the file 
    return ed,gc1





# # GET n MOST CONNECTED GENES IN THE NETWORK



def get_most_connected_nodes(network,n,entrez_to_names=False):
    name_degrees=[]
    for vertex in network.vs():
        name_degrees.append((vertex['name'],vertex.degree()))
    sorted_degrees=sorted(name_degrees, key=lambda tup: tup[1],reverse=True)
    lista_genes=[]
    lista_degrees=[]
    for element in sorted_degrees[0:n]:
        lista_genes.append(element[0])
        lista_degrees.append(element[1])
    if entrez_to_names:
        return pd.DataFrame({'id':lista_genes,
                             'degree':lista_degrees,
                            'names':[wr.gene_mapping(int(x),'entrez','symbol') for x in lista_genes]})
        
    return pd.DataFrame({'id':lista_genes,
                  'degree':lista_degrees})




