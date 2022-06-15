#Francesco Gualdi

""" This script is for creating GUILD inputs from DisGeNet and HIPPIE by keeping
    only the largest connected component it will create in the folder input a nodes.sif 
    file and a interactome.sif that are the inputs of guild """

import pandas as pd
import argparse 
import numpy as np
import igraph as ig
import os


parser=argparse.ArgumentParser(description='parse the HIPPIE PPI dataset and DisGeNet database')
parser.add_argument('-hs','--hippie_score',type=float,help='Score to Filter PPI reliability')
parser.add_argument('-ds','--dis_score',type=float,help='Score threshold for Disgenet')
parser.add_argument('-c','--curated',type=str,help='if True uses the curated version of DisGeNet')
parser.add_argument('-m','--scoring_mode',type=str,default='dis',help='scoring mode for the nodes, available modes: dis, norm_dis, binary')
parser.add_argument('-d','--disease',help='CUI of the Disease to create the nodes for')
args=parser.parse_args()



def ParseHippie():
    
    hippie=pd.read_csv("../datasets/hippie_current.txt",sep='\t',low_memory=False)
    #Parse the columns headers
    hippie.columns=[x.replace(" ", "_") for x in hippie.columns]
    edges=hippie[hippie.Confidence_Value >= args.hippie_score ] 
    edges['ID_Interactor_A']=edges['ID_Interactor_A'].apply(lambda x: x.replace('entrez gene:',''))
    edges['ID_Interactor_B']=edges['ID_Interactor_B'].apply(lambda x: x.replace('entrez gene:',''))

    genes=pd.read_csv("../datasets/hs_genes.txt",sep='\t')

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
    ed.to_csv('../inputs/hippie_interactome.sif',index=False,header=False, sep=' ')#save the file 


def CreateNodes(Dis=args.disease,curated=args.curated,GdaScore=args.dis_score,mode=args.scoring_mode):
    
    
    ppi_network_edges=pd.read_csv('../inputs/hippie_interactome.sif',sep=' ')
    ppi_network_edges.columns=['source','score','target']
    
    if curated=='True':
        GdaDataset=pd.read_csv('../datasets/curated_gda.tsv',sep='\t')
        GeneColumn='Gene_id'
        ScoreColumn='Score_gda'
        DisColumn='Disease_id'
    else:
        GdaDataset=pd.read_csv('../datasets/all_gda.tsv',sep='\t')
        GeneColumn='geneId'
        ScoreColumn='score'
        DisColumn='diseaseId'
    
    GdaDataset=GdaDataset[(GdaDataset[ScoreColumn]>=GdaScore)&(GdaDataset[DisColumn]==Dis)]
    print(GdaDataset.iloc[0:11])
    seeds=GdaDataset[[GeneColumn, ScoreColumn]]
    print(seeds.iloc[0:11])
    seeds=seeds[(seeds[GeneColumn].isin(ppi_network_edges.source) | (seeds[GeneColumn].isin(ppi_network_edges.target)))]
    nodes=pd.DataFrame(set(ppi_network_edges.source.tolist() + ppi_network_edges.target.tolist()),columns=[GeneColumn])
    nodes_wo_seeds=pd.DataFrame([node for node in nodes[GeneColumn].tolist() if node not in seeds[GeneColumn].tolist()],columns=[GeneColumn])
    nodes_wo_seeds[ScoreColumn]=0
    if mode=='dis':
        nodes=pd.concat([seeds,nodes_wo_seeds])
        nodes.to_csv('../inputs/nodes.sif', sep=' ', index=False, header=False)
    elif mode=='norm_dis':
        nodes=pd.concat([seeds,nodes_wo_seeds])
        ScoreMax=nodes[ScoreColumn].max()
        ScoreMin=nodes[ScoreColumn].min()
        nodes[ScoreColumn]=nodes[ScoreColumn].apply(lambda x: (x-ScoreMin)/(ScoreMax-ScoreMin))
        nodes.to_csv('../inputs/nodes.sif', sep=' ', index=False, header=False)
    elif mode == 'binary':
        nodes=pd.concat([seeds,nodes_wo_seeds])
        nodes[ScoreColumn]=[1 if score > 0 else 0 for score in nodes[ScoreColumn].tolist()]
        nodes.to_csv('../inputs/nodes.sif', sep=' ', index=False, header=False)







if __name__=='__main__':
    if not os.path.isdir('../inputs'):
        os.mkdir('../inputs')
    ParseHippie()
    CreateNodes()





