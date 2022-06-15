#!/usr/bin/env python
# coding: utf-8

# # Plot degree distribution
import matplotlib.pyplot as plt
import numpy as np
import re
import math
import matplotlib as mpl
import networkx as nx
import igraph as ig
import numpy as np
from collections import Counter

def plot_degree_distribution(network,save_fig=False):
    dd=list(network.degree_distribution().bins())
    degree=[]
    frequency=[]
    frequency0=[]
    degree0=[]
    for x in dd:
        if x[0]==0:
            degree0.append(int(x[0]))
            frequency0.append(int(x[2]))
        else:
            degree.append(int(x[0]))
            frequency.append(int(x[2]))
            
    fig, ax= plt.subplots(figsize=(20,10))
    datax=np.log10(degree)
    datay=np.log10(frequency)
    ax.scatter(datax,datay, alpha=0.5,s=300,edgecolor='b')
    
    ##CHOOSE THE DEGREES TICKS TO DISPLAY                    
    ax.set_yticks(np.log10(np.geomspace(1,max(frequency),20,dtype=int)))
    ax.set_yticklabels(np.geomspace(1,max(frequency),20,dtype=int),fontsize=12)

    
    ##CHOOSE THE DEGREES TICKS TO DISPLAY
    ax.set_xticks(np.log10(np.geomspace(1,max(degree),20,dtype=int)))
    ax.set_xticklabels(np.geomspace(1,max(degree),20,dtype=int),fontsize=12)
    
    #PLOT DEGREE == 0  NODES  IF THERE ARE ANY
    if len(degree0)>0:
        ax.scatter(-0.1,np.log10(frequency0),alpha=0.5,s=300,edgecolor='red')
    ax.grid()
    ax.set_xlabel('Degree',fontsize=20)
    ax.set_ylabel('Frequency',fontsize=20)
    ax.set_xlim(-0.112)
    if save_fig:
        plt.savefig('dd.jpeg',dpi=400,bbox_inches='tight')
    plt.show()
    

def plot_degree_distribution_nx(graph,save_fig=False):
    degree=[val for (node, val) in graph.degree()]
    freqdict=Counter(degree)
    frequency=[]
    frequency0=[]
    degree0=[]
    for x in degree:
        if x==0:
            degree0.append(0)
            frequency0.append(freqdict[x])
        else:
            frequency.append(freqdict[x])

    fig, ax= plt.subplots(figsize=(20,10))
    datax=np.log10([elem for elem in degree if elem!=0])
    datay=np.log10(frequency)

    #PLOT DEGREE == 0  NODES  IF THERE ARE ANY
    if len(frequency0)>0:
        ax.scatter(datax,datay, alpha=0.5,s=300,edgecolor='b')
        ax.scatter(-0.1,np.log10(freqdict[0]),alpha=0.5,s=300,edgecolor='red')
    else:
        ax.scatter(datax,datay, alpha=0.5,s=300,edgecolor='b')

    ##CHOOSE THE DEGREES TICKS TO DISPLAY                    
    ax.set_yticks(np.log10(np.geomspace(1,max(frequency),20,dtype=int)))
    ax.set_yticklabels(np.geomspace(1,max(frequency),20,dtype=int),fontsize=12)


    ##CHOOSE THE DEGREES TICKS TO DISPLAY
    ax.set_xticks(np.log10(np.geomspace(1,max(degree),20,dtype=int)))
    ax.set_xticklabels(np.geomspace(1,max(degree),20,dtype=int),fontsize=12)


    ax.grid()
    ax.set_xlabel('Degree',fontsize=20)
    ax.set_ylabel('Frequency',fontsize=20)
    ax.set_xlim(-0.150)
    plt.show()
    if save_fig:
        plt.savefig('dd',dpi=300)

# # Plot pathway enrichment gprofiler 

# ### Plot single graph with annotations

#takes in input data coordinates and labels it preventing the overlapping

def labelling_without_overlapping(x,y,list_of_annotations,ax,**kwargs):
    
    class Point:
        def __init__(self, x, y):
            self.x = x
            self.y = y
    
    
    def doOverlap(ret1,ret2):
        l1 = Point(ret1[0,0],ret1[1,1])
        r1 = Point(ret1[1,0],ret1[0,1])
        l2 = Point(ret2[0,0],ret2[1,1])
        r2 = Point(ret2[1,0],ret2[0,1])

        # If one rectangle is on left side of other
        if l1.x >= r2.x or l2.x >= r1.x:
            return False

        # If one rectangle is above other
        if(r1.y >= l2.y or r2.y >= l1.y):
            return False

        return True

    annotations_coord=[]
    for i, dot in enumerate(y):
        x_coords=x[i]
        y_coords=y[i]
        annotation=ax.annotate(str(list_of_annotations[i]),
                                xy=(x[i],y[i]),
                                 xytext=(x_coords,y_coords),
                                    **kwargs)

        ax.figure.canvas.draw()
        bbox=matplotlib.text.Text.get_window_extent(annotation)
        bbox_data = ax.transData.inverted().transform(bbox)
        factor=0.2*(bbox_data[0,0]-bbox_data[1,0])
        annotations_coord.append(bbox_data)
        ##BUILD THE SPIRAL##
        theta=np.radians(np.linspace(1,360*200,500))
        r=np.linspace(0,max(max(zip(x,y))),len(theta))
        x_2 = r*np.cos(theta)+x_coords#move the spiral onto the data point
        y_2 = r*np.sin(theta)+y_coords
        n=0
        keep_cycling=True
        while keep_cycling:
            keep_cycling=False
            #print('start checking box %s'% i)
            for ind, box in enumerate (annotations_coord[0:-1]):
                #print('checking %s and %s' % (i,ind))
                if doOverlap(box,bbox_data):
                    #print('%s and %s overlap' % (i,ind))
                    annotation.set_x(x_2[n])
                    annotation.set_y(y_2[n])
                    n+=1
                    ax.figure.canvas.draw()
                    bbox=matplotlib.text.Text.get_window_extent(annotation)
                    bbox_data = ax.transData.inverted().transform(bbox)
                    annotations_coord.pop()
                    annotations_coord.append(bbox_data)

                    #print('new coords (x=%i,y=%i)'%(x_coords,y_coords))
                    #print('new bbox data',bbox_data)
                    #print('annotation coordinates',box)
                    #print('restart iteration')
                    keep_cycling=True
                    break






def plot_enrichment_analisys(df,p_value):
    cmap = plt.get_cmap("tab10")
    colors=cmap(np.arange(len(df.source.value_counts())))
    ylim=max(-np.log10(df.p_value.tolist()))
    for i, (s,v) in enumerate(zip(df.source.value_counts().index,df.source.value_counts())):
        x=sorted(np.random.uniform(low=0.0, high=1, size=v))
        y=-np.log10(df[df.source==s].p_value.tolist())
        fig,ax=plt.subplots(figsize=(15,10))
        scat=ax.scatter(x,y,color=colors[i],alpha=0.4,
                    s=np.clip(df[df.source==s].term_size.tolist(),200,2000))
        ax.set_xlabel(s,fontsize=20)
        ax.set_ylabel('$-\log_{10}(p\ value)$',fontsize=20)
        ax.spines['bottom'].set_color(colors[i])
        ax.spines['bottom'].set_linewidth(4)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.yaxis.grid(color=colors[i], linestyle='-', linewidth=2,alpha=0.05)
        ax.set_xticks([])
        ax.tick_params(direction='out',which= 'major' ,length=8, width=2, left=False)
        ax.set_ylim(0,max(-np.log10(df.p_value.tolist())))
        
        
        
        ## ADD LABELS
        bbox_props = dict(boxstyle="round", fc="w", ec="k", lw=1)
        handles=[]
        lables=[]
        linea_y=np.linspace(ax.get_ylim()[1],0,
                            sum(map(lambda x: -np.log10(x)>=p_value,df[df.source==s].p_value.tolist())))
        h=0
        for ind, pv in enumerate(y):
            factory=np.random.uniform(-2,5,1)
            factorx=np.random.uniform(-0.1,0.1,1)
            if pv > p_value:
                handles.append(str(h))
                lables.append(df[df.source==s].name.tolist()[ind])
                an=ax.annotate(str(df[df.source==s].name.tolist()[ind]),
                                 xy=(x[ind],y[ind]),
                                 xytext=(ax.get_xlim()[1]+0.2, linea_y[ind]),
                                 #arrowprops=dict(arrowstyle="-",shrinkB=10),
                               bbox=bbox_props,
                                    fontsize=12)
                ax.annotate('',
                                 xy=(x[ind],y[ind]),
                                 xytext=(ax.get_xlim()[1]+0.2, linea_y[ind]),
                                 arrowprops=dict(arrowstyle="-",shrinkB=10,linewidth=0.2),
                               #bbox=bbox_props,
                                    fontsize=12)
                
                h+=1
        
        
        plt.show()
        
    
        
def plot_enrichment_analisys_network(df,pvalue,colormap='cividis',edgecolor='red',mkcolor='grey',mkfsize=15000,layout='spring',
                                     mklinewidths=4,alpha=1,figsize=(40,20),savefig=False,factor=1,k=10,
                                     cbarfontsize=10,labelling=True,**kwargs):
    maxpv=max([-np.log10(p) for p in df.p_value.tolist()])
    for i, (s,v) in enumerate(zip(df.source.value_counts().index,df.source.value_counts())):
        data=df[(df.source==s)&(-np.log10(df.p_value)>pvalue)].reset_index()
        if data.shape[0]==0:
            continue
        else:
            
            nxen=nx.Graph()
            #add nodes
            for i,r in data.iterrows():
                 nxen.add_node(r['name'],size=r['intersection_size'],pvalue=-np.log10(r['p_value']))

            #add edges
            for i,r in data.iterrows():
                for index,row in data.iloc[i+1:].reset_index().iterrows():
                    if len(set(r['intersections']).intersection(set(row['intersections'])))>0:
                        nxen.add_edge(r['name'],
                                    row['name'], 
                                    weight= len(set(r['intersections']).intersection(set(row['intersections']))))
            # Get positions for the nodes in G
            if layout=='spring':
                pos_ = nx.spring_layout(nxen,k)
            
            elif layout=='auto':
                ig_subgraph=ig.Graph.from_networkx(nxen)
                pos_= dict(zip([v['_nx_name'] for v in ig_subgraph.vs],[coord for coord in ig_subgraph.layout_auto()]))
                    
                    


            

            #Normalize connections
            connections=[]
            for edge in nxen.edges(data=True):
                connections.append(edge[2]['weight'])
            if len(connections)!=0:
                
                if ((max(connections)-min(connections)==0) | (len(connections)==0)):
                    norm_connections=[x/100 for x in connections]
                else:
                    norm_connections=[(x-min(connections))/(max(connections)-min(connections)) for x in connections]
            else:
                connections=norm_connections


            #Normalize sizes
            markers=[]
            for node in nxen.nodes(data=True):
                markers.append(node[1]['size'])
            if len(markers)!=0:
                
                if ((max(markers)-min(markers)==0) | (len(markers)==0)):
                    norm_markers=[x/100 for x in markers]
                else:
                    norm_markers=[(x-min(markers))/(max(markers)-min(markers)) for x in markers]
            else:
                markers=norm_markers
            
               
            norm_markers=np.clip(norm_markers,0.3, 1)
            
            
            fig,ax=plt.subplots(figsize=figsize)
            
            ##Plot the nodes
            xses,yses=[],[]
            lab=[]
            colors=[]
            for node in nxen.nodes(data=True):
                xses.append(pos_[node[0]][0])
                yses.append(pos_[node[0]][1])
                lab.append(node[0])
                colors.append(node[1]['pvalue'])
            nodez=ax.scatter(xses,yses,s=[mkfsize*size for size in norm_markers],
                           c=colors,cmap=colormap,vmax=maxpv,alpha=alpha,edgecolors=mkcolor,
                             linewidths=mklinewidths,clip_on=False,zorder=1)


            ##Mark the labels
            if labelling:
                labelling_without_overlapping(xses,yses,lab,ax,**kwargs)
            


            ##Plot the edges
            for indx, edge in enumerate(nxen.edges(data=True)):
                if edge[2]['weight'] > 0:
                    path_1 = edge[0]#prepare the data to insert in make edge
                    path_2 = edge[1]
                    x0, y0 = pos_[path_1]
                    x1, y1 = pos_[path_2]
                    edgez=ax.plot(np.linspace(x0,x1),np.linspace(y0,y1),
                            color=edgecolor,
                            linewidth = 3*norm_connections[indx]**4,
                                 zorder=0)


            cbar=plt.colorbar(nodez,ax=ax,orientation='horizontal',panchor=(0.5, 0.5))
            cbar.set_label(r'$-log_{10}(p-value)$',fontsize=cbarfontsize+4)
            cbar.ax.tick_params(labelsize=cbarfontsize)
            


            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.set_title(s,fontsize=20)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            axis = plt.gca()
            # maybe smaller factors work as well, but 1.1 works fine for this minimal example
            axis.set_xlim([factor*x for x in axis.get_xlim()])
            axis.set_ylim([factor*y for y in axis.get_ylim()])
            if savefig:
                plt.savefig(str(s)+'enrichment_analysis.jpeg', dpi=300)
            
            
            plt.tight_layout()
            plt.show()



        
        
             
        
    
    
    

# ### Plot global graph with legends 


def plot_enrichment_analisys_global(df,pvalue,legend=False,figsize=(40,10),title='',**kwargs):
######################################################################   
    def split_in_lines(string,n,words=False):
        if words:
            w = string.split()
            grouped_words = [' '.join(w[i: i + n]) for i in range(0, len(w), n)]
            return '\n'.join(grouped_words)
        grouped_chars=[string[i:i+n] for i in range(0, len(string), n)]
        return '\n'.join(grouped_chars)
####################################################################
    cmap = plt.get_cmap("tab10")
    colors = cmap(np.arange(len(df.source.value_counts())))
    fig, axes = plt.subplots(nrows=1, ncols=len(df.source.value_counts()), 
                             sharey=True,
                             #gridspec_kw={'width_ratios': np.log(df.source.value_counts().values)+5}, 
                             figsize=figsize)
    fig.subplots_adjust(wspace=0)
    a=0 
    shifted_xses=np.linspace(-2,2,len(df.source.value_counts()))
    for i,v in zip(df.source.value_counts().index,df.source.value_counts()):
        x=np.random.uniform(low=0.0, high=100, size=v)
        y=-np.log10(df[df.source==i].p_value.tolist())
        axes[a].scatter(x,y,color=colors[a],alpha=0.4,
                        s=np.clip(df[df.source==i].term_size.tolist(),150,500),label=y)
        axes[a].set_xlabel(i,fontsize=20)
        axes[a].spines['bottom'].set_color(colors[a])
        axes[a].spines['bottom'].set_linewidth(4)
        if a == 0:
            axes[a].spines['right'].set_visible(False)
            axes[a].spines['top'].set_visible(False)
            axes[a].xaxis.set_visible(True)
            axes[a].patch.set_alpha(0.5)
            axes[a].set_xticks([])

        else:
            axes[a].spines['right'].set_visible(False)
            axes[a].spines['top'].set_visible(False)
            axes[a].spines['left'].set_visible(False)
            axes[a].yaxis.set_visible(False)
            axes[a].set_xticks([])
            axes[a].patch.set_alpha(0)

        ##DEFINE HANDLES LABLES
        handles=[]
        lables=[]
        data_x=[]
        data_y=[]
        
        h=0
        for ind, pv in enumerate(y):
            if pv >= pvalue:
                handles.append(str(h))
                lables.append(split_in_lines(df[df.source==i].name.tolist()[ind],2,words=True))
                data_x.append(x[h])
                data_y.append(y[h])
                h+=1
        
        ##ADD LABLES
        if legend:
            labelling_without_overlapping(data_x,data_y,handles,axes[a],**kwargs)




            def create_proxy(label):
                line = matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='black',
                                               mec='none', marker=r'$\mathregular{{{}}}$'.format(label),markersize=20)
                return line



            proxies = [create_proxy(item) for item in handles]
            if len(proxies)>0:
                axes[a].legend(proxies,lables,bbox_to_anchor=(shifted_xses[a], 1.2), loc='lower left', borderaxespad=0,
                        ncol=1,fontsize=20,frameon=False,title=df.source.value_counts().index[a],
                               title_fontsize=20)

        a+=1
    fig.suptitle(title, fontsize="x-large")
    plt.show()
# # Plot enrichment analisys last



from utilities import average_list
import matplotlib
def plot_enrichment_analisys_1(df,p_value,figsize=(15,10),size_factor=80,save_fig=False):
    cmap = plt.get_cmap("tab10")
    colors=cmap(np.arange(len(df.source.value_counts())))
    #ylim=max(-np.log10(df.p_value.tolist()))
    for i, (s,v) in enumerate(zip(df.source.value_counts().index,df.source.value_counts())):
        y= list(reversed(-np.log10(list(filter(lambda x: -np.log10(x)>=p_value,df[df.source==s].p_value.tolist())))))
        x=list(reversed(df[(df.source==s) & (-np.log10(df.p_value)>p_value)].name.tolist()))
        sizes=[si for si in df[(df.source==s) & (-np.log10(df.p_value)>=p_value)].intersection_size.tolist()]
        if len(y)==0:
            continue
        else:
            fig,ax=plt.subplots(figsize=figsize)
            scat=ax.scatter(y,x,color=colors[i],alpha=0.4,s=[size_factor*s for s in sizes])
            ax.set_xlabel('$-\log_{10}(p-value)$',fontsize=15)
            ax.set_yticks([int(l) for l in ax.get_yticks()])
            ax.set_title(s,fontsize=20)
            

            ax.spines['bottom'].set_color(colors[i])
            ax.spines['left'].set_color(colors[i])
            ax.spines['bottom'].set_linewidth(4)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.grid(color=colors[i], linestyle='-', linewidth=2,alpha=0.05)
            ax.tick_params(direction='out',which= 'major' ,length=8, width=2, left=False,labelsize=20)
            
            handles, labels = scat.legend_elements(prop='sizes',alpha=0.6)
            lab=[int(re.search('\d+',l).group(0))//size_factor for l in labels] 
            
            legend_markers=(handles[0],handles[len(handles)//2],handles[-1])
            legend_labels=(lab[0],lab[len(handles)//2],lab[-1])
                        
            leg = ax.legend(legend_markers, legend_labels, 
                            loc="lower right", 
                            title="Number of Genes",
                            prop={'size': 30},
                            title_fontsize=20,
                           ncol=2,
                           markerscale=0.8)
            
            
            [hand.set_color(colors[i]) for hand in leg.legendHandles]
            if save_fig:
                title='enrich_%s'%(str(s))
                plt.savefig(title,dpi=400,bbox_inches='tight')
            plt.show()
           
            










