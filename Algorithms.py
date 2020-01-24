import networkx as nx
import numpy as np
#import cvxopt as cvx
#import picos as pic
import scipy.stats as stat
import scipy.sparse as sparse
import matplotlib.pyplot as plt
import operator
from collections import defaultdict
import itertools
import random
import math
import time
from heapq import nlargest
#import tensorflow as tf



def Binary_Quadratic(paramP, k):
    prob = pic.Problem()
    x = prob.add_variable('x', (nr_weak_edges, 1), vtype='binary')
#    x = prob.add_variable('x', (n, 1))
#    J = pic.new_param('J', np.ones((n,n)))
    e=pic.new_param('e', cvx.matrix(1,((nr_weak_edges,1))))
    P = pic.new_param('P', paramP.todense())

    prob.add_constraint(e|x <= k)
 #   prob.add_constraint(x.T*P*x <= up)
    prob.set_objective('max', x.T*P*x)
    prob.solve(verbose=0, solver='cplex')


    print('bound from the SDP relaxation: {0}'.format(prob.obj_value()))
    return x.value, prob.obj_value()

def greedy(Prow, Pcol, Pdiag, k, nr_weak_edges):
    selected=[]
    unselected=list(range(nr_weak_edges))
    cost=0
    while cost<=k:
        max_f=0
        max_j=-1
        for j in unselected:
            rowsum=0
            colsum=0
            if len(selected)>0:
#                rowsum=np.sum(Prow[j,:][:,selected])
#                colsum=np.sum(Pcol[:,j][selected,:])
                rowsum=np.sum(Prow[j,selected])
                colsum=np.sum(Pcol[selected,j])
            f=rowsum+colsum+Pdiag[j]
            if max_f < f:
                max_f=f
                max_j=j
        selected.append(max_j)
        unselected.remove(max_j)
        cost+=1
    x=np.zeros(nr_weak_edges)
    x[selected]=1
    return x

def randomized_greedy(P, w, k, x):
    selected=[]
    unselected=[]
    cost=0
    for j in range(nr_weak_edges):
        if x[j]==1:
            selected.append(j)
            cost+=w[j]
        else:
            unselected.append(j)
    
    obj_value=0
    iterations=1
    best_obj=0
    best_indices=[]
    for i in range(iterations):
#        print("=============ITERATION========")
        #fill up
        while cost<=k:
            max_f=0
            max_j=-1
            for j in range(nr_weak_edges):
                if j not in selected and cost+w[j]<=k:
                    rowsum=0
                    colsum=0
                    if len(selected)>0:
                        rowsum=np.sum(P[j,selected])
                        colsum=np.sum(P[selected,j])
                    f=(rowsum+colsum+P[j,j])/float(w[j])
                    if max_f < f:
                        max_f=f
                        max_j=j
            if max_f<=0:
                break
            obj_value+=max_f
            selected.append(max_j)
            cost+=w[max_j]
#            print("added: "+str(max_j))
#            print("the contribution is: "+str(max_f))

        #throw out
        while cost > k:
            min_f=10000000000
            min_j=-1
            for j in selected:
                rowsum=0
                colsum=0
                if len(selected)>0:
                    rowsum=np.sum(P[j,selected])
                    colsum=np.sum(P[selected,j])
                f=(rowsum+colsum-P[j,j])/float(w[j])
                if f< min_f:
                    min_f=f
                    min_j=j
            cost=cost-w[min_j]
            selected.remove(min_j)
            obj_value-=min_f
#            print("thrown out: "+str(min_j))
#            print("the contribution was: "+str(min_f))
            
        if obj_value > best_obj:
            best_obj=obj_value
            best_indices=list(selected)
#        print("objective value: "+str(obj_value))
        
        ind=random.randint(0, nr_weak_edges-1)
#        ind=random.choice(selected)
        if ind in selected:
            rowsum=np.sum(P[ind,selected])
            colsum=np.sum(P[selected,ind])
            f=(rowsum+colsum-P[ind,ind])/float(w[j])
            cost=cost-w[ind]
            selected.remove(ind)
            obj_value-=f
#        print("removed randomly: "+str(ind))
#        print("the contribution was: "+str(f))
        else:
            rowsum=np.sum(P[ind,selected])
            colsum=np.sum(P[selected,ind])
            f=(rowsum+colsum+P[ind,ind])/float(w[j])
            cost=cost+w[ind]
            selected.append(ind)
            obj_value+=f
#            print("added randomly: "+str(ind))
#            print("the contribution is: "+str(f))
            
    print(best_obj)
    x=np.zeros(nr_weak_edges)
    x[best_indices]=1
    return x        
            
                
                
def heuristic(gains, k, nr_weak_edges):
    f_best=0
    relative_gains=gains
#    start_time=time.time()
    max_indices=np.argpartition(relative_gains, -k)[-k:]
#    print("sorted the indices using argpartition in: "+str(time.time()-start_time))
#    start_time=time.time()
#    values,max_indices=tf.math.top_k(relative_gains, k=k, sorted=False, name=None)
#    print("ellapsed seconds for tensorflow: "+str(time.time()-start_time))
    x=np.zeros(nr_weak_edges)
    x[max_indices]=1
    
    return x


def quadratic(P, s):
    return np.transpose(s).dot(P.dot(s))



def find_wedges_fast(G, A):
    wedges=[]
    cou=0
    Asquared=A.dot(A).tocoo()
    sets_of_neighbours={}
    for i in range(A.shape[0]):
        sets_of_neighbours[i]=A[i,:]
    for x,y in zip(Asquared.row, Asquared.col):
        if not G.has_edge(x,y):
            neighborhood1=sets_of_neighbours[x]
            neighborhood2=sets_of_neighbours[y]
#            print(neighborhood1)
            prod=neighborhood1.multiply(neighborhood2)
            common_neighbors=np.nonzero(prod)[1]
            for n in common_neighbors:
                cou+=1
                tup1=tup(x,n)
                tup2=tup(n,y)
                wedges.append([tup1,tup2])
    return wedges
        
def find_wedges(A):#,c):
    wedges=[]
    cou=0
    for n in range(A.shape[0]):
        for (x,y) in itertools.combinations(np.nonzero(A[n,:])[1],2):
                if A[x,y]!=1:
                    cou+=1
#                    if c[n]==c[x] and c[x]==c[y]:
                    tup1=tup(x,n)
                    tup2=tup(n,y)
                    wedges.append([tup1,tup2])
    return wedges


def find_wedges_graph(G):#,c):
    wedges=[]
    cou=0
    for n in range(len(G.nodes())):
        for (x,y) in itertools.combinations(G.neighbors(n),2):
                if not G.has_edge(x,y):
                    cou+=1
#                    if c[n]==c[x] and c[x]==c[y]:
                    tup1=tup(x,n)
                    tup2=tup(n,y)
                    wedges.append([tup1,tup2])
    return wedges

def construct_P(G,labeling,ind):
    I=[]
    J=[]
    V=[]
    diagonal=np.array(m*[0])
    m=len(ind.keys())
    for n in range(len(G.nodes())):
        for (x,y) in itertools.combinations(G.neighbors(n),2):
            if not G.has_edge(x,y):
                e1=tup(x,n)
                e2=tup(n,y)
                if labeling[e1]==0 and labeling[e2]==0:
                    continue
                #e1 is weak
                if labeling[e1]==0 and labeling[e2]==1:
                    diagonal[ind[e1]]+=1
                #e2 is weak
                if labeling[e1]==1 and labeling[e2]==0:
                    diagonal[ind[e2]]+=1
    return diagonal


def tup(v1,v2):
    if int(v1)<int(v2):
        return (v1,v2)
    else:
        return (v2,v1)
    


def find_nr_triangles(G,A):
    tr={}
    
    A2=np.dot(A,A)
    for e in G.edges():
        v=int(G.node[e[0]]["old"])
        u=int(G.node[e[1]]["old"])
        tupl=tup(u,v)
        tr[tupl]=A2[e[0],e[1]]
           

    return tr

def find_nr_wedges(tr, G):
    we={}
    
    for e in G.edges():
        v=int(G.node[e[0]]["old"])
        u=int(G.node[e[1]]["old"])
        tupl=tup(u,v)
        we[tupl]=G.degree(e[0])+G.degree(e[1])-2*tr[tupl]
    
    return we


def read_graph(dataset_directory, graphfile):
    fh = open(dataset_directory + graphfile, 'rb')
    Gc = nx.read_edgelist(fh, delimiter=',', nodetype=int, data=(('weight', float),))
    O = max(nx.connected_component_subgraphs(Gc), key=len)
    print("Extracted largest connected component")
    O = nx.convert_node_labels_to_integers(O, first_label=0, label_attribute="real")
    return nx.Graph(O)



def read_communities(O, positive_file, negative_file):
    with open(positive_file) as f:
        positive_list = set(f.read().splitlines())
    with open(negative_file) as f:
        negative_list = set(f.read().splitlines())
    print("Positive nodes:")
    print(len(positive_list))
    print("Negative nodes:")
    print(len(negative_list))
    weights=[]

    for e in O.edges():
        e0=e[0]
        e1=e[1]
        if str(O.node[e0]["real"]) in positive_list and str(O.node[e1]["real"]) in positive_list:
            O[e0][e1]['weight']=0
            weights.append(0)
        elif str(O.node[e0]["real"]) in negative_list and str(O.node[e1]["real"]) in negative_list:
            O[e0][e1]['weight']=0
            weights.append(0)
        else:
            O[e0][e1]['weight']=1.0
            weights.append(1.0)
    return weights

def full_run(diag, nr_weak_edges, k):
    timestring=''
    objstring=''
#    print("##########################################################")
    
    start_time=time.time()
    x=heuristic(diag, k, nr_weak_edges)
    timestring=str(round(time.time()-start_time,3))+" \\\ "+timestring
#    print("ellapsed seconds for heuristic: "+str(time.time()-start_time))
    obj_value=quadratic(np.diag(diag),x)
#    print("violations: "+str(obj_value))
    objstring=str(int(obj_value))+" \\\ "+objstring
    
#    print("##########################################################")

    """
    start_time=time.time()
    x=greedy(Prow, Pcol, diag, k, nr_weak_edges)
#    x=greedy(P.todense(),w,k)
#    x=full_greedy(P.todense(), w, k, np.zeros(n))
    timestring=str(round(time.time()-start_time,3))+" & "+timestring
#    print("ellapsed seconds for greedy: "+str(time.time()-start_time))
    obj_value=quadratic(Prow,x)
#    print("violations: "+str(obj_value))
    objstring=str(int(obj_value))+" & "+objstring
    
    print("##########################################################")


    start_time=time.time()
    x,bound=Binary_Quadratic(P, k)
    timestring=str(round(time.time()-start_time,3))+" & "+timestring
    print("ellapsed seconds for binary Quadratic: "+str(time.time()-start_time))
    obj_value=quadratic(P,x)
    print("violations: "+str(obj_value))
    objstring=str(int(obj_value[0,0]))+" & "+objstring
    """

    
    


    print("objectives:")
    print(objstring)
    print("times:")
    print(timestring)

def read_edgelist_graph(dataset_directory, dataset):
    
    G = nx.Graph()
    weights=[]
    #n=len(G.nodes())
    edgelist=open(dataset_directory+dataset, 'r')
    for line in edgelist:
        splitted=line.strip().split(delimiter)
        v1=int(splitted[0])
        v2=int(splitted[1])
        tupl=tup(v1,v2)
        if len(splitted)>2:
            G.add_edge(tupl[0],tupl[1], weight=float(splitted[2]))
            weights.append(float(splitted[2]))
    
    return G,weights
def main_code(k):
    
    G,weights=read_edgelist_graph(dataset_directory,dataset)
#    G= read_graph(dataset_directory, 'edgesBig.txt')
    n=len(G)
#    weights=read_communities(G, dataset_directory+'Community1.txt', dataset_directory+'Community2.txt')
    G = max(nx.connected_component_subgraphs(G), key=len)
    G=nx.convert_node_labels_to_integers(G, label_attribute="old")


#    print('nodes: ' +str(len(G.nodes())))
#    print('edges: ' +str(len(G.edges())))
    splitting=np.percentile(weights, 70)
#    print(splitting)
    nr_strong_edges=0
#    global nr_weak_edges
    nr_weak_edges=0
    #dictionary that maps weak edges to indices in P    
    edge_to_ind={}  
    
    labeling={}
    for e in G.edges():
        if G[e[0]][e[1]]['weight']>splitting:
            #assign to strong
            labeling[e]=1
            nr_strong_edges+=1
        else:
            #assign to weak
            labeling[e]=0
            edge_to_ind[e]=nr_weak_edges
            nr_weak_edges+=1
#    print('strong edges: '+str(nr_strong_edges))
#    print('weak edges: '+str(nr_weak_edges))



#    A = nx.adjacency_matrix(G, weight=None)

    
    #tr=find_nr_triangles(G,A)
    #we=find_nr_wedges(tr, G)
    
#    start_time=time.time()
#    wedges=find_wedges(A)
#    print("ellapsed seconds for slow wedges: "+str(time.time()-start_time))
    start_time=time.time()
#    wedges=find_wedges_graph(G)
#    start_time=time.time()
#    wedges=find_wedges_fast(G, A)
#    print("ellapsed seconds for (supposedly) fast wedges: "+str(time.time()-start_time))
#    print("found wedges")
    diagonal=construct_P(G,labeling,edge_to_ind)
#    print("ellapsed seconds for P: "+str(time.time()-start_time))


    kappa=int(round(k*nr_weak_edges))
#    print("k is: "+str(kappa))
    full_run(diagonal, nr_weak_edges, k=kappa)


dataset='edgelist'

dataset_directory='/home/cloud-user/Datasets/stc/LesMis/'
delimiter=','
print("#######################LesMis###################")
main_code(0.1)
main_code(0.2)

#dataset_directory='/home/cloud-user/Datasets/stc/Bigdatasets/bitcoinotc/'
#delimiter=','
#print("#######################bitcoin otc###################")
#main_code()
#dataset_directory='/home/cloud-user/Datasets/stc/Bigdatasets/authors/'
#delimiter='\t'
#dataset_directory='/home/cloud-user/Datasets/stc/Bigdatasets/facebooksmall2/'
#delimiter=';'

dataset_directory='/home/cloud-user/Datasets/stc/Bigdatasets/KDD/'
delimiter=' '
print("#######################kdd###################")
main_code(0.01)
main_code(0.02)
main_code(0.1)

dataset_directory='/home/cloud-user/Datasets/stc/Bigdatasets/facebook/'
delimiter=';'
print("#######################facebook###################")
main_code(0.01)
main_code(0.02)
main_code(0.1)

dataset_directory='/home/cloud-user/Datasets/stc/Bigdatasets/twitter/'
delimiter=';'
print("#######################twitter###################")
main_code(0.01)
main_code(0.02)
main_code(0.1)
#dataset_directory='/home/cloud-user/Datasets/stc/Bigdatasets/Retweetssmall2/'
#delimiter=';'
dataset_directory='/home/cloud-user/Datasets/stc/Bigdatasets/Telecoms/'
delimiter=';'
print("#######################telecoms###################")
main_code(0.01)
main_code(0.02)
main_code(0.1)
dataset_directory='/home/cloud-user/Datasets/stc/Bigdatasets/bitcoinalpha/'
delimiter=','
print("#######################bitcoin alpha###################")
main_code(0.01)
main_code(0.02)
main_code(0.1)

#dataset_directory='/home/cloud-user/Datasets/TwitterFollowers/'



#main_code()


#assigned_c=open(dataset_directory+datafile+'_c.txt', 'r')






