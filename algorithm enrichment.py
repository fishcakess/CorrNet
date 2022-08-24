import scipy.stats as stats
import json
import pickle
import copy

#create pathway list
mypathways = []

with open('Metabolism.txt') as file:
    content = file.readlines()[1:]
    meta = []
    for i in content:
        mypathways.append(i.strip())
        meta.append(i.strip())
    file.close()
    
with open('Cellular Processes.txt') as file:
    content = file.readlines()[1:]
    cp = []
    for i in content:
        mypathways.append(i.strip())
        cp.append(i.strip())
    file.close()
    
with open('Environmental Information Processing.txt') as file:
    content = file.readlines()[1:]
    eip = []
    for i in content:
        mypathways.append(i.strip())
        eip.append(i.strip())
    file.close()
    
with open('Genetic Information Processing.txt') as file:
    content = file.readlines()[1:]
    gip = []
    for i in content:
        mypathways.append(i.strip())
        gip.append(i.strip())
    file.close()    

def load_clusters(filepath):
    clusters = {}
    with open(filepath) as file:
        content = file.readlines()
        for i in content:
            i = i.strip()
            i = i.split(',')
            pathway = i[0].strip('"')
            cluster = i[1]
            if cluster not in clusters:
                clusters[cluster] = []
            clusters[cluster].append(pathway) 
    return clusters
            
#network clusters come as dictionary - think I dont need the other function to tidy names 

#make a function to get counts per cluster of each pathway type 

def pathway_count(clustersinput):
    
    pathtypes = [meta,cp,eip,gip]
    
    result = {}
    
    lenmeta = len(meta)
    lencp = len(cp)
    leneip =len(eip)
    lengip = len(gip)
    
    for clustersgroups,pathwaysgroups in clustersinput.items():
        lencluster = len(pathwaysgroups)

        metcount = 0
        cpcount = 0
        eipcount = 0 
        gipcount = 0

        for pathways in pathwaysgroups:
            if pathways in meta:
                metcount+=1
            if pathways in cp:
                cpcount+=1
            if pathways in eip:
                eipcount+=1
            if pathways in gip:
                gipcount+=1
        #got all counts now calculate
        result[clustersgroups] = [metcount,cpcount,eipcount,gipcount,(metcount+cpcount+eipcount+gipcount)]
    print('Met,',lenmeta,'CP,',lencp,'EIP,',leneip,'GIP,',lengip,'Total',144)
    return result

### enrichment scores

def pathway_enrichment_all(clusters):
    
    pathtypes = [meta,cp,eip,gip]
    
    result = {}
    
    lenmeta = len(meta)
    lencp = len(cp)
    leneip =len(eip)
    lengip = len(gip)
    
    sortedkeys = sorted(clusters.keys())
    
    
    for key in sortedkeys:
        
        pathwaysgroups = clusters[key]
        lencluster = len(pathwaysgroups)

            
        metcount = 0
        cpcount = 0
        eipcount = 0 
        gipcount = 0

        for pathways in pathwaysgroups:
            if pathways in meta:
                metcount+=1
            if pathways in cp:
                cpcount+=1
            if pathways in eip:
                eipcount+=1
            if pathways in gip:
                gipcount+=1
        #got all counts now calculate
        meta_r = (metcount/lencluster)/(lenmeta/144)
        cp_r = (cpcount/lencluster)/(lencp/144)
        eip_r = (eipcount/lencluster)/(leneip/144)
        gip_r = (gipcount/lencluster)/(lengip/144)
        result[key] = [meta_r,cp_r,eip_r,gip_r]

    print('Meta, cp, eip, gip') 
    return result

def pathway_enrichment_nometa(clusters):
    
    pathtypes = [cp,eip,gip]
    
    result = {}
    lencp = len(cp)
    leneip =len(eip)
    lengip = len(gip)
    
    sortedkeys = sorted(clusters.keys())
    
    
    for key in sortedkeys:
        
        pathwaysgroups = clusters[key]
        lencluster = len(pathwaysgroups)

            
        cpcount = 0
        eipcount = 0 
        gipcount = 0

        for pathways in pathwaysgroups:
            if pathways in cp:
                cpcount+=1
            if pathways in eip:
                eipcount+=1
            if pathways in gip:
                gipcount+=1
        #got all counts now calculate
        cp_r = (cpcount/lencluster)/(lencp/144)
        eip_r = (eipcount/lencluster)/(leneip/144)
        gip_r = (gipcount/lencluster)/(lengip/144)
        result[key] = [cp_r,eip_r,gip_r]

    print('cp, eip, gip') 
    return result

def pathway_enrichment_meta_eip(clusters):
    
    pathtypes = [meta,eip]
    
    result = {}
    
    lenmeta = len(meta)
    leneip =len(eip)
    sortedkeys = sorted(clusters.keys())
    
    
    for key in sortedkeys:
        
        pathwaysgroups = clusters[key]
        lencluster = len(pathwaysgroups)

            
        metcount = 0
        eipcount = 0 

        for pathways in pathwaysgroups:
            if pathways in meta:
                metcount+=1

            if pathways in eip:
                eipcount+=1

        #got all counts now calculate
        meta_r = (metcount/lencluster)/(lenmeta/144)
        eip_r = (eipcount/lencluster)/(leneip/144)
        result[key] = [meta_r,eip_r]

    print('Meta, eip') 
    return result
#enrichment signinifance 

def distribution(clustersinput):
    distribution = {'Met':[],'CP':[],'EIP':[],'GIP':[]}
    sortedkeys = sorted(clustersinput.keys())
    #now the output is cluster 1, 2, 3, 4

    for key in sortedkeys:
        clusterlist = clustersinput[key]

        for classification in distribution.keys():

            M = 144
            N = len(clusterlist)
            if classification == 'Met':
                index = 0
                n = len(meta)
            elif classification == 'CP': 
                index = 1
                n = len(cp)
            elif classification == 'EIP':
                index = 2
                n = len(eip)
            elif classification == 'GIP':
                index = 3
                n = len(gip)

            x = clusterlist[index]
            pval = stats.hypergeom(M = M, n = n, N = N).sf(x-1)
            distribution[classification].append(pval)

    return distribution


#add check to remove distribution key if all  == 1 that means that the classification was removed from corrmap


def clusters_enrichment_scores(clustersfile):
            
    with open(clustersfile) as f:
        clusters = json.load(f)

    pathwaycount = pathway_count(clusters)
    pathwayenrichment = pathway_enrichment_nometa(clusters)
    pathwaydistribution = distribution(pathwaycount)

    return pathwayenrichment, pathwaydistribution
    

pathwayenrichment, pathwaydistribution = clusters_enrichment_scores('no_met_hca.json')

iterable = copy.deepcopy(pathwaydistribution)
for key,value in iterable.items():
    x = pathwaydistribution[key]
    if x.count(x[0]) == len(x):
        del pathwaydistribution[key] 

with open('no_met_hca_pathwayenrichment','w') as f:
    json.dump(pathwayenrichment,f)

with open('no_met_hca_pathwaydistribution','w') as f:
    json.dump(pathwaydistribution,f)