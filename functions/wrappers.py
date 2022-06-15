#!/usr/bin/env python
# coding: utf-8

# # Maps identifiers UNIPROT
import pandas as pd
import requests
import numpy as np
import re
import os
import sys


ncbidb=pd.read_csv('https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz',sep='\t')
tab_uni=pd.read_csv('https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz', sep='\t',names=['uniprot','mapper','id'])


# # Wrapper to UNIPROT

def retrieve_protein_info(typology,ID,column):
    base_url='https://www.uniprot.org/uniprot/'
    payload={'query':'%s:%s AND reviewed:yes AND organism:9606' % (typology,ID),
            'format':'tab',
            'columns':'%s' % column}
    result= requests.get(base_url,params=payload)
    if result.ok:
        return result.text
    else:
        return 'something went wrong'


def gene_mapping(query,source,target):
    
    if source=='entrez' and target=='ensembl':
        output=[]
        for st in ncbidb[ncbidb.GeneID==query].dbXrefs.values[0].split('|'):
            if re.match('Ensembl',st):
                ensembl=re.sub(r'^.*?:', '', st)
                output.append(ensembl)
        return output
    elif source=='ensembl' and target=='entrez':
        try:
            return ncbidb[ncbidb.dbXrefs.str.contains(query)].GeneID.values[0]
        except:
            return 'none'
    elif source=='symbol' and target=='ensembl':
        for st in ncbidb[ncbidb.Symbol==query].dbXrefs.values[0].split('|'):
            if re.match('Ensembl',st):
                ensembl=re.sub(r'^.*?:', '', st)
                output.append(ensembl)
        return output
    elif source == 'ensembl' and target=='symbol':
        return ncbidb[ncbidb.dbXrefs.str.contains(query)].Symbol.values[0]
    elif source == 'entrez' and target=='symbol':
        return ncbidb[ncbidb.GeneID==query].Symbol.values[0]
    elif source=='symbol' and target=='entrez':
        return ncbidb[ncbidb.Symbol==query].GeneID.values[0]
        

def gene_mapping_many(query_list,source,target):

    if source=='ensembl' and target=='symbol':
        ense=ncbidb['dbXrefs'].apply(lambda x : re.sub(r'.+?(?=ENS)', '', x))
        dictio=dict((x.split('|')[0], y) for x, y in list(zip(ense,ncbidb['Symbol'])) if 'ENS' in x)
        return list(map(dictio.get,query_list))
    
    elif source=='ensembl' and target=='entrez':
        ense=ncbidb['dbXrefs'].apply(lambda x : re.sub(r'.+?(?=ENS)', '', x))
        dictio=dict((x.split('|')[0], y) for x, y in list(zip(ense,ncbidb['GeneID'])) if 'ENS' in x)
        return list(map(dictio.get,query_list))
    
    elif source=='ensembl' and target=='uniprot':
        dictio=dict(zip(tab_uni[tab_uni.mapper=='Ensembl'].id.tolist(),
                        tab_uni[tab_uni.mapper=='Ensembl'].uniprot.tolist()))
        return list(map(dictio.get,query_list))
    
    elif source=='entrez' and target=='ensembl':
        ense=ncbidb['dbXrefs'].apply(lambda x : re.sub(r'.+?(?=ENS)', '', x))
        dictio=dict((y, x.split('|')[0]) for x, y in list(zip(ense,ncbidb['GeneID'])) if 'ENS' in x)
        return list(map(dictio.get,query_list))
    
    elif source == 'entrez' and target == 'symbol':
        dictio=dict((x, y) for x, y in list(zip(ncbidb['GeneID'],ncbidb['Symbol'])))
        return list(map(dictio.get,query_list))
    
    elif source=='entrez' and target=='uniprot':
        dictio=dict(zip(tab_uni[tab_uni.mapper=='GeneID'].id.tolist(),
                        tab_uni[tab_uni.mapper=='GeneID'].uniprot.tolist()))
        return list(map(dictio.get,query_list))
    
    elif source=='symbol' and target=='ensembl':
        ense=ncbidb['dbXrefs'].apply(lambda x : re.sub(r'.+?(?=ENS)', '', x))
        dictio=dict((y, x.split('|')[0]) for x, y in list(zip(ense,ncbidb['Symbol'])) if 'ENS' in x)
        return list(map(dictio.get,query_list))
    
    elif source == 'symbol' and target == 'entrez':
        dictio=dict((y, x) for x, y in list(zip(ncbidb['GeneID'],ncbidb['Symbol'])))
        return list(map(dictio.get,query_list))
    
    elif source=='symbol' and target=='uniprot':
        dictio=dict(zip(tab_uni[tab_uni.mapper=='Gene_Name'].id.tolist(),
                        tab_uni[tab_uni.mapper=='Gene_Name'].uniprot.tolist()))
        return list(map(dictio.get,query_list))
    
    elif source=='uniprot' and target=='ensembl':
        dictio=dict(zip(tab_uni[tab_uni.mapper=='Ensembl'].uniprot.tolist(),
                        tab_uni[tab_uni.mapper=='Ensembl'].id.tolist()))
        return list(map(dictio.get,query_list))
    
    elif source=='uniprot' and target=='entrez':
        dictio=dict(zip(tab_uni[tab_uni.mapper=='GeneID'].uniprot.tolist(),
                        tab_uni[tab_uni.mapper=='GeneID'].id.tolist()))
        return list(map(dictio.get,query_list))
    
    elif source=='uniprot' and target=='symbol':
        dictio=dict(zip(tab_uni[tab_uni.mapper=='Gene_Name'].uniprot.tolist(),
                        tab_uni[tab_uni.mapper=='Gene_Name'].id.tolist()))
        return list(map(dictio.get,query_list))


# # RETRIEVE GENE DISEASE ASSOCIATIONS AND VARIANT DISEASE ASSOCIATIONS



def retrieve_gda(disease):
    auth_params = {"email":"francesco.gualdi@live.com","password":"nidodelcuculo0"}
    api_host='https://www.disgenet.org/api'
    req=requests.Session()
    url = api_host+'/auth/'
    response = req.post(url, data=auth_params)
    token=response.json()['token']
    req.headers.update({"Authorization": "Bearer %s" % token}) 
    return req.get(api_host+'/gda/disease/{}'.format(disease)).json()
    
def retrieve_vda(disease):
    auth_params = {"email":"francesco.gualdi@live.com","password":"nidodelcuculo0"}
    api_host='https://www.disgenet.org/api'
    req=requests.Session()
    url = api_host+'/auth/'
    response = req.post(url, data=auth_params)
    token=response.json()['token']
    req.headers.update({"Authorization": "Bearer %s" % token}) 
    return req.get(api_host+'/vda/disease/{}'.format(disease)).json()

# # PANTHER GENE CLASSIFICATION # #

def retrieve_from_panther(gene):
    url='http://pantherdb.org/services/oai/pantherdb/geneinfo?geneInputList=%s&organism=9606'%(gene)
    response=requests.get(url)
    return response.json()

def retrieve_from_panther_many(gene_list):
    url='http://pantherdb.org/services/oai/pantherdb/geneinfo?geneInputList=%s&organism=9606'%('%2C'.join(gene_list))
    response=requests.get(url)
    results=[]
    i=0
    for gene in response.json()['search']['mapped_genes']['gene']:
        try:
            if type(gene['annotation_type_list']['annotation_data_type'][0]['annotation_list']['annotation'])==list:
                tmp=[]
                for el in gene['annotation_type_list']['annotation_data_type'][0]['annotation_list']['annotation']:
                    tmp.append(el['name'])
                results.append(tuple(tmp))
            else:
                results.append(gene['annotation_type_list']['annotation_data_type'][0]['annotation_list']['annotation']['name'])
            i+=1
        except:
            results.append(None)
            print(gene_list[i])
            i+=1
    return results
            

    
    

class gwas_api():
    
    ##Retrieves associations to an efo trait
    def get_associations(efotrait):
        df=pd.DataFrame(columns=['variantid','p-value','risk_allele','RAF','beta','CI'])
        http= 'https://www.ebi.ac.uk/gwas/rest/api/efoTraits/%s/associations' %(efotrait)
        associ=requests.get(http).json()
        for i,element in enumerate(associ['_embedded']['associations']):
            try:
                df.at[i,'variantid']=''.join(element['loci'][0]['strongestRiskAlleles'][0]['riskAlleleName'].split('-')[0:1])
                df.at[i,'risk_allele']=element['loci'][0]['strongestRiskAlleles'][0]['riskAlleleName'].split('-')[-1]
                df.at[i,'p-value']=int(element['pvalueMantissa'])*10**int(element['pvalueExponent'])
                df.at[i,'RAF']=float(element['riskFrequency'])
                df.at[i,'beta']=[float(element['betaNum']) if type(element['betaNum'])==float else None][0]
                df.at[i,'SE']=float(element['standardError'])
                df.at[i,'CI']=element['range']
            except:
                pass
        return df
    
    #Retrieves position of a given variant
    def get_variant_position(variant):
        http= "https://rest.ensembl.org/variation/human/%s" %(variant)
        try:
            associ=requests.get(http,headers={ "Content-Type" : "application/json"}).json()
            pos=associ['mappings'][0]['location'].split(':')[0],associ['mappings'][0]['location'].split(':')[1].split('-')[0]
            return pos
        except:
            pass
    
    
    #retrieve the coordinates of many genes
    def get_gene_position_many(idlist,chunked=False,chunksize=200):
        http="https://rest.ensembl.org/lookup/id"
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
        if chunked:
            chunked_idlist=[]
            print('total number of chunks: %s' %(int(len(idlist)/chunksize)+1))
            for i in range(0,len(idlist),chunksize):
                chunked_idlist.append(idlist[i:i+chunksize])
            results=[]
            for i,chunk in enumerate(chunked_idlist):
                response = requests.post(http, headers=headers, 
                                         data="{" + '"ids" : {}'.format(str(chunk).replace("'",'"'))+"}").json()

                ListOfTuples=[]
                for k,v in response.items():
                    try:
                        ListOfTuples.append((k,v['seq_region_name'],v['start'],v['end']))

                    except:

                        print('error') 
                        continue

                results.append(ListOfTuples)
                print('chunk %s processed' % (i))
            return sum(results,[])
        else:
            response = requests.post(http, headers=headers, 
                                         data="{" + '"ids" : {}'.format(str(idlist).replace("'",'"'))+"}").json()

            ListOfTuples=[]

            for k,v in response.items():
                try:
                    ListOfTuples.append((k,v['seq_region_name'],v['start'],v['end']))
                except:

                    print('error')
                    pass
            return ListOfTuples

    def gene_coordinates(gene_name,Id='symbol'):
        http="https://rest.ensembl.org//lookup/symbol/homo_sapiens/%s?expand=1"%(gene_name)
        risposta=requests.get(http,headers={ "Content-Type" : "application/json"}).json()
        return risposta['start'],risposta['end']
    #Retrieve position for many variants
    def get_variant_position_many(idlist,chunked=False,chunksize=200):
        http="https://rest.ensembl.org/variation/homo_sapiens"
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
        if chunked:
            chunked_idlist=[]
            print('total number of chunks: %s' %(int(len(idlist)/chunksize)+1))
            for i in range(0,len(idlist),chunksize):
                chunked_idlist.append(idlist[i:i+chunksize])
            results=[]
            for i,chunk in enumerate(chunked_idlist):
                response = requests.post(http, headers=headers, 
                                         data="{" + '"ids" : {}'.format(str(chunk).replace("'",'"'))+"}").json()
                results.append(list(zip(response.keys(),
                            [Map['mappings'][0]['location'].split(':')[0] for Map in response.values()],
                           [Map['mappings'][0]['start'] for Map in response.values()])))
                print('chunk %s processed' % (i))
            return sum(results,[])
        else:
            response = requests.post(http, headers=headers, 
                                         data="{" + '"ids" : {}'.format(str(idlist).replace("'",'"'))+"}").json()


            return list(zip(response.keys(),
                            [Map['mappings'][0]['location'].split(':')[0] for Map in response.values()],
                           [Map['mappings'][0]['start'] for Map in response.values()]))
    
    
    def VEP_for_rsid(idlist, chunked=False,chunksize=200):
        '''Variants must be fed in HGVS notation that is cr:glocREF>ALT'''
        http="https://rest.ensembl.org/vep/human/hgvs"
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
        chunked_idlist=[]
        if chunked:
            print('total number of chunks: %s' %(int(len(idlist)/chunksize)+1))
            for i in range(0,len(idlist),chunksize):
                chunked_idlist.append(idlist[i:i+chunksize])
            results=[]
            for i,chunk in enumerate(chunked_idlist):
                response = requests.post(http, headers=headers, 
                                         data="{" + '"hgvs_notations" : {}'.format(str(chunk).replace("'",'"'))+"}")
                results.append(response.json())
                print('chunk %s processed' % (i))
            return results

        else:
            response = requests.post(http, headers=headers, 
                                         data="{" + '"hgvs_notations" : {}'.format(str(idlist).replace("'",'"'))+"}")
            return response.json()


    #Retrieves variant in LD with a given variant
    def get_variants_in_LD(variant,r2,pop='EUR'):
        http= "https://rest.ensembl.org/ld/human/%s/1000GENOMES:phase_3:%s?r2=%s" %(variant,pop,r2)
        try:
            variants=requests.get(http,headers={ "Content-Type" : "application/json"}).json()
      
            return [x['variation2'] for x in variants if float(x['r2'])>=r2]
        except:
            pass

    #Retrieve summary statistic of a given study    
    def get_summary_statistic(study):
        http= 'https://www.ebi.ac.uk/gwas/summary-statistics/api/studies/%s/associations' %(study)
        ss=requests.get(http).json()
        return ss
        
    
    #Retrieves the list of studies with summary statistics available for a given trait
    def get_summary_statistic_list():
        http= 'https://www.ebi.ac.uk/gwas/summary-statistics/api/associations'
        ss=requests.get(http).json()
        return ss
    
    #Retrieves annotations of a given genomic region
    def get_phenotypes(chromosome,start,stop,feature_type='Genes',only_phenotypes=1):
        http="https://rest.ensembl.org/phenotype/region/homo_sapiens/%s:%s-%s?only_phenotypes=%s;feature_type=%s"%(chromosome,start,stop,only_phenotypes,feature_type)
        annot=requests.get(http, headers={ "Content-Type" : "application/json"})
        return annot.json()
    
    #Retrieves overlapping elements of a given region 
    def get_ov_region(chromosome,start,stop,features=list):
        str_features=';'.join(['feature='+x for x in features])
        http="https://rest.ensembl.org/overlap/region/human/%s:%s-%s?%s"%(chromosome,start,stop,str_features)
        risposta=requests.get(http,headers={ "Content-Type" : "application/json"}).json()
        return risposta
    
    #Retrieves the nucleotide sequence of a given position 
    def get_sequence(chromosome,start,stop):
        http="https://rest.ensembl.org/sequence/region/human/%s:%s..%s?" %(chromosome,start,stop)
        risposta=requests.get(http,headers={ "Content-Type" : "application/json"}).json()
        return risposta['seq']
    
    ## lift from grch37 to 38
    def grch_liftover(chromosome,start,end):
        url="https://rest.ensembl.org/map/human/GRCh37/%s:%i..%i/GRCh38?"%(chromosome,start,end)
        r = requests.get(url, headers={ "Content-Type" : "application/json"}).json()
        try:
            return (chromosome,r['mappings'][0]['mapped']['start'],r['mappings'][0]['mapped']['end'])
        except:
            return None
    
    ## function to get variant list of eqtls
    def get_eqtl_variant(rsid):
        url='http://www.ebi.ac.uk/eqtl/api/associations/%s'%(rsid)
        risp=requests.get(url).json()
        return risp
    ## return single variant info ##
    def get_variant_info(variant):
        http= "https://rest.ensembl.org/variation/human/%s?"%(variant)
        associ=requests.get(http,headers={ "Content-Type" : "application/json"}).json()
        return associ
    def get_eqtl_df(rsid,p_value=0.005,increase_index=False):
        url='http://www.ebi.ac.uk/eqtl/api/associations/%s?size=1000'%(rsid)
        eqtls=requests.get(url).json()
        genes_eqtls=[]
        eqtl_df=pd.DataFrame(columns=['variantid','p_value','gene'])
        for ass in eqtls['_embedded']['associations'].keys():
            eqtl_df.loc[ass]=[rsid,eqtls['_embedded']['associations'][ass]['pvalue'],eqtls['_embedded']['associations'][ass]['gene_id']]
        eqtl_df=eqtl_df.loc[eqtl_df.p_value<=p_value]
        eqtl_df['-log10_pval']=eqtl_df.p_value.apply(lambda x: -np.log10(x))
        eqtl_df=eqtl_df.reset_index(drop=True)
        if increase_index:
            eqtl_df.index+=1
        return eqtl_df



    
    ##get genes in a window centered in a genomic position and compute the distance between the position and all the genes
    def get_genes(cr,start,stop,window=10000,pop='EUR',features=['gene'],mode='all'):
        winstart=start-window//2
        winend=stop+window//2
        str_features=';'.join(['feature='+x for x in features])
        http="https://rest.ensembl.org/overlap/region/human/%s:%s-%s?%s"%(cr,winstart,winend,str_features)
        risposta=requests.get(http,headers={ "Content-Type" : "application/json"}).json()
        if mode=='complete_data':
            return risposta
        elif mode=='all':
            elements={}
            for el in risposta:
                if el['biotype']=='protein_coding':
                    try:
                        elements[el['external_name']]=int(el['start']-start)
                    except:
                        pass
            return elements
        elif mode=='closest_forward':
            elements={}
            for el in risposta:
                if el['biotype']=='protein_coding':
                    try:
                        elements[el['external_name']]=int(el['start']-start)
                    except:
                        pass
            try:
                return min([(k,v) for (k,v) in elements.items() if v>0], key=lambda x:x[1])
            except:
                return 'no_genes_forward'
        elif mode=='closest_backward':
            elements={}
            for el in risposta:
                if el['biotype']=='protein_coding':
                    try:
                        elements[el['external_name']]=int(el['start']-start)
                    except:
                        pass
            try:
                return max([(k,v) for (k,v) in elements.items() if v<0], key=lambda x:x[1])
            except:
                return 'no_genes_backward'
        elif mode=='closest_overall':
            elements={}
            for el in risposta:
                if el['biotype']=='protein_coding':
                    try:
                        elements[el['external_name']]=int(el['start']-start)
                    except:
                        pass
            try:
                return min([(k,np.absolute(v)) for (k,v) in elements.items()], key=lambda x:x[1])
            except:
                return 'no_genes'
    
   
    
    
    def VEP(idlist,chunked=False,all_data=False,output_type='gene_symbol'):
        http="https://rest.ensembl.org/vep/human/id/"
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
        sev_var_gen=[]
        nsev_var_gen=[]
        severe_effects=['transcript_ablation','splice_acceptor_variant','splice_donor_variant','stop_gained',
                            'frameshift_variant','stop_lost','start_lost','transcript_amplification','missense_variant']
        reg_consequence=['protein_altering_variant','coding_sequence_variant',
                               'TFBS_ablation','TFBS_amplification','TF_binding_site_variant','regulatory_region_ablation'
                                 'regulatory_region_amplification','regulatory_region_variant','feature_truncation']


        non_severe_effects=['5_prime_UTR_variant','3_prime_UTR_variant','inframe_insertion','inframe_deletion','intron_variant','upstream_gene_variant','downstream_gene_variant','non_coding_transcript_exon_variant']

        if len(idlist)<=200:
            response = requests.post(http, headers=headers, data="{" + '"ids" : {}'.format(str(idlist).replace("'",'"'))+"}")
            if response.ok:
                r=response.json()
                if all_data:
                    return r
                else:
                    ## Severe Variants
                    for e in r:
                        if e['most_severe_consequence'] in severe_effects:
                            for transcript_consequences in e['transcript_consequences']:
                                if any(item in transcript_consequences['consequence_terms'] for item in severe_effects):
                                    if output_type=='gene_symbol':
                                        try:
                                            tupla=(e['id'],transcript_consequences['gene_symbol'])
                                            if tupla not in sev_var_gen:
                                                sev_var_gen.append(tupla)
                                        except:
                                            tupla=(e['id'],transcript_consequences['gene_id'])
                                            if tupla not in sev_var_gen:
                                                sev_var_gen.append(tupla)
                                    else:
                                        tupla=(e['id'],transcript_consequences['gene_id'])
                                        if tupla not in sev_var_gen:
                                            sev_var_gen.append(tupla)
                                        
                    ## Non Severe Effects
                        elif e['most_severe_consequence'] in non_severe_effects:
                            nosev=[]
                            for tc in e['transcript_consequences']:
                                #check wether consequence term is present in non_severe_effect
                                if ((any(x in tc['consequence_terms'] for x in non_severe_effects)) & (tc['biotype']=='protein_coding')):
                                    if output_type=='gene_symbol':
                                        try:
                                            nosev.append(tc['gene_symbol'])
                                        except: 
                                            nosev.append(tc['gene_id'])
                                    else:
                                        nosev.append(tc['gene_id'])

                            nsev_var_gen.append((e['id'],list(set(nosev))))

                    ### Transcription factors
                    
                    try:
                        tf={x['id']:[(','.join(ele['transcription_factors']).replace('::',':'),ele['motif_score_change']) for ele in x['motif_feature_consequences']] for x in r if x['most_severe_consequence']=='TF_binding_site_variant'}
                    except:
                        tf={x['id']:[(','.join(ele['transcription_factors']).replace('::',':'),0) for ele in x['motif_feature_consequences']] for x in r if x['most_severe_consequence']=='TF_binding_site_variant'}




                     ## Return all the effects ##
                    return sev_var_gen,nsev_var_gen,tf


            ## Return the response if there was an error in the query ##
            else:
                ##Code for Chunked case##
                print(response.text)
        elif (len(idlist)>200 | chunked):

            chunked_idlist=[]
            for i in range(0,len(idlist),200):
                chunked_idlist.append(idlist[i:i+200])
            for variant_lst in chunked_idlist:
                tmp_sev_var_gen=[]
                tmp_nsev_var_gen=[]
                tf={}
                response = requests.post(http, headers=headers, data="{" + '"ids" : {}'.format(str(variant_lst).replace("'",'"'))+"}")
                if response.ok:
                    r=response.json()
                    if all_data:
                        return r
                    else:
                        ## Severe Variants
                        for e in r:
                            if e['most_severe_consequence'] in severe_effects:
                                for transcript_consequences in e['transcript_consequences']:
                                    if any(item in transcript_consequences['consequence_terms'] for item in severe_effects):
                                        if output_type=='gene_symbol':
                                            try:
                                                tupla=(e['id'],transcript_consequences['gene_symbol'])
                                                if tupla not in tmp_sev_var_gen:
                                                    tmp_sev_var_gen.append(tupla)
                                            except:
                                                tupla=(e['id'],transcript_consequences['gene_id'])
                                                if tupla not in tmp_sev_var_gen:
                                                    tmp_sev_var_gen.append(tupla)
                                        else:
                                            tupla=(e['id'],transcript_consequences['gene_id'])
                                            if tupla not in tmp_sev_var_gen:
                                                tmp_sev_var_gen.append(tupla)
                                            
                        ## Non Severe Effects
                            elif e['most_severe_consequence'] in non_severe_effects:
                                nosev=[]
                                for tc in e['transcript_consequences']:
                                    #check wether consequence term is present in non_severe_effect
                                    if ((any(x in tc['consequence_terms'] for x in non_severe_effects)) & (tc['biotype']=='protein_coding')):
                                        if output_type=='gene_symbol':
                                            try:
                                                nosev.append(tc['gene_symbol'])
                                            except: 
                                                nosev.append(tc['gene_id'])
                                        else:
                                            nosev.append(tc['gene_id'])
                                            

                                tmp_nsev_var_gen.append((e['id'],list(set(nosev))))

                        sev_var_gen=sev_var_gen+tmp_sev_var_gen
                        nsev_var_gen=nsev_var_gen+tmp_nsev_var_gen
                        ### Transcription factors
                        try:
                            tf.update({x['id']:[(','.join(ele['transcription_factors']).replace('::',':'),ele['motif_score_change']) for ele in x['motif_feature_consequences']] for x in r if x['most_severe_consequence']=='TF_binding_site_variant'})
                        except:
                            tf.update({x['id']:[(','.join(ele['transcription_factors']).replace('::',':'),0) for ele in x['motif_feature_consequences']] for x in r if x['most_severe_consequence']=='TF_binding_site_variant'})


                    ## Return all the effects ##
                    return sev_var_gen,nsev_var_gen,tf

            else:
                print(response.text)



def VEP_for_rsid(idlist, chunked=False,chunksize=200):
    http="https://rest.ensembl.org/vep/human/hgvs"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    chunked_idlist=[]
    if chunked:
        print('total number of chunks: %s' %(int(len(idlist)/chunksize)+1))
        for i in range(0,len(idlist),chunksize):
            chunked_idlist.append(idlist[i:i+chunksize])
        results=[]
        for i,chunk in enumerate(chunked_idlist):
            response = requests.post(http, headers=headers, 
                                     data="{" + '"hgvs_notations" : {}'.format(str(chunk).replace("'",'"'))+"}")
            results.append(response.json())
            print('chunk %s processed' % (i))
        return results
            
    else:
        response = requests.post(http, headers=headers, 
                                     data="{" + '"hgvs_notations" : {}'.format(str(idlist).replace("'",'"'))+"}")
        return response.json()

    
 ###################################################################################################################

                                                #SQL# 

####################################################################################################################

import sqlite3

class sql:
    def list_tables(con):
        cur=con.cursor()
        cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
        return cur.fetchall()

    def display_table(table,con):
        cursor=con.cursor()
        cursor.execute('SELECT * FROM {}'.format(table))
        index= [i[0] for i in cursor.description]
        dataframe=pd.DataFrame(cursor.fetchall(),columns=index)
        return dataframe.set_index(dataframe.columns[0])

    def drop_table(table,con):
        cursor=con.cursor()
        cursor.execute('DROP TABLE {}'.format(table))
        con.commit()
    
    def drop_all_tables(con):
        for c in sql.list_tables(con):
            sql.drop_table(c[0],con)


    def executeScriptsFromFile(filepath, connessione):
        fd = open(filepath, 'r')
        sql = fd.read()
        fd.close()
        one=''.join(sql.split('\n'))
        reso=[x +';' for x in one.split(';')]
        reso.pop(-1)
        for i,x in enumerate(reso):
            try:
                cursor=connessione.cursor()
                cursor.execute(x)
                cursor.close()
            except:
                return(i,x)






