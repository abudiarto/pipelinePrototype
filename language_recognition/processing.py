import pandas as pd
import requests, sys
from pandas.io.json import json_normalize

import pandas as pd
import numpy as np

from pyvis.network import Network

from bs4 import BeautifulSoup

def countInput(inputData):
	data = inputData.splitlines()
	# data = pd.DataFrame(data)
	# data.columns = ['RSID']
	return(data)


def vepQuery(inputData):
	import requests, sys
	from pandas.io.json import json_normalize
	import pandas as pd


	snps = countInput(inputData)
	myobj = {"ids": snps}
	myobj = str(myobj).replace("'", "/")
	myobj = myobj.replace("/", '"')
	    
	server = "http://grch37.rest.ensembl.org"
	ext = "/vep/human/id"
	headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
	r = requests.post(server+ext, headers=headers, data=str(myobj))
	 
	if not r.ok:
	  r.raise_for_status()
	  sys.exit()
	 
	decoded = r.json()
	
	decoded = [x for x in decoded if x['most_severe_consequence'] == 'intron_variant']

	fields = ['id', 'seq_region_name', 'start', 'end', 'biotype', 'consequence_terms', 'gene_symbol', 'gene_id', ] 
	SNP_list = pd.DataFrame(columns=fields)

	if len(decoded) > 0:
		SNP_list = json_normalize(data=decoded, record_path='transcript_consequences', 
	                            meta=['id', 'start','end', 'allele_string','seq_region_name'],
	                            errors='ignore')
		SNP_list = SNP_list.drop_duplicates(subset = ['gene_id', 'gene_symbol'], keep = 'first')
		
		SNP_list = SNP_list[fields]
		SNP_list.reset_index(inplace = True, drop = True)
		SNP_list["consequence_terms"] = SNP_list.apply(lambda x: x["consequence_terms"][0], axis=1)
		SNP_GENE = SNP_list[["id","gene_symbol", "gene_id", "consequence_terms"]]
	
		return(SNP_list.to_html(), list(SNP_list.gene_symbol), SNP_GENE)
	else:
		error = '<h3> All input variants are in intergenic regions </h3>'
		return(error, error, 0)

def wikiPathwayQuery(listGene):
	separator = '&ids='
	ids = (separator.join(listGene))
	APIendPoint = 'https://webservice.wikipathways.org/findPathwaysByXref?ids='
	output = '&format=json'
	requestString = APIendPoint + ids + output
	reqWiki = requests.get(requestString)
	genePathway = reqWiki.json()

	pathwayResult = json_normalize(data=genePathway['result'])
	field = ["fields.x.id.values", "id", "name", "score.0", "species","revision","url"]
	pathwayResult = pathwayResult[field]
	pathwayResult = pathwayResult.drop_duplicates("id")
	pathwayResult = pathwayResult.reset_index(drop=True)
	pathwayResult.columns = ["GeneEnsembleID", "PathwayID", "PathwayName", "PathwayScore", "Species", "Revision", "URL"]
	pathwayResult["GeneEnsembleID"] = pathwayResult.apply(lambda x: x["GeneEnsembleID"][0], axis=1)
	pathwayResult = pathwayResult[pathwayResult.Species=="Homo sapiens"]
	GENE_PATHWAY = pathwayResult[["GeneEnsembleID","PathwayID","PathwayName"]]


	return(pathwayResult.to_html(), GENE_PATHWAY)

def buildNetwork(GeneSNP, GenePathway):
	GeneSNP = GeneSNP
	GenePathway = GenePathway
	a = GeneSNP[["gene_symbol", "id", "id","consequence_terms"]]
	a.columns = ["var1","var2", "var3", "var4"]
	b = GenePathway[["GeneEnsembleID", "PathwayID", "PathwayName"]]
	b.columns = ["var1","var2","var3"]
	b.var4 = 99
	dataInput = pd.concat([a,b],axis=0)
	dataInput["value"] = .5


	got_net = Network(height="720", width="1280", bgcolor="#FFFFFF", font_color="black")

	# set the physics layout of the network
	#got_net.barnes_hut(central_gravity=0.8)
	got_net.force_atlas_2based(gravity=-100, central_gravity=0.01, spring_length=100,
	spring_strength=0.08, damping=0.4, overlap=0)
	sources = dataInput['var1']
	targets = dataInput['var2']
	title = dataInput["var3"]
	consTerms = dataInput["var4"]
	weights = dataInput['value']
	edge_data = zip(sources, targets, title, consTerms, weights)
	for e in edge_data:
	    src = e[0]
	    dst = e[1]
	    ttl = e[2]
	    w = e[4]
	    if (e[1][0:2]=='rs'):
	        dstCol = "#996633"
	        dstshape ="triangle"
	        dstlevel = 1
	        dstmass = 2
	        titleEdge = e[3]
	    else:
	        dstCol = "#336666"
	        dstshape ="square"
	        dstlevel = 3
	        dstmass = 1	        
	    
	    got_net.add_node(src, src, title=src, 
	                     size= 20, 
	                     color='#660000', 
	                     level=2,
	                     mass = 1)
	    got_net.add_node(dst, dst, title=ttl, 
	                     size= 20, 
	                     color=dstCol,
	                     level = dstlevel,
	                     shape = dstshape,
	                     mass = dstmass)
	    if (e[1][0:2]=='rs'):
	    	got_net.add_edge(src, dst, value=w, color='#75D0E6', title = titleEdge)
	    else:
	    	got_net.add_edge(src, dst, value=w, color='#75D0E6')
	    neighbor_map = got_net.get_adj_list()

	got_net.save_graph("network.html")

	with open("network.html", "r") as f:
		contents = f.read()
		soup = BeautifulSoup(contents, 'html')
	scriptNetwork = soup.body.script
	scriptNetwork = '<div id = "mynetwork"></div>' + str(scriptNetwork)
	
	return(scriptNetwork)


