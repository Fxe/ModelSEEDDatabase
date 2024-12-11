#!/usr/bin/env python
import re, sys
#sys.path.append('../../../Libs/Python/')
#from BiochemPy import Compounds

###########################################################
# Parse pathways
###########################################################

reaction_file = "data/pathway"
pathway_dict = dict()
pathways_list = list()
with open(reaction_file) as cfh:
	for line in cfh.readlines():
		line=line.strip('\r\n')
		tmp_list=re.split('\s+',line)
		data=" ".join(tmp_list[1:])
		if(tmp_list[0] != ""):
			field = tmp_list[0]

		if(field == "ENTRY"):
			pathway_dict = {'id':tmp_list[1],
					 'name':"",'description':"",
					 'reactions':[]}

		elif(field == "NAME"):
			data=data.rstrip(';')
			pathway_dict['name']=data
		elif(field == "DESCRIPTION"):
			pathway_dict['description']=data
		elif(field == "REACTION"):
			rxn=data.split(" ")[0]
			pathway_dict['reactions'].append(rxn)
		else:
			#print(field)
			pass

		if(field == "///"):

			if(len(pathway_dict['reactions']) == 0):
				continue

			pathways_list.append(pathway_dict)

with open("KEGG_pathways.tsv",'w') as kpfh:
	kpfh.write("\t".join(["id","name","description","reactions"])+"\n")
	for pwy in pathways_list:
		kpfh.write("\t".join([pwy['id'], pwy['name'], pwy['description'], "|".join(pwy['reactions'])])+"\n")
