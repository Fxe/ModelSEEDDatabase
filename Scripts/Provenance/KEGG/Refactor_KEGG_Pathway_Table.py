#!/usr/bin/env python
import re, sys
#sys.path.append('../../../Libs/Python/')
#from BiochemPy import Compounds

###########################################################
# Parse pathways
###########################################################

pathway_file = "data/pathway"
pathways_dict = dict()
pathway_id=None
with open(pathway_file) as cfh:
	for line in cfh.readlines():
		line=line.strip('\r\n')
		tmp_lst=re.split('\s+',line)
		data=" ".join(tmp_lst[1:])
		if(tmp_lst[0] != ""):
			field = tmp_lst[0]

		if(field == "ENTRY"):
			pathway = tmp_lst[1]
			pathways_dict[pathway] = {'name':"",'description':"",
						  'reactions':[],
						  'parent':[]}

		elif(field == "NAME"):
			data=data.rstrip(';')
			pathways_dict[pathway]['name']=data
		elif(field == "DESCRIPTION"):
			pathways_dict[pathway]['description']=data
		elif(field == "REACTION"):
			rxn=data.split(" ")[0]
			pathways_dict[pathway]['reactions'].append(rxn)
		else:
			#print(field)
			pass

		if(field == "///"):

			if(len(pathways_dict[pathway]['reactions']) == 0):
				del(pathways_dict[pathway])

module=None
module_file = "data/module"
with open(module_file) as cfh:
	for line in cfh.readlines():
		line=line.strip('\r\n')
		tmp_lst=re.split('\s+',line)
		data=" ".join(tmp_lst[1:])
		if(tmp_lst[0] != ""):
			field = tmp_lst[0]

		if(field == "ENTRY"):
			module = tmp_lst[1]
			pathways_dict[module] = {'name':"",'description':"",
						 'reactions':[],
						 'parent':[]}

		elif(field == "NAME"):
			data=data.rstrip(';')
			pathways_dict[module]['name']=data
		elif(field == "DESCRIPTION"):
			pathways_dict[module]['description']=data
		elif(field == "REACTION"):
			rxn=data.split(" ")[0]
			pathways_dict[module]['reactions'].append(rxn)
		elif(field == "PATHWAY"):
			pwy=data.split(" ")[0]
			pwy = pwy.replace('map','rn')
			if(pwy not in pathways_dict):
				print(module,pwy)
			else:
				pathways_dict[module]['parent'].append(pwy)
		else:
			#print(field)
			pass

		if(field == "///"):

			if(len(pathways_dict[module]['reactions']) == 0):
				del(pathways_dict[module])

with open("KEGG_pathways.tsv",'w') as kpfh:
	kpfh.write("\t".join(["id","name","description","reactions","parent"])+"\n")
	for pwy in pathways_dict:
		kpfh.write("\t".join([pwy, pathways_dict[pwy]['name'], pathways_dict[pwy]['description'], "|".join(pathways_dict[pwy]['reactions']), "|".join(pathways_dict[pwy]['parent'])])+"\n")
