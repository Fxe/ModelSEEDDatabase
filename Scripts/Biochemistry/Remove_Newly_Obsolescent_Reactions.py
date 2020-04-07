#!/usr/bin/env python
import os, sys
temp=list();
header=1;

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions, Compounds, InChIs

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()
rxns_aliases_dict = reactions_helper.loadMSAliases()
rxns_names_dict = reactions_helper.loadNames()
rxns_ecs_dict = reactions_helper.loadECs()

# We actually don't want obsolete reactions and compounds in our database
# So we're striving to remove any 'new' ones that are obsolete
# Any information attached to them should be associated with their linked counterpart

last_rxn_str='rxn40000'
last_rxn_int=int(last_rxn_str[3:])

delete_rxns=list()
for rxn in reactions_dict:
    rxn_int = int(rxn[3:])
    if(rxn_int > last_rxn_int and reactions_dict[rxn]['is_obsolete']):
        delete_rxns.append(rxn)

for rxn in delete_rxns:

    # I need to update the linked_reaction field to make sure that the removed reaction is deleted
    for lnkd_rxn in reactions_dict[rxn]['linked_reaction']:
        reactions_dict[lnkd_rxn]['linked_reaction'].remove(rxn)

        # I need to move the names, aliases, and ec numbers to linked reactions
        # This should already have been done when merging reactions, but doing it here to double-check
        if(rxn in rxns_aliases_dict):
            for source in rxns_aliases_dict[rxn]:
                if(source not in rxns_aliases_dict[lnkd_rxn]):
                    rxns_aliases_dict[lnkd_rxn][source]=list()
                for alias in rxns_aliases_dict[rxn][source]:
                    if(alias not in rxns_aliases_dict[lnkd_rxn][source]):
                        rxns_aliases_dict[lnkd_rxn][source].append(alias)
                        print("Warning: adding "+alias+" for "+source+" to "+lnkd_rxn)

        if(rxn in rxns_names_dict):
            for name in rxns_names_dict[rxn]:
                if(name not in rxns_names_dict[lnkd_rxn]):
                    rxns_names_dict[lnkd_rxn].append(name)
                    print("Warning: adding "+name+" to "+lnkd_rxn)

        if(rxn in rxns_ecs_dict):
            for ec in rxns_ecs_dict[rxn]:
                if(ec not in rxns_ecs_dict[lnkd_rxn]):
                    rxns_ecs_dict[lnkd_rxn].append(ec)
                    print("Warning: adding "+ec+" to "+lnkd_rxn)

    del(reactions_dict[rxn])
    del(rxns_aliases_dict[rxn])
    del(rxns_ecs_dict[rxn])
    del(rxns_names_dict[rxn])

if(len(delete_rxns)>0):
    print("Removing "+str(len(delete_rxns))+" newly obsolete reactions")
    reactions_helper.saveReactions(reactions_dict)