from CGRtools.files import RDFRead, SMILESRead
from CGRtools.reactor import CGRReactor 
from CGRtools.containers import QueryContainer, MoleculeContainer, ReactionContainer
from io import StringIO
from collections import deque

getboroncount = lambda mol : sum(1 for atom in mol.atoms() if atom[1].atomic_number==5)

boron = QueryContainer()
boron.add_atom('C')
boron.add_atom('B')
boron.add_bond(1, 2, 1)
deborylated = QueryContainer()
deborylated.add_atom('C')
template = ReactionContainer([boron], [deborylated])
reactor = CGRReactor(template)

def deborylation(mol):
    result = set([mol])
    d = deque([mol])
    while len(d)>0:
        current_deborylated = set(r for r in reactor(d.popleft()))
        for prod in current_deborylated:
            if prod not in result:
                result.add(prod)
                d.append(prod)
    result.discard(mol)    
    return list(result)

#takes in ReactionContainers from a successfully parsed RDF file
def substrate_product(rxn, boron=False):
    products = [prod for prod in rxn.products if getboroncount(prod)>0]
    reactants = rxn.reactants
    for product in products:
        deborylated = deborylation(product)
        for deborylated_prod in deborylated:
            if deborylated_prod in reactants:
                if boron==False:
                    return deborylated_prod, product
                boron_source_list = [rct for rct in rxn.reactants if (getboroncount(rct)>0 and rct!=sub)]
                if len(boron_source_list)<1:
                    return None
                for mol in boron_source_list:
                    if str(mol) in ['CC1(C)OBOC1(C)C', 'CC1(C)OB(OC1(C)C)B2OC(C)(C)C(C)(C)O2']:
                        return deborylated_prod, mol, product
                    return deborylated_prod, boron_source_list[0], product
    return None

def is_borylation(rxn):
    return substrate_product(rxn) is not None