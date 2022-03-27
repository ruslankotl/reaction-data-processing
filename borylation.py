# determine whether the reaction is a single-step C–H borylation

from rdkit import Chem
from rdkit.Chem import rdChemReactions
from collections import deque
from itertools import chain


#determines if molecules are identical
def molecules_identical(mol1, mol2):
    return mol1.HasSubstructMatch(mol2) & mol2.HasSubstructMatch(mol1)


deborylation = rdChemReactions.ReactionFromSmarts('[#6:1][B]>>[#6:1]')

getboroncount = lambda mol : sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum()==5)

def unique_deborylation_single_step(mol):
    deborylation_products = list(chain.from_iterable(deborylation.RunReactants([mol])))
    deborylation_unique_products = []
    
    for prod in deborylation_products:
        unique = True
        for u_prod in deborylation_unique_products:
            unique ^= molecules_identical(prod, u_prod)
        if unique:
            deborylation_unique_products.append(prod)
    
    return deborylation_unique_products
    

def deque_deborylation(mol):
    
    if getboroncount(mol)==1:
        return unique_deborylation_single_step(mol)
    
    result = [mol]
    d = deque(result)
    
    while (len(d)>0):
        current_deborylated = unique_deborylation_single_step(d.popleft())
        for prod in current_deborylated:
            
            unique = True
            for u_prod in result:
                unique ^= molecules_identical(prod, u_prod)
            if unique:
                result.append(prod)
                d.append(prod)
        
    
    return result[1:]

def substrate_product(rxn, sanitize=True):
    
    borylation = False
    products = [prod for prod in rxn.GetProducts() if getboroncount(prod)>0]
    reactants = rxn.GetReactants()
    if sanitize:
        [Chem.SanitizeMol(prod) for prod in products]
        [Chem.SanitizeMol(reac) for reac in reactants]
    for product in products:
        deborylated_products = deque_deborylation(product)
        for deborylated_product in deborylated_products:
            for reactant in reactants:
                borylation = molecules_identical(reactant, deborylated_product)
                if borylation:
                    return reactant, product
    return None
        
def is_borylation(rxn):
    return substrate_product(rxn) is not None
