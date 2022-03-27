from CGRtools.files import RDFRead, RDFWrite
from CGRtools.reactor import CGRReactor 
from CGRtools.containers import QueryContainer, MoleculeContainer, ReactionContainer
from collections import deque
from tqdm import tqdm
from sys import argv
from wrapt_timeout_decorator import *

def main(rdf_files, destination):
    
    getboroncount = lambda mol : sum(1 for atom in mol.atoms() if atom[1].atomic_number==5)

    boron = QueryContainer()
    boron.add_atom('C')
    boron.add_atom('B')
    boron.add_bond(1, 2, 1)
    deborylated = QueryContainer()
    deborylated.add_atom('C')
    template = ReactionContainer([boron], [deborylated])
    reactor = CGRReactor(template)
    
    @timeout(0.1)
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

    def is_borylation(rxn):
        borylation = False
        products = [prod for prod in rxn.products if getboroncount(prod)>0]
        reactants = rxn.reactants
        for product in products:
            deborylated = deborylation(product)
            for deborylated_prod in deborylated:
                borylation |= (deborylated_prod in reactants)
                
        return borylation
    def rdf_borylation_filter(filename):
        if '.rdf' not in filename:
            yield
        else:
            logname = filename.rsplit('.', maxsplit=1)[0]+'.log'
            log = open(logname, 'w')
        
        with RDFRead(filename) as f:
            for entry in tqdm(f, mininterval=0.5, desc=f'{filename}'):
                try:
                    entry.canonicalize()
                except:
                    continue
                try:
                    if is_borylation(entry):
                        yield entry
                except TypeError:
                    continue
                except TimeoutError:
                    log.write(str(entry)+'\n')
        log.close()
    
    with RDFWrite(destination) as f:
        for num, piece in enumerate(rdf_files):
            print('File {0} out of {1}'.format(num+1, len(rdf_files)))
            for rxn in rdf_borylation_filter(piece):
                f.write(rxn)

if __name__ == "__main__":
    main(argv[1:-1], argv[-1])