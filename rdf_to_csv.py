from sys import argv
import csv
from CGRtools import RDFRead
from tqdm import tqdm

def main(filename, destination):
  
    with RDFRead(filename) as rdf:
        data = [rxn for rxn in tqdm(rdf, desc=f'Parsing from {filename}')]
        [rxn.canonicalize() for rxn in tqdm(data, desc='Canonicalising')]
    
    possible_keys=set()
    for rxn in data:
        for key in rxn.meta.keys():
            if 'RX' in key:
                possible_keys.add(key.rsplit(':', maxsplit=1)[1])
    
    fields = ['RX_ID', 'RX_SMI', 'RX_NVAR', 'RX_RANK', 'RX_RTYP', 'RX_RCT', 'RX_PRO', 'RX_BIN', 'CL', 'LB', 'TI', 'TXT', 'STP', 'SNR', 'MID', 'MTEXT', 'TIM', 'T', 'P', 'PH', 'COND', 'TYP', 'SUB', 'PRT', 'RXDES', 'LCN', 'COM', 'YPRO', 'YD', 'NYD', 'YDO', 'SRCT', 'RGT', 'CAT', 'SOL', 'citation']
    catch_fields = list(possible_keys - set(fields))
    if len(catch_fields)>0:
        for elem in catch_fields:
            fields.append(elem)
    
    with open(destination, 'w+') as csvfile:
        writer = csv.DictWriter(csvfile, restval='',fieldnames=fields, quoting=csv.QUOTE_NONNUMERIC)
    
        writer.writeheader()
    
        for rxn in tqdm(data, desc = f'Writing to {destination}'):
            try:
                rxn_dict = {'RX_SMI': str(rxn)}
                for key in rxn.meta:
                    if 'ROOT:RX_' in key:
                        rxn_dict.update({key.rsplit(':', maxsplit=1)[1]: rxn.meta[key]})
                try:
                    num_variants = int(rxn.meta['ROOT:RX_NVAR'])
                except KeyError:
                    num_variants = 1
                for variant in range(num_variants):
                    var_dict = dict()
                    var_string = 'ROOT:RXD({0})'.format(variant+1)
                    for var_key in rxn.meta.keys():
                        if var_string in var_key:
                            var_dict.update({var_key.rsplit(':', maxsplit=1)[1]: rxn.meta[var_key]})
                    var_dict.update(rxn_dict.copy())
                    writer.writerow(var_dict)
            except KeyError:
                continue

if __name__ == "__main__":
    main(argv[1], argv[2])