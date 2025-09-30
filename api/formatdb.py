
# coding: utf-8

# In[1]:

import pandas as pd
#from classyfire import *
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rdMD
from urllib.request import urlopen
from urllib.request import HTTPError
import json
import re
import os
#from api.npclassify import classify_structure


def query_inchikey(keylist):
    count_failed = 0
    dlist = []
    dmetadatalist = []
    for i in range(len(keylist)):
        key = keylist[i]
        url = 'http://classyfire.wishartlab.com/entities/' + key + '.json'

        try:
            r = urlopen(url)
            fn = r.read().decode('utf8')
        except  HTTPError as err:
            fn = ''
        #try:
        if fn != '':
            mdict = {}
            data = json.loads(fn)

            mdict['inchikey'] = key
            if 'kingdom' not in data.keys():
                mdict['kingdom'] = ''
                mdict['superclass'] = ''
                mdict['class'] = ''
                mdict['subclass'] = ''
                mdict['direct_parent'] = ''
                mdict['molecular_framework'] = ''
                dmetadatalist.append(mdict)
                continue

            if data['kingdom']!=None:
                mdict['kingdom'] = data['kingdom']['name']
            else:
                mdict['kingdom'] = ''
            if data['superclass']!=None:
                mdict['superclass'] = data['superclass']['name']
            else:
                mdict['superclass'] = ''
            if data['class']!=None:
                mdict['class'] = data['class']['name']
            else:
                mdict['class'] = ''
            if data['subclass']!=None:
                mdict['subclass'] = data['subclass']['name']
            else:
                mdict['subclass'] = ''
            if data['direct_parent'] is not None:
                mdict['direct_parent'] = data['direct_parent']['name']
            else:
                mdict['direct_parent'] = ''
            if data['molecular_framework'] is not None:
                mdict['molecular_framework'] = data['molecular_framework']
            else:
                mdict['molecular_framework'] = ''
            dmetadatalist.append(mdict)
        #except HTTPError as err:
        else:
            count_failed += 1
            dmetadatalist.append({})

    df_metares = pd.DataFrame.from_dict(dmetadatalist)
    return(df_metares)


def formatdb(smiles, tool, oformat):
    df = pd.read_csv(smiles,  sep='\t', header=None)
    os.remove(smiles)

    smi = list(df[0])
    m = [Chem.MolFromSmiles(x) for x in smi]
    inchi = []
    ikeys = []
    ikey1 = []
    ikey2 = []
    form = []
    exmass = []
    for i in range(len(m)):
        try:
            inchi.append(Chem.rdinchi.MolToInchi(m[i])[0])
            ikey = Chem.rdinchi.InchiToInchiKey(inchi[i])
            ikeys.append(ikey)
            if ikey!='':
                ikey1.append(ikey.split('-')[0])
                ikey2.append(ikey.split('-')[1])
                form.append(rdMD.CalcMolFormula(m[i]))
                exmass.append(rdMD.CalcExactMolWt(m[i]))
            else:
                ikey1.append('')
                ikey2.append('')
                form.append('')
                exmass.append('')
        except:
            ikeys.append('')
            inchi.append('')
            ikey1.append('')
            ikey2.append('')
            form.append('')
            exmass.append('')

    data = {'inchikey': ikeys, 'MonoisotopicMass': exmass, 'InChI': inchi, 'SMILES': list(df[0]),
                  'Identifier': list(df[1]), 'InChIKey2': ikey2, 'InChIKey1': ikey1, 'MolecularFormula': form}

    cn = ["inchikey", "MonoisotopicMass", "InChI", "SMILES", "Identifier", "InChIKey2", "InChIKey1", "MolecularFormula"]
    formdata = pd.DataFrame(data, columns=cn)
    if tool=='ClassyFire':
        classy = query_inchikey(ikeys)
    elif tool=='NPClassifier':
        classy =[]
        for i in formdata.index:
            isglycoside, class_results, superclass_results, pathway_results, path_from_class, path_from_superclass, n_path, fp1, fp2 = classify_structure(formdata.loc[i, 'SMILES'])
            if len(superclass_results):
                 supclass = superclass_results[0]
            else:
                 supclass = ''

            if len(class_results):
                 classe = class_results[0]
            else:
                 classe = ''

            classy.append([formdata.loc[i, 'inchikey'], '', supclass, classe, ''])
        classy = pd.DataFrame(classy)
        classy.columns = ['inchikey', 'kingdom', 'superclass', 'class', 'subclass']
    else:
        raise Exception("Invalid tool option.")

    # If the structure do not show a classification, try query
    #in_process = get_class(list(df[0]), chunksize=100)
    #classy = poll(in_process)

    classy = classy[['inchikey', 'kingdom', 'superclass', 'class', 'subclass']]
    classy.columns = ['inchikey', 'kingdom_name', 'superclass_name', 'class_name', 'subclass_name']

    formfinal = pd.merge(formdata, classy, how='left', on=['inchikey'])

    formfinal = formfinal.fillna('')
    formfinal.drop('inchikey', axis=1, inplace=True)

    id = [x for x in range(len(ikeys)) if ikeys[x]=='']
    formfinal.drop(formfinal.index[id], inplace=True)

    if oformat=='TSV':
        formfinal.to_csv(smiles+'_FORMATED.txt', index=False, sep='\t')
    elif oformat=='PSV':
        formfinal.to_csv(smiles+'_FORMATED.txt', index=False, sep='|')
    else:
        raise Exception("Format should be either TSV or PSV.")

    return 'Done'
