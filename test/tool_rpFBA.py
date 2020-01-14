#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpFBA REST service

"""
import requests
import argparse
import json


##
#
#
def rpFBAUpload(inputTar,
        pathway_id,
        dontMerge,
        inSBML,
        server_url,
        outputTar,
        compartment_id,
        fillOrphanSpecies):
    #convert dontMerge)
    if dontMerge=='False':
        dontMerge = False
    elif dontMerge=='True':
        dontMerge = True
    else:
        print('Cannot recognise dontMerge')
    #convert fillOrphanSpecies
    if fillOrphanSpecies=='False':
        fillOrphanSpecies = False
    elif fillOrphanSpecies=='True':
        fillOrphanSpecies = True
    else:
        print('Cannot recognise fillOrphanSpecies')
    # Post request
    data = {'pathway_id': pathway_id,
            'dontMerge': dontMerge,
            'compartment_id': compartment_id,
            'fillOrphanSpecies': fillOrphanSpecies}
    files = {'inputTar': open(inputTar, 'rb'),
             'inSBML': open(inSBML, 'rb'),
             'data': ('data.json', json.dumps(data))}
    r = requests.post(server_url+'/Query', files=files)
    r.raise_for_status()
    with open(outputTar, 'wb') as ot:
        ot.write(r.content)


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to calculate FBA to generate rpFBA collection')
    parser.add_argument('-inputTar', type=str)
    parser.add_argument('-inSBML', type=str)
    parser.add_argument('-dontMerge', type=str)
    parser.add_argument('-pathway_id', type=str)
    parser.add_argument('-server_url', type=str)
    parser.add_argument('-outputTar', type=str)
    parser.add_argument('-compartment_id', type=str)
    parser.add_argument('-fillOrphanSpecies', type=str)
    params = parser.parse_args()
    rpFBAUpload(params.inputTar,
            params.pathway_id,
            params.dontMerge,
            params.inSBML,
            params.server_url,
            params.outputTar,
            params.compartment_id,
            params.fillOrphanSpecies)
