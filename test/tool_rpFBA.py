#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpFBA REST service


python tool_rpFBA.py -inputTar test_input.tar -inSBML test_inSBML.sbml -sim_type fraction -reactions biomass RP1_sink -coefficients 0 0 -isMax True -fraction_of 0.95 -outputTar test_output.tar -dontMerge False -pathway_id rp_pathway -compartment_id MNXC3 -fill_orphan_species False

python tool_rpFBA.py -sbml test_rpCofactors.sbml -inSBML test_inSBML.sbml -sim_type fraction -reactions biomass RP1_sink -coefficients 0 0 -isMax True -fraction_of 0.95 -outputTar test_output.tar -dontMerge False -pathway_id rp_pathway -compartment_id MNXC3 -fill_orphan_species False

"""
import requests
import argparse
import json
import logging
import tempfile
import tarfile

##
#
#
def rpFBAUpload(inputTar, 
        inSBML,
        sim_type,
        reactions,
        coefficients,
        isMax,
        fraction_of,
        dontMerge,
        pathway_id,
        server_url, 
        outputTar,
        compartment_id,
        fill_orphan_species):
    #convert dontMerge)
    if dontMerge=='False' or dontMerge==False:
        dontMerge = False
    elif dontMerge=='True' or contMerge==True:
        dontMerge = True
    else:
        logging.error('Cannot recognise dontMerge')
    #convert fill_orphan_species
    if fill_orphan_species=='False' or fill_orphan_species==False:
        fill_orphan_species = False
    elif fill_orphan_species=='True' or fill_orphan_species==True:
        fill_orphan_species = True
    else:
        logging.error('Cannot recognise fill_orphan_species')
    # Post request
    data = {'sim_type': sim_type,
            'reactions': reactions,
            'coefficients': coefficients,
            'isMax': isMax,
            'fraction_of': fraction_of,
            'dontMerge': dontMerge, 
            'pathway_id': pathway_id,
            'compartment_id': compartment_id, 
            'fill_orphan_species': fill_orphan_species}
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
    parser.add_argument('-sbml', type=str)
    parser.add_argument('-inSBML', type=str)
    parser.add_argument('-sim_type', type=str)
    parser.add_argument('-reactions', type=str, nargs='*')
    parser.add_argument('-coefficients', type=str, nargs='*')
    parser.add_argument('-isMax', type=str)
    parser.add_argument('-fraction_of', type=str)
    parser.add_argument('-dontMerge', type=str)
    parser.add_argument('-pathway_id', type=str)
    parser.add_argument('-server_url', type=str)
    parser.add_argument('-outputTar', type=str)
    parser.add_argument('-compartment_id', type=str)
    parser.add_argument('-fill_orphan_species', type=str)
    params = parser.parse_args()
    if params.sbml=='None' or params.sbml==None or params.sbml=='':
        if params.inputTar=='None' or params.inputTar==None or params.inputTar=='':
            logging.error('Cannot have no SBML and no TAR input')
            exit(0)
        rpFBAUpload(params.inputTar, 
                params.inSBML,
                params.sim_type,
                params.reactions,
                params.coefficients,
                params.isMax,
                params.fraction_of,
                params.dontMerge,
                params.pathway_id,
                params.server_url,
                params.outputTar,
                params.compartment_id,
                params.fill_orphan_species)
    else:
        #make the tar.xz 
        with tempfile.TemporaryDirectory() as tmpOutputFolder:
            inputTar = tmpOutputFolder+'/tmp_input.tar.xz'
            with tarfile.open(inputTar, mode='w:xz') as tf:
                tf.add(params.sbml)
            rpFBAUpload(inputTar, 
                    params.inSBML,
                    params.sim_type,
                    params.reactions,
                    params.coefficients,
                    params.isMax,
                    params.fraction_of,
                    params.dontMerge,
                    params.pathway_id,
                    params.server_url,
                    params.outputTar,
                    params.compartment_id,
                    params.fill_orphan_species)
