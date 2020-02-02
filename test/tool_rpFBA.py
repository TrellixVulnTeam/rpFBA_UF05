#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpFBA REST service


python tool_rpFBA.py -input test_input.tar -full_sbml test_full_sbml.sbml -sim_type fraction -reactions biomass RP1_sink -coefficients 0 0 -is_max True -fraction_of 0.95 -output test_output.tar -dont_merge False -pathway_id rp_pathway -compartment_id MNXC3 -fill_orphan_species False

python tool_rpFBA.py -sbml test_rpCofactors.sbml -full_sbml test_full_sbml.sbml -sim_type fraction -reactions biomass RP1_sink -coefficients 0 0 -is_max True -fraction_of 0.95 -output test_output.tar -dont_merge False -pathway_id rp_pathway -compartment_id MNXC3 -fill_orphan_species False

"""
import requests
import argparse
import json
import logging
import tempfile
import tarfile
import glob
import shutil

##
#
#
def rpFBAUpload(input_tar, 
        full_sbml,
        sim_type,
        source_reaction,
        target_reaction,
        source_coefficient,
        target_coefficient,
        is_max,
        fraction_of,
        dont_merge,
        pathway_id,
        server_url, 
        output,
        compartment_id,
        fill_orphan_species):
    #convert dont_merge)
    if dont_merge=='False' or dont_merge==False:
        dont_merge = False
    elif dont_merge=='True' or dont_merge==True:
        dont_merge = True
    else:
        logging.error('Cannot recognise dont_merge')
    #convert fill_orphan_species
    if fill_orphan_species=='False' or fill_orphan_species==False:
        fill_orphan_species = False
    elif fill_orphan_species=='True' or fill_orphan_species==True:
        fill_orphan_species = True
    else:
        logging.error('Cannot recognise fill_orphan_species')
    # Post request
    data = {'sim_type': sim_type,
            'source_reaction': source_reaction,
            'target_reaction': target_reaction,
            'source_coefficient': source_coefficient,
            'target_coefficient': target_coefficient,
            'is_max': is_max,
            'fraction_of': fraction_of,
            'dont_merge': dont_merge, 
            'pathway_id': pathway_id,
            'compartment_id': compartment_id, 
            'fill_orphan_species': fill_orphan_species}
    files = {'input_tar': open(input_tar, 'rb'),
             'full_sbml': open(full_sbml, 'rb'),
             'data': ('data.json', json.dumps(data))}
    r = requests.post(server_url+'/Query', files=files)
    r.raise_for_status()
    with open(output, 'wb') as ot:
        ot.write(r.content)


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to calculate FBA to generate rpFBA collection')
    parser.add_argument('-input', type=str)
    parser.add_argument('-full_sbml', type=str)
    parser.add_argument('-sim_type', type=str)
    parser.add_argument('-source_reaction', type=str)
    parser.add_argument('-target_reaction', type=str)
    parser.add_argument('-source_coefficient', type=int)
    parser.add_argument('-target_coefficient', type=int)
    parser.add_argument('-is_max', type=str)
    parser.add_argument('-fraction_of', type=float)
    parser.add_argument('-dont_merge', type=str)
    parser.add_argument('-pathway_id', type=str)
    parser.add_argument('-server_url', type=str)
    parser.add_argument('-output', type=str)
    parser.add_argument('-compartment_id', type=str)
    parser.add_argument('-fill_orphan_species', type=str)
    parser.add_argument('-format', type=str)
    params = parser.parse_args()
    if params.format=='tar': 
        rpFBAUpload(params.input, 
                params.full_sbml,
                params.sim_type,
                params.source_reaction,
                params.target_reaction,
                params.source_coefficient,
                params.target_coefficient,
                params.is_max,
                params.fraction_of,
                params.dont_merge,
                params.pathway_id,
                params.server_url,
                params.output,
                params.compartment_id,
                params.fill_orphan_species)
    elif params.format=='sbml': 
        #make the tar.xz 
        with tempfile.TemporaryDirectory() as tmpOutputFolder:
            input_tar = tmpOutputFolder+'/tmp_input.tar.xz'
            output_tar = tmpOutputFolder+'/tmp_output.tar.xz'
            with tarfile.open(input_tar, mode='w:xz') as tf:
                tf.add(params.input)
            rpFBAUpload(input_tar, 
                    params.full_sbml,
                    params.sim_type,
                    params.source_reaction,
                    params.target_reaction,
                    params.source_coefficient,
                    params.target_coefficient,
                    params.is_max,
                    params.fraction_of,
                    params.dont_merge,
                    params.pathway_id,
                    params.server_url,
                    output_tar,
                    params.compartment_id,
                    params.fill_orphan_species)
            with tarfile.open(output_tar) as outTar:
                outTar.extractall(tmpOutputFolder)
            out_file = glob.glob(tmpOutputFolder+'/*.rpsbml.xml')
            if len(out_file)>1:
                logging.warning('There are more than one output file...')
            shutil.copy(out_file[0], params.output)
    else:
        logging.error('Cannot have no SBML and no TAR input')
        exit(0)
