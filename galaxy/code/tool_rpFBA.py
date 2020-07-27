#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpFBA REST service

"""
import argparse
import sys
import shutil
import logging
import tempfile
import tarfile
import glob
import os

sys.path.insert(0, '/home/')
import rpToolServe

##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to calculate FBA to generate rpFBA collection')
    parser.add_argument('-input', type=str)
    parser.add_argument('-input_format', type=str)
    parser.add_argument('-gem_sbml', type=str)
    parser.add_argument('-output', type=str)
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    parser.add_argument('-sink_species_group_id', type=str, default='rp_sink_species')
    parser.add_argument('-species_group_id', type=str, default='central_species')
    parser.add_argument('-objective_id', type=str, default='None')
    parser.add_argument('-compartment_id', type=str, default='MNXC3')
    parser.add_argument('-sim_type', type=str, default='fraction')
    parser.add_argument('-source_reaction', type=str, default='biomass')
    parser.add_argument('-target_reaction', type=str, default='RP1_sink')
    parser.add_argument('-source_coefficient', type=float, default=1.0)
    parser.add_argument('-target_coefficient', type=float, default=1.0)
    parser.add_argument('-num_workers', type=int, default=10)
    parser.add_argument('-is_max', type=str, default='True')
    parser.add_argument('-fraction_of', type=float, default=0.75)
    parser.add_argument('-dont_merge', type=str, default='True')
    params = parser.parse_args()
    if params.num_workers>=20 or params.num_workers<=0:
        logging.error('Cannot have more than 20 and 0 or less workers: '+str(params.num_workers))
        exit(1)
    if params.fraction_of<=0.0:
        logging.error('Cannot have -fraction_of less or equal than 0: '+str(params.fraction_of))
        exit(1)
    if params.is_max==True or params.is_max=='True' or params.is_max=='true':
        is_max = True
    elif params.is_max==False or params.is_max=='False' or params.is_max=='false':
        is_max = False
    else:
        logging.error('Cannot interpret '+str(params.is_max))
        exit(1)
    if params.dont_merge==True or params.dont_merge=='True' or params.dont_merge=='true':
        dont_merge = True
    elif params.dont_merge==False or params.dont_merge=='False' or params.dont_merge=='false':
        dont_merge = False
    else:
        logging.error('Cannot interpret '+str(params.dont_merge))
        exit(1)
    if params.objective_id=='None':
        objective_id = 'obj_'+params.sim_type
    else:
        objective_id = params.objective_id
    if params.input_format=='tar':
        rpToolServe.main(params.input,
                         params.gem_sbml,
                         params.output,
                         params.sim_type,
                         params.source_reaction,
                         params.target_reaction,
                         params.source_coefficient,
                         params.target_coefficient,
                         is_max,
                         params.fraction_of,
                         dont_merge,
                         params.num_workers,
                         params.pathway_id,
                         objective_id,
                         params.compartment_id,
                         None, #this is fillorphanspecies
                         params.species_group_id,
                         params.sink_species_group_id)
    elif params.input_format=='sbml':
        #make the tar.xz 
        with tempfile.TemporaryDirectory() as tmpOutputFolder:
            input_tar = tmpOutputFolder+'/tmp_input.tar'
            output_tar = tmpOutputFolder+'/tmp_output.tar'
            with tarfile.open(input_tar, mode='w:gz') as tf:
                #tf.add(params.input)
                info = tarfile.TarInfo('single.rpsbml.xml') #need to change the name since galaxy creates .dat files
                info.size = os.path.getsize(params.input)
                tf.addfile(tarinfo=info, fileobj=open(params.input, 'rb'))
            rpToolServe.main(input_tar,
                             params.gem_sbml,
                             output_tar,
                             params.sim_type,
                             params.source_reaction,
                             params.target_reaction,
                             params.source_coefficient,
                             params.target_coefficient,
                             is_max,
                             params.fraction_of,
                             dont_merge,
                             params.num_workers,
                             params.pathway_id,
                             objective_id,
                             params.compartment_id,
                             None, #this is filloorphanspecies
                             params.species_group_id,
                             params.sink_species_group_id)
            with tarfile.open(output_tar, mode='r') as outTar:
                outTar.extractall(tmpOutputFolder)
            out_file = glob.glob(tmpOutputFolder+'/*.sbml.xml')
            if len(out_file)>1:
                logging.warning('There are more than one output file...')
            shutil.copy(out_file[0], params.output)
    else:
        logging.error('Cannot have no SBML and no TAR input')
        exit(0)
