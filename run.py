#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Extract the sink from an SBML into RP2 friendly format

"""
import argparse
import tempfile
import os
import logging
import shutil
import docker


##
#
#
def main(inputfile,
         input_format,
         gem_sbml,
         output,
         pathway_id='rp_pathway',
         species_group_id='central_species',
         sink_species_group_id='rp_sink_species',
		 objective_id='None',
         compartment_id='MNXC3',
         sim_type='fraction',
         source_reaction='biomass',
         target_reaction='RP1_sink',
         source_coefficient=1.0,
         target_coefficient=1.0,
		 num_workers=10,
         is_max=True,
         fraction_of=0.75,
         dont_merge=True):
    docker_client = docker.from_env()
    image_str = 'brsynth/rpfba-standalone:dev'
    try:
        image = docker_client.images.get(image_str)
    except docker.errors.ImageNotFound:
        logging.warning('Could not find the image, trying to pull it')
        try:
            docker_client.images.pull(image_str)
            image = docker_client.images.get(image_str)
        except docker.errors.ImageNotFound:
            logging.error('Cannot pull image: '+str(image_str))
            exit(1)
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        shutil.copy(inputfile, tmpOutputFolder+'/input.dat')
        shutil.copy(gem_sbml, tmpOutputFolder+'/gem_sbml.dat')
        command = ['/home/tool_rpFBA.py',
                   '-input',
                   '/home/tmp_output/input.dat',
                   '-gem_sbml',
                   '/home/tmp_output/gem_sbml.dat',
                   '-sim_type',
                   str(sim_type),
                   '-source_reaction',
                   str(source_reaction),
                   '-target_reaction',
                   str(target_reaction),
                   '-source_coefficient',
                   str(source_coefficient),
                   '-target_coefficient',
                   str(target_coefficient),
                   '-is_max',
                   str(is_max),
                   '-num_workers',
                   str(num_workers),
                   '-fraction_of',
                   str(fraction_of),
                   '-dont_merge',
                   str(dont_merge),
                   '-pathway_id',
                   str(pathway_id),
                   '-objective_id',
                   str(objective_id),
                   '-species_group_id',
                   str(species_group_id),
                   '-sink_species_group_id',
                   str(sink_species_group_id),
                   '-output',
                   '/home/tmp_output/output.dat',
                   '-compartment_id',
                   str(compartment_id),
                   '-input_format',
                   str(input_format)]
        container = docker_client.containers.run(image_str, 
                                                 command, 
                                                 detach=True, 
                                                 stderr=True,
                                                 volumes={tmpOutputFolder+'/': {'bind': '/home/tmp_output', 'mode': 'rw'}})
        container.wait()
        err = container.logs(stdout=False, stderr=True)
        err_str = err.decode('utf-8')
        if not 'ERROR' in err_str:
            if os.path.exists(tmpOutputFolder+'/output.dat'):
                shutil.copy(tmpOutputFolder+'/output.dat', output)
            else:
                logging.error('\n'+err_str)
                logging.error('rpFBA did not generate an output file')
        else:
            logging.error('\n'+err_str)
            logging.error('rpFBA did not generate an output file')
        container.remove()



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
    #TODO: check that the files exist
    main(params.input,
         params.input_format,
         params.gem_sbml,
         params.output,
         params.pathway_id,
         params.species_group_id,
         params.sink_species_group_id,
		 params.objective_id,
         params.compartment_id,
         params.sim_type,
         params.source_reaction,
         params.target_reaction,
         params.source_coefficient,
         params.target_coefficient,
		 params.num_workers,
         params.is_max,
         params.fraction_of,
         params.dont_merge)
