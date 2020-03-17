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
         full_sbml,
         output,
         pathway_id,
         compartment_id,
         sim_type,
         source_reaction,
         target_reaction,
         source_coefficient,
         target_coefficient,
         is_max,
         fraction_of,
         dont_merge):
    docker_client = docker.from_env()
    image_str = 'brsynth/rpfba-standalone'
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
        shutil.copy(full_sbml, tmpOutputFolder+'/input_sbml.dat')
        command = ['/home/tool_rpFBA.py',
                   '-input',
                   '/home/tmp_output/input.dat',
                   '-full_sbml',
                   '/home/tmp_output/input_sbml.dat',
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
                   '-fraction_of',
                   str(fraction_of),
                   '-dont_merge',
                   str(dont_merge),
                   '-pathway_id',
                   str(pathway_id),
                   '-output',
                   '/home/tmp_output/output.dat',
                   '-compartment_id',
                   str(compartment_id),
                   '-input_format',
                   str(input_format)]
        docker_client.containers.run(image_str, 
                command, 
                auto_remove=True, 
                detach=False, 
                volumes={tmpOutputFolder+'/': {'bind': '/home/tmp_output', 'mode': 'rw'}})
        shutil.copy(tmpOutputFolder+'/output.dat', output)


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to calculate FBA to generate rpFBA collection')
    parser.add_argument('-input', type=str)
    parser.add_argument('-input_format', type=str, default='tar')
    parser.add_argument('-full_sbml', type=str)
    parser.add_argument('-output', type=str)
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    parser.add_argument('-compartment_id', type=str, default='MNXC3')
    parser.add_argument('-sim_type', type=str, default='fraction')
    parser.add_argument('-source_reaction', type=str, default='biomass')
    parser.add_argument('-target_reaction', type=str, default='RP1_sink')
    parser.add_argument('-source_coefficient', type=float, default=1.0)
    parser.add_argument('-target_coefficient', type=float, default=1.0)
    parser.add_argument('-is_max', type=bool, default=True)
    parser.add_argument('-fraction_of', type=float, default=0.75)
    parser.add_argument('-dont_merge', type=bool, default=True)
    params = parser.parse_args()
    main(params.input,
         params.input_format,
         params.full_sbml,
         params.output,
         params.pathway_id,
         params.compartment_id,
         params.sim_type,
         params.source_reaction,
         params.target_reaction,
         params.source_coefficient,
         params.target_coefficient,
         params.is_max,
         params.fraction_of,
         params.dont_merge)
