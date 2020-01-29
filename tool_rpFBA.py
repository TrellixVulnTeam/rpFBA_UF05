#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpFBA REST service

"""
import argparse
import sys

sys.path.insert(0, '/home/')
import rpToolServe

##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to calculate FBA to generate rpFBA collection')
    parser.add_argument('-inputTar', type=str)
    parser.add_argument('-inSBML', type=str)
    parser.add_argument('-sim_type', type=str)
    parser.add_argument('-reactions', type=str, nargs='*')
    parser.add_argument('-coefficients', type=str, nargs='*')
    parser.add_argument('-isMax', type=str)
    parser.add_argument('-fraction_of', type=str)
    parser.add_argument('-outputTar', type=str)
    parser.add_argument('-dontMerge', type=str)
    parser.add_argument('-pathway_id', type=str)
    parser.add_argument('-compartment_id', type=str)
    parser.add_argument('-fill_orphan_species', type=str)
    params = parser.parse_args()
    rpToolServe.main(params.inputTar, 
            params.inSBML,
            params.outputTar,
            params.sim_type,
            params.reactions,
            params.coefficients,
            params.isMax,
            params.fraction_of,
            params.dontMerge,
            params.pathway_id,
            params.fill_orphan_species,
            params.compartment_id)
