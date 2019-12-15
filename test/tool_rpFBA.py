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
    parser.add_argument('-outputTar', type=str)
    parser.add_argument('-dontMerge', type=str)
    parser.add_argument('-pathway_id', type=str)
    parser.add_argument('-fillOrphanSpecies', type=str)
    parser.add_argument('-compartment_id', type=str)
    params = parser.parse_args()
    rpToolServe.main(params.inputTar,
                     params.inSBML,
                     params.outputTar,
                     params.dontMerge,
                     params.pathway_id,
                     params.fillOrphanSpecies,
                     params.compartment_id)
