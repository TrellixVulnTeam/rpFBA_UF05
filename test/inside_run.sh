#!/bin/bash

python tool_rpFBA.py -inputTar test_input.tar -inSBML test_inSBML.sbml -sim_type fraction -reactions biomass RP1_sink -coefficients 0 0 -isMax True -fraction_of 0.95 -outputTar test_output.tar -dontMerge False -pathway_id rp_pathway -compartment_id MNXC3 -fill_orphan_species False
mv test_output.tar results/
