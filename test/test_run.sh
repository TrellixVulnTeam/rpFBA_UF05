#!/bin/sh

python tool_rpFBA.py -inputTar test_input.tar -outputTar test_output.tar -inSBML test_inSBML.sbml -dontMerge True -pathway_id pathway_id  -compartment_id MNXC3 -fillOrphanSpecies False
mv test_output.tar results/
