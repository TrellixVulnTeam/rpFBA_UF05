#!/bin/sh

docker run --network host -d -p 8888:8888 --name test_rpFBA brsynth/rpfba
sleep 10
python tool_rpFBA.py -inputTar test_input.tar -inSBML test_inSBML.sbml -sim_type fraction -reactions biomass RP1_sink -coefficients 0 0 -isMax True -fraction_of 0.95 -outputTar test_output.tar -dontMerge True -pathway_id rp_pathway -compartment_id MNXC3 -fill_orphan_species False -server_url http://0.0.0.0:8888/REST
python tool_rpFBA.py -sbml test_rpCofactors.sbml -inSBML test_inSBML.sbml -sim_type fraction -reactions biomass RP1_sink -coefficients 0 0 -isMax True -fraction_of 0.95 -outputTar test_output_single.tar -dontMerge True -pathway_id rp_pathway -compartment_id MNXC3 -fill_orphan_species False -server_url http://0.0.0.0:8888/REST
docker kill test_rpFBA
docker rm test_rpFBA
