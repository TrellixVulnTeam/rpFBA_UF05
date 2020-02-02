#!/bin/sh

docker run --network host -d -p 8888:8888 --name test_rpFBA brsynth/rpfba
sleep 10
python tool_rpFBA.py -inputTar test_input.tar -inSBML test_inSBML.sbml -sim_type fraction -reactions biomass RP1_sink -coefficients 0 0 -isMax True -fraction_of 0.95 -outputTar test_output.tar -dontMerge True -pathway_id rp_pathway -compartment_id MNXC3 -fill_orphan_species False -server_url http://0.0.0.0:8888/REST

python tool_rpFBA.py -input rp_25_1.rpsbml.xml -full_sbml test_inSBML.sbml -sim_type fraction -source_reaction biomass -target_reaction RP1_sink -is_max True -fraction_of 0.95 -output test_output_single.sbml -dont_merge True -pathway_id rp_pathway -compartment_id MNXC3 -fill_orphan_species False -server_url http://0.0.0.0:8888/REST -format sbml

python tool_rpFBA.py -input test_input.tar -full_sbml test_inSBML.sbml -sim_type fraction -source_reaction biomass -target_reaction RP1_sink -source_coefficient 1 -target_coefficient 1 -is_max True -fraction_of 0.95 -output test_output.tar -dont_merge True -pathway_id rp_pathway -compartment_id MNXC3 -fill_orphan_species False -server_url http://0.0.0.0:8884/REST -format tar

docker kill test_rpFBA
docker rm test_rpFBA
