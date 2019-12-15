
docker run -v ${PWD}/test_run.sh:/home/test_run.sh -v ${PWD}/tool_rpFBA.py:/home/tool_rpFBA.py -v ${PWD}/test_input.tar:/home/test_input.tar -v ${PWD}/test_inSBML.sbml:/home/test_inSBML.sbml -v ${PWD}/results/:$/home/results/ --rm --user root brsynth/rpfba /bin/sh test_run.sh


docker run -v ${PWD}/test_run.sh:/home/test_run.sh -v ${PWD}/tool_rpFBA.py:/home/tool_rpFBA.py -v ${PWD}/test_input.tar:/home/test_input.tar -v ${PWD}/test_inSBML.sbml:/home/test_inSBML.sbml -v ${PWD}/results/:$/home/results/--rm --user root brsynth/rpfba /bin/sh test_run.sh


docker run "GALAXY_SLOTS=$GALAXY_SLOTS" -v ${PWD}/test_input.tar:/home/mdulac/test_input.tar:ro --rm --user root brsynth/rpoptbiodes /bin/sh test_run.sh 
