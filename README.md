# rpFBA

* Docker Image: [brsynth/rpfba-standalone](https://hub.docker.com/r/brsynth/rpfba-standalone)

Perform FBA on a single or collection of SBML files containing heterologous pathways, as tar.xz archives. The package performs the following steps: 1) it merges a user defined GEM SBML model to a given heterologous pathway. 2) it performs FBA using the [cobrapy](https://opencobra.github.io/cobrapy/) package using a user defined mathod that include, FBA, parsimonious FBA or fraction of optimum of another reaction. For the first two, the user must know the reaction name that the model will optimise to, while the latter the use must provide the target reaction but also another reaction that will be restricted. The first step involves performing FBA using the "source" reaction as the objective. Then the flux of that reaction has its upper and lower bounds set to the same value, determined as a fraction of its FBA flux value. Thereafter the objective is set to the initial target reaction and FBA is performed once again. The tool uses the [FBC](https://co.mbine.org/specifications/sbml.level-3.version-1.fbc.version-2.release-1) package to manage the objective and flux bounds.

## Input

Required:
* **-input**: (string) Path to the input file
* **-input_format**: (string) Valid options: tar, sbml. Format of the input file
* **-full_sbml**: (string) Path to the GEM SBML model

Advanced options:
* **-pathway_id**: (string, default=rp_pathway) ID of the heterologous pathway
* **-compartment_id**: (string, default=MNXC3 (i.e. cytoplasm)) ID of the compartment ID that contains the heterologous pathway
* **-sim_type**: (string, default=fraction) Valid options include: fraction, fba, pfba. The type of constraint based modelling method 
* **-source_reaction**: (string, default=biomass) Name of the source reaction that will be restricted in the "fraction" simulation type. This parameter is ignored for "fba" and "pfba"
* **-target_reaction**: (string, default=RP1_sink) Heterologous pathway flux sink reaction. This parameters is required in all simulation type
* **-source_coefficient**: (float, default=1.0) Objective coefficient for the source reaction. This parameter is ignored for "fba" and "pfba"
* **-target_coefficient**: (float, default=1.0) Objective coefficient for the target reaction. 
* **-is_max**: (boolean, default=True) Maximise or minimise the objective function
* **-fraction_of**: (float, default=0.75) Portion of the maximal flux used to set the maximal and minimal bounds for the source reaction of the "fraction" simulation type
* **-dont_merge**: (boolean, default=True) Return the merged GEM+heterologous pathway SBML or only the heterologous pathway SBML files

## Output

* **-output**: (string) Path to the output file

## Building the docker

```
docker build -t brsynth/rpfba-standalone -f Dockerfile .
```

## Running the tests

To test the model extract the test.tar and run the following command:

```
python run.py -input test/test_rpCofactors.tar -input_format tar -full_sbml test/e_coli_model.sbml -output test/test_rpFBA.tar
```

## Dependencies

* Base docker image: [brsynth/rpBase](https://hub.docker.com/r/brsynth/rpbase)

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

v0.1

## Authors

* **Melchior du Lac** 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thomas Duigou
* Joan Hérisson
