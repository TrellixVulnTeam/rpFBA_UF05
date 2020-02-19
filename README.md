# rpFBA

Reads a single or collection of rpSBML files (as a tar.xz) and a GEM SBML model, merges the two and perform FBA. The addition of the heterologous pathway to the organism of choice. FBA is performed with all defined objectives defined by the [FBC](https://co.mbine.org/specifications/sbml.level-3.version-1.fbc.version-2.release-1) SBML package and saved inside custom (BRSynth) annotations. The output is also a collection of rpSBML files in a tar.xz. The merged version of the model may be kept instead of the heterologous pathway alone (default). 

Currently, the following FBA methods have been implemented. FBA, parsimonious FBA and fraction of reaction. The latter forces a source reaction to a fraction of optimum, and then perform FBA using another target reaction.

## Information Flow

### Input

Required information:
* Fraction of Reaction:
    * **Fraction of the source reaction**: (default: 0.75) Float between 0.0 and 1.0, with 0.0 being no flux to 1.0 maximum of optimum
    * **Source Reaction**: (default: biomass) ID of the source reaction
    * **Source Coefficient**: (default: 1.0) Float determining the FBA coefficient
    * **Target Reaction**: (default: RP1_sink) ID of the target reaction
    * **Target Coefficient**: (default: 1.0) Float determining the FBA coefficient
* Parsimonious FBA:
    * **Fraction of optimal**: (default: 0.95) 
    * **Target Reaction**: (default: RP1_sink) ID of the target reaction
    * **Target Coefficient**: (default: 1.0) Float determining the FBA coefficient
* FBA:
    * **Target Reaction**: (default: RP1_sink) ID of the target reaction
    * **Target Coefficient**: (default: 1.0) Float determining the FBA coefficient
* **Input/output format**: (default: TAR) Either use TAR collection of files or single rpSBML file
* **Input GEM SBML model**: SBML model

Advanced options:
* **Name of the heterologous pathway**: (default: rp_pathway) Groups ID of the heterologous pathway
* **Maximise the objective?**: (default: Yes) If set to No, then minimise the objective
* **Don't output the merged model?**: (default: Yes) If Not then the full merged model is output
* **IP address of the rpFBA REST service**: IP address of the REST service
* **Cytoplasm compartment id**: (default: MNXC3) Compartment ID of interest
* **If there are orphan source species, create creation reaction**: (default: No) NOT SUPPORTED: Upon the merge, if within the heterologous species some are orphaned, then create a reaction that produces that species.
 
### Output

* **rpFBA**: The output is a tar.xz archive containing a list of rpSBML files or a single SBML file

## Algorithm

Three different implementations of constraint based simulation are supported with this tool:
* FBA
* Parsimonious FBA
* Fraction of source reaction: In this method, the flux of a source reaction (For example BIOMASS) is calculated using FBA. Thereafter, the maximal and minimal bounds of that reaction is set from a fraction of that obtimum (default is 0.75) and another FBA is performed for a target reaction. In the pipeline the later would be the reaction that produces the target molecule of interest.

## Installing

To build the image using the Dockerfile, use the following image:

```
docker build -t brsynth/rpfba-rest:dev -f Dockerfile .
```

To run the service on localhost, use the following command:

```
docker run -p 8883:8888 brsynth/rpfba-rest:dev
```

## Prerequisites

* Docker - [Install](https://docs.docker.com/v17.09/engine/installation/)
* libSBML - [Anaconda library](https://anaconda.org/SBMLTeam/python-libsbml)
* cobrapy - [CobaraPy](https://github.com/opencobra/cobrapy)

## Contributing

TODO

## Versioning

Version 0.1

## Authors

* **Melchior du Lac** 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thomas Duigou
* Joan Hérisson

### How to cite rpFBA?

TODO
