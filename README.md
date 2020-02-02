# REST rpFBA

REST tool that reads a collection of rpSBML files (in a tar.xz) and a GEM SBML model, merges the the two to simulate the introduction of a heterologous pathway to the organism of choice. FBA is then performed with all defined objectives defined by the [FBC](https://co.mbine.org/specifications/sbml.level-3.version-1.fbc.version-2.release-1) SBML package and saved inside custom (IBISBA) annotations. The output is also a collection of rpSBML files in a tar.xz. The merged version of the model may be kept instead of the heterologous pathway alone (default).


### Prerequisites

* Docker - [Install](https://docs.docker.com/v17.09/engine/installation/)
* libSBML - [Anaconda library](https://anaconda.org/SBMLTeam/python-libsbml)
* cobrapy - [CobaraPy](https://github.com/opencobra/cobrapy)

### Compiling and running

```
docker build -t brsynth/rpfba-rest -f Dockerfile .
```

Run the service

```
docker run -p 8884:8888 brsynth/rpfba-rest
```

## Running the tests

TODO

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Galaxy](https://galaxyproject.org) - The Galaxy project

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

TODO

## Authors

* **Melchior du Lac** 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thomas Duigou
* Joan Hérisson
