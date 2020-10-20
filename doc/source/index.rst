rpFBA's documentation
=====================

Indices and tables
##################

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Introduction
############

.. _CobraPy: https://opencobra.github.io/cobrapy/
.. _rpBase: https://github.com/Galaxy-SynBioCAD/rpBase
.. _rpCache: https://github.com/Galaxy-SynBioCAD/rpCache
.. _libSBML: http://sbml.org/Software/libSBML

Welcome to the documentation of rpFBA. The project uses the rpBase_ project to merge a rpSBML file and a SBML files to simulate FBA using CobraPy_. The project supports the "fraction of reaction", standard FBA and parsimonious FBA. 

Usage
#####

First build the rpBase_ and rpCache_ dockers before building the local docker:

.. code-block:: bash

   docker build -t brsynth/rpfba-standalone:v2 -f Dockerfile .

To call the docker locally you can use the following command:

.. code-block:: bash

   python run.py -input /path/to/file -input_format tar -gem_sbml /path/to/file -output /path/to/output

API
###

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. currentmodule:: rpToolServe

.. currentmodule:: rpTool

.. autoclass:: rpFBA
    :show-inheritance:
    :members:
    :inherited-members:

.. currentmodule:: run

.. autoclass:: main
    :show-inheritance:
    :members:
    :inherited-members:

