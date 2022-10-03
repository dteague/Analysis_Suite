.. Analysis Suite documentation master file, created by
   sphinx-quickstart on Fri Jan 29 08:21:39 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree** directive.

Main Page
=========

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   self      
   running
   analysis_suite


******************
Installation Guide
******************

This area requires python3 for all python scripts as well as ROOT, so you need to be in a CMSSW area with ROOT 6.20 or greater. At time of writing, ``CMSSW_11_2_0_pre9`` fits the bill. You will also have to have certain packages not default for

.. code-block:: bash

   # Generic Area Setup
   cmsrel CMSSW_12_3_1
   cd CMSSW_12_3_1/src
   cmsenv
   git clone https://github.com/dteague/Analysis_suite analysis_suite
   cd analysis_suite

   # Setting up scale factors
   git submodule init
   git submodule update
   cp -r jsonpog-integration/POG data

   # Building
   python3 -m pip install --user -r requirements.txt
   scram b -j 4

**Note**: After you scram, the python will compile, but this is not needed, it is also compiled at run time. Scram just sets up the links, so scram is not needed whenever python code is changed after the first run. HOWEVER, scram is needed to compile the cpp code!

************
Quick Start
************

There are two main scripts that can be run in the code:

- ``analyze.py`` for skimming NanoAOD files and turning them in ntuples
- ``run_suite.py`` for modifying and changing the ntuples

For ``analyze.py``, you must run on an xrootd accessed file or to run locally, add the ``--local`` flag. The analysis is the TSelector that one wants the the file to be run with. In hte help text, a list oof the choices are given for convience. A typical example of how the code is run is

.. code-block:: bash

   ./analyze.py -a ThreeTop -v 1 -i <remote file>

Which creates a file called ``output.root``, or if one wants to run over a local file:

.. code-block:: bash

   ./analyze.py -a ThreeTop -v 1 -i <local file> --local -o outName.root

As for the ``run_suite.py``, it takes the output ntuple created by ``analyze.py`` or the farmout jobs. The mode are

- flatten: Turn jagged ntuple into a rectilinear dataframe. Useful for machine learning
- mva: Run the flatten ntuples through one of the many machine learning tools. Returns the same flattened ntuple with the discriminate added and the overall dataframe shortened for the test/train sets needed
- plot: Plots both ntuples and flattened ntuples
- combine: Take input file and creates histograms and cards to be used by combine

A more detailed account of how each module works is given in the running section of this documentation along with a full analysis runthrough in the workflow section. An example of how this code works is:

.. code-block:: bash

   ./run_suite.py <mode> -d <Work Directory> -y <year>

Or if you need more information, just run ``./run_suite.py --help``
