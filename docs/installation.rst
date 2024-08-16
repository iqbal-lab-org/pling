Installation
============

Conda
-----

To install with conda, just run::

	conda config --add channels bioconda
	conda config --add channels conda-forge
	conda install pling


From source
-----------

Start with::

	git clone https://github.com/iqbal-lab-org/pling.git && cd pling
	
From here you can install dependancies using conda with the ``env.yaml`` file, or whichever method you prefer. To install dependancies with conda, do::

	conda env create -f env.yaml
	conda activate pling
	
and then finally install pling with::

	python -m pip install .
	
	
Installing optional dependancies
--------------------------------

By default, the DCJ-Indel distances are calculated with the open source solver GLPK. Optionally, pling can use Gurobi for its calculations, but this dependancy is not included in the usual pling installation. This is because Gurobi is a commercial software which requires a license, which is free for academic purposes. Information on obtaining this license can be found here: https://www.gurobi.com/academia/academic-program-and-licenses/. Gurobi can then be installed via conda as follows::

	conda config --add channels https://conda.anaconda.org/gurobi
	conda install gurobi=10.0.1
