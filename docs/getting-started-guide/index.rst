.. currentmodule:: cratermaker

################
Getting Started
################

As Cratermaker is still early in its development, installing it is not yet straightforward. To begin, you will need to install 
conda, as ``mpas_tools`` is currently only available via a conda package. Once conda is installed, you can create a new conda 
environment from the provided ``environment.yml`` file::

    conda env create -f environment.yml

This will create a new environment called ``cratermaker``. Activate this environment with::
   
      source activate cratermaker

You can then install cratermaker with pip::
   
      pip install -e .

You can verify the installation was successful by running the tests::

      pytest tests/



