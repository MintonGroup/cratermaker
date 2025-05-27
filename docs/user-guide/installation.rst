.. currentmodule:: cratermaker

.. _ug-installation:

################################
Getting Started with Cratermaker
################################

Cratermaker is designed to be both user-friendly and highly configurable, ensuring that users can start simulations quickly with reasonable defaults or tailor their simulations to specific needs. It is built to be highly modular, and so many of the core components of Cratermaker can be run independently, outside of its core landscape evolution simulation.

Normal installation
===================

Cratermaker is a Python package that is designed to work with Python 3.10 or later on MacOS, Linux, and Windows. It contains some Rust components for performance. You may need to install Rust before installing Cratermaker. Installation instructions for your operating system can be found in `the Rust documentation <https://www.rust-lang.org/tools/install>`__. 

You can install Cratermaker using pip into the Python environment of your choice, provided it is at least Python 3.10 or above.

.. code-block:: bash 

      pip install cratermaker

.. note::

    Cratermaker is undergoing active development, as it is still in alpha stage. You may need to update your installation frequently to get the latest features and bug fixes. You can do this using the following command:
    .. code-block:: bash

        pip install --upgrade cratermaker


Installing Cratermaker from Source
==================================

If you plan to contribute to Cratermaker, or want to try out the latest features that have not yet been released, you can install Cratermaker from source. You can find the latest source code from our `GitHub repository <https://github.com/MintonGroup/cratermaker>`__.  

The first time you install Cratermaker from source, you will need to generate the version file. Due to limitations of the Maturin build system, this must be done manually. You do this by running the following command in the root directory of the Cratermaker source code:

.. code-block:: bash

    python buildscripts/set_version.py    

For development, it's best to build and install the package in "editable" mode. To build and install an editable version you can install Cratermaker from source using the following command:

.. code-block:: bash

    pip install -e .      

Building documentation pages
----------------------------

The documentation pages are built using Sphinx. If you wish to build local versions of the documentation for testing, you will need to install additional dependencies that are not installed by default. We maintain a two requirements files for this purpose, one for pip, and the other for conda, which is used by the Read the Docs build system. You can install the documentation dependencies using the following command:

.. code-block:: bash

    pip install -r environments/requirements-docdev.txt

or, you can generate a conda environment using the conda requirements file:

.. code-block:: bash

    conda env create -f environments/cratermaker-docs.yml

Once you have installed the documentation dependencies, be sure that you also have installed Cratermaker as an editable install into the same environment. Then you can build the documentation using the following command:

.. code-block:: bash

    cd docs
    make clean
    make rtdhtml


The documentation pages will be built in the ``docs/_build/html`` directory, which you can then open in your web browser to view.

Automatically rebuild Rust components
-------------------------------------

The build backend for the Rust components is Maturin. Maturin provides an import hook that can allow the Rust libraries to be automatically recompiled if they are changed, which is useful when Cratermaker installed in editable configuration. To use this, first install the package:

.. code-block:: bash

    pip install maturin_import_hook

In the same virtual environment that you installed Cratermaker in, activate the import hook:

.. code-block:: bash

    python -m maturin_import_hook site install

For more information, see the the `Import Hook <https://www.maturin.rs/import_hook.html>`__ section of the Maturin documentation.
