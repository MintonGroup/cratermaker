.. currentmodule:: cratermaker

.. _ug-installation:

################################
Getting Started with Cratermaker
################################

Cratermaker is designed to be both user-friendly and highly configurable, ensuring that users can start simulations quickly with reasonable defaults or tailor their simulations to specific needs. It is built to be highly modular, and so many of the core components of Cratermaker can be run independently, outside of its core landscape evolution simulation.

Cratermaker is a Python package that is designed to work with Python 3.10 or later on MacOS, Linux, and Windows. It contains some Rust components for performance. You may need to install Rust before installing Cratermaker. Installation instructions for your operating system can be found in `the Rust documentation <https://www.rust-lang.org/tools/install>`__. 

Install dependencies
====================

Some components of Cratermaker use OpenBLAS for linear algebra operations. On some systems, you may need to install OpenBLAS separately before installing Cratermaker. Installation instructions for your operating system can be found in `the OpenBLAS documentation <http://www.openmathlib.org/OpenBLAS/docs/install/>`__.  You may also need to install OpenSSL development libraries for your system.

Cratermaker uses a number of packages which must be installed on the system before it can be built. If any are missing, they could cause a failure to build. Try one of the following, depending on what OS you are using:

Debian Linux (e.g. Ubuntu)
--------------------------

.. code-block:: console

    sudo apt-get update 
    sudo apt-get install pkg-config libssl-dev libopenblas-dev

Redhat Linux (e.g. Fedora) with yum
-----------------------------------

.. code-block:: console

    yum install epel-release
    yum install pkgconfig openssl-devel openblas-devel


Redhat Linux (e.g. Fedora) with dnf
-----------------------------------

.. code-block:: console
    
   dnf install pkgconf-pkg-config openssl-devel flexiblas-devel 


MacOS with Homebrew
--------------------

.. code-block:: console

    brew install openssl@3 openblas


Windows with vcpkg
------------------

.. code-block:: console

    vcpkg install openblas intel-mkl --triplet x64-windows

Normal Installation
===================

You can install Cratermaker using pip into the Python environment of your choice, provided it is at least Python 3.10 or above.

.. code-block:: console 

      pip install cratermaker

.. note::

    Cratermaker is undergoing active development, as it is still in alpha stage. You may need to update your installation frequently to get the latest features and bug fixes. You can do this using the following command:

    .. code-block:: console

        pip install --upgrade cratermaker


Installing Cratermaker from Source
==================================

If you plan to contribute to Cratermaker, or want to try out the latest features that have not yet been released, you can install Cratermaker from source. You can find the latest source code from our `GitHub repository <https://github.com/MintonGroup/cratermaker>`__.  

Some components of Cratermaker are built using Rust, and so you will need to have Rust installed on your system to build Cratermaker from source. Installation instructions for your operating system can be found in `the Rust documentation <https://rust-lang.org/learn/get-started/>`__.

The first time you install Cratermaker from source, you will need to generate the version file. Due to limitations of the Maturin build system, this must be done manually. You do this by running the following command in the root directory of the Cratermaker source code:

.. code-block:: console

    python buildscripts/set_version.py    

For development, it's best to build and install the package in "editable" mode. To build and install an editable version you can install Cratermaker from source using the following command:

.. code-block:: console

    pip install -e .      

Building documentation pages
----------------------------

The documentation pages are built using Sphinx. If you wish to build local versions of the documentation for testing, you will need to install additional dependencies that are not installed by default. We maintain a two requirements files for this purpose, one for pip, and the other for conda, which is used by the Read the Docs build system. You can install the documentation dependencies using the following command:

.. code-block:: console

    pip install -r environments/requirements-docdev.txt

or, you can generate a conda environment using the conda requirements file:

.. code-block:: console

    conda env create -f environments/cratermaker-docs.yml

Once you have installed the documentation dependencies, be sure that you also have installed Cratermaker as an editable install into the same environment. Then you can build the documentation using the following command:

.. code-block:: console

    cd docs
    make clean
    make rtdhtml


The documentation pages will be built in the ``docs/_build/html`` directory, which you can then open in your web browser to view.

Automatically rebuild Rust components
-------------------------------------

The build backend for the Rust components is Maturin. Maturin provides an import hook that can allow the Rust libraries to be automatically recompiled if they are changed, which is useful when Cratermaker installed in editable configuration. To use this, first install the package:

.. code-block:: console

    pip install maturin_import_hook

In the same virtual environment that you installed Cratermaker in, activate the import hook:

.. code-block:: console

    python -m maturin_import_hook site install

For more information, see the the `Import Hook <https://www.maturin.rs/import_hook.html>`__ section of the Maturin documentation.

Common problems 
===============

Here are a few things you can try if you have trouble installing the Cratermaker project from source.


Make sure Rust is up-to-date
----------------------------

Be sure your rust install is up-to-date 

.. code-block:: console

    rustup update --nightly


If you see errors related to linalg functions, you may have an old and incompatable linalg crate installed. The most expedient way I've found to deal with this is to remove the Cargo.lock file. You can do this in a clean way with the cargo-unlock crate:

.. code-block:: console

    cargo install cargo-unlock
    cargo unlock

Installing polars on legacy hardware
------------------------------------

If you ever see an error like this::

    illegal hardware instruction (core dumped)

Or a message about polars and incompatible hardware, you may have to install the runtime-compatability version of polars.

.. code-block:: console

    pip install "polars[rtcompat]"


    

