.. currentmodule:: cratermaker

.. _ug-installation:

################################
Getting Started with Cratermaker
################################

Cratermaker is designed to be both user-friendly and highly configurable, ensuring that users can start simulations quickly with reasonable defaults or tailor their simulations to specific needs. It is built to be highly modular, and so many of the core components of Cratermaker can be run independently, outside of its core landscape evolution simulation.

Normal installation
===================
To begin, simply install Cratermaker using pip into the Python environment of your choice, provided it is at least Python 3.10 or above.

.. code-block:: bash 

      pip install cratermaker

.. note::

    Cratermaker is undergoing active development, as it is still in alpha stage. You may need to update your installation frequently to get the latest features and bug fixes. You can do this using the following command:
    .. code-block:: bash

        pip install --upgrade cratermaker


Installing Cratermaker from Source
==================================

Cratermaker is mostly Python, but it does have some Rust components. To install Cratermaker from source, you will need to have Rust installed. You can install Rust using the following command:

.. code-block:: bash

    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

The build and install an editable version you can install Cratermaker from source using the following command:

.. code-block:: bash

    pip install -e .      

The build backend for the Rust components is Maturin. Maturin provides an import hook that can allow the Rust libraries to be automatically recompiled if they are changed, which is useful when Cratermaker installed in editable configuration. To use this, first install the package:

.. code-block:: bash

    pip install maturin_import_hook

In the same virtual environment that you installed Cratermaker in, activate the import hook:

.. code-block:: bash

    python -m maturin_import_hook site install

For mor information, see the the `Import Hook <https://www.maturin.rs/import_hook.html>`__ section of the Maturin documentation.



