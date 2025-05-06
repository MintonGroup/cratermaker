.. currentmodule:: cratermaker

Target
======

The Target class contains basic information about a celestial body on which craters are emplaced. It contains information about the body's size, material properties, and surface gravity. To create a standalone target body, you use the :func:`cratermaker.Target.maker` method. 

.. ipython:: python
    :okwarning:

    from cratermaker import Target
    target = Target.maker("Mars")
    print(target)
