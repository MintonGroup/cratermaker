.. currentmodule:: cratermaker

.. image:: ../_static/scaling_icon.svg
    :alt: Production
    :align: center
    :width: 600px


.. _ug-scaling: 

Scaling
=======

Cratermaker's :ref:`Scaling <api-scaling>` component provides tools to convert :ref:`Projectile <ug-projectile>` parameters (e.g., projectile diameter, velocity, angle, density), and :ref:`Target <ug-target>` properties (material type, surface gravity, density) into a final crater diameter. It also categorizes craters by their morphology: **simple**, **transitional**, or **complex**. There are two scaling models available: "montecarlo" and "ctem". The "ctem" model closely resembles the model used in the CTEM code, a Fortran-based ancestor of CTEM. The "montecarlo" model is a newer non-deterministic model that includes the intrinsic variability in projectile to crater size scaling relationships. We have also included updated simple-to-complex transition diameter values from Schenk et al. (2021) that incorporate data from Vesta and Ceres from the Dawn mission.   

Setting up a Scaling object
---------------------------

As discussed in the :ref:`defaults <ug-defaults>` section, all Cratermaker components have a default configuration that is set when no arguments are passed to its ``maker`` method.  We can inspect the defaults for Scaling:

.. ipython:: python
   :okwarning:

   from cratermaker import Scaling
   scaling = Scaling.maker()
   print(scaling)

We can see that default scaling model is "montecarlo" and that the Scaling object contains a :ref:`Target <ug-target>` (which defaults to "Moon") and a :ref:`Projectile <ug-projectile>` (which defaults to "asteroids"). The other properties, such as target and projectile densities, target material, and the scaling parameters, were set based on the specific target and projectile types that were chosen. We can explore what happens when we pass other options to the model:


.. ipython:: python
   :okwarning:

   scaling = Scaling.maker("ctem", target="Mars", projectile="comets")
   print(scaling)


Typically a scaling object would not be used on its own, but instead as an object that gets passed to the ``maker`` method of a :ref:`Crater <ug-crater>`. However, we can still use it on its own and see how it is constructed and used.

You can also override the default materials that are associated with target bodies. For instance the default material for Mars is "Soft Rock", but we could override this with something else, like "sand" or "ice":

.. ipython:: python
   :okwarning:

   print(Scaling.maker(target="Mars", material="sand"))
   print()
   print(Scaling.maker(target="Mars", material="ice"))


You can also adjust individual scaling parameters manually:

.. ipython:: python
    :okwarning:

   print(Scaling.maker(target="Moon", Ybar=1e7))


You are not limited to using pre-defined materials. For instance, you can specify your own materials with unique properties, but you must specify all of the parameters required by the given model. Both scaling models available in Cratermaker use the same set of required parameters (K1, mu, Ybar, and density), and these must all be set. See Richardson (2009) and Holsapple (1993) for details on the scaling model parameters. 

.. ipython:: python
    :okwarning:
    
    print(Scaling.maker(material="Flubber", K1=3.8, mu=0.1, Ybar=1e7, density=1500.0))

Using the Scaling object
------------------------

Once you have a scaling object, you can use it to compute between projectile, transient, and final crater sizes. Any :ref:`Scaling <api-scaling>`  model will come with the following methods:

.. currentmodule:: cratermaker.Scaling


- :meth:`projectile_to_transient`: Takes the projectile diameter and returns the transient crater diameter.
- :meth:`transient_to_projectile`: Takes the transient crater diameter and returns the projectile diameter.
- :meth:`transient_to_final`: Takes the transient crater diameter and returns both the final crater diameter and the morphology type.
- :meth:`final_to_transient`: Takes the final crater diameter and returns the transient crater diameter. Optionally, you can also provide the morphology type, but if you don't provide it, it will be computed (though it might not be the same as the one used to create the final crater!)


.. ipython:: python
    :okwarning:


    scaling = Scaling.maker()
    dt = scaling.projectile_to_transient(1000.0) 
    print(f"1 km projectile -> Transient crater: {dt*1e-3:.2f} km")
    df, mt = scaling.transient_to_final(dt) 
    print(f"1 km projectile -> Final crater: {df*1e-3:.2f} km, Morphology type: {mt}")
    dt, mt = scaling.final_to_transient(df, mt)
    print(f"Final crater -> Transient: {dt*1e-3:.2f} km")
    pd = scaling.transient_to_projectile(dt)
    print(f"Final crater -> projectile: {pd*1e-3:.3f} km")


Because the scaling model is probabilistic, the results will vary slightly each time you run it. This is particualrly true when the crater size is near the simple-to-complex transition. Remember how passing the morphology type to :meth:`final_to_transient` is optional? In the above example we used the morphology type that was returned from the :meth:`transient_to_final` method. If we had not passed it, the :meth:`final_to_transient` method would have computed it again, and it might not be the same as the one used to create the final crater. Here is an example showing this effect:

.. ipython:: python
    :okwarning:

    scaling = Scaling.maker(rng_seed=3098351)
    pd_in = 800.0
    dt = scaling.projectile_to_transient(pd_in) 
    df, mt_in = scaling.transient_to_final(dt) 
    dt, mt_out = scaling.final_to_transient(df)
    pd_out = scaling.transient_to_projectile(dt)
    print(f"Transient -> Final morphology: {mt_in}")
    print(f"Final -> Transient morphology: {mt_out}")
    print(f"Input projectile diameter  : {pd_in*1e-3:.3f} km")
    print(f"Output projectile diameter : {pd_out*1e-3:.3f} km")


More Scaling examples
---------------------

See more examples at  :ref:`gal-scaling`

