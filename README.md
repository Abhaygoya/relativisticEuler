relativisticEuler
=================

Computational Physics Final Project; Relativistic Hydrodynamics

Final project for class.

I chose to advance my Euler code to account for relativistic factors, because it's an interesting application of hydrodynamics to astrophysics.

All the equations change, but the main difficulty is that we need to account for the variables being coupled due to the Lorentz factor.
This creates problems in solving for the primitive variables.

Requires use of Newton-Rapheson to determine prim vars to some arbitrary precision. 

2 problems were solved:

1) 1d Reimann with no transverse velocity. This is just a simple shock tube but accounting for relativistic effects.

2) 1d Reimann with a y velocity. This wouldn't really change anything in a non-relativistic system, but due to Lorentz
factor coupling and universal speed limit (light speed) this has a HUGE effect. The initial conditions I chose are particularly
easy to resolve at low resolution. 

3) Higher res version of #2 that, as predicted, can't be resolved at this resolution. Would need access to better computer.
Also considering attempting the adaptive mesh MacFadyen used to make this doable. Might be very difficult to implement. 
