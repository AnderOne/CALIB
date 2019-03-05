# CALIB (Contour Advection Library)

# Dependencies

GCC (version >= 5.5.0) or other compiler supporting C++14.

Blitz++:

https://github.com/blitzpp/blitz

fftw (version >= 3.3.6):

http://www.fftw.org/

...

# Synopsis

The Contour Advection Method is one of the Lagrangian approaches to simulation of scalar field transport in quasi-two-dimensional inviscid incompressible flows that are considered in numerous applications of geophysical fluid dynamics.
The main idea is to use the contours of scalar fields as fluid elements transferred by flow.
This approach allows to represent fine-scale structures of the flow in numerical models.
The algorithm also includes contour editing procedure ("surgery") for conserving computational efficiency, when tracer field structures become very complex.

The Contour Advective semi-Lagrangian method (CASL) is a hybrid of the Eulerian and Lagrangian approaches to simulate of two-dimension flows of inviscid incompressible fluid.
It's based on other Lagrangian approach call as Contour Dynamics which found wide use in simulation of geophysical flows.
Here the contours of relative or potential vorticity are used as fluid elements whose distribution determines the velocities responsible for their transfer.

The Contour Dynamics uses inversion procedure based on Green's function calculate.
The CASL proposes another inversion procedure that includes conversion from contours to Eulerian grid and calculate velocity by FFT-based solver of Poisson equation.

This library implemets simple version of the CASL algorithm to single-layer shallow-water model.
Such model is governed of the potential vorticity equation:  
...

More details:  
[...](doc/slides.pdf)

# References

* Zabusky N.J., Hughes M.H., Roberts K.V. Contour Dynamics for the Euler equations in two dimensions // Journal of Computational Physics. 1979, Vol.30, P.96-106.

* Kozlov V.F. The Contour Dynamics Method in Model Problems of the Ocean Topographic Cyclogenesis // Izv. Akad. Nauk SSSR, Fiz. Atmos. Okean. 1983, Vol.19, N.8, P.635-640.

* Dritschel D.G. Contour Dynamics and Contour Surgery: Numerical algorithms for extended, high-resolution modelling of vortex dynamics in two-dimensional, inviscid, incompressible flows // Computer Physics Reports. 1989, Vol.10, N.3, P.77-146.

* Makarov V. G. The Numerical Algorithm of the Contour Dynamics Method with Varying Topology of Investigated Domains // Model. Mekh. 1991, Vol.5, N.4, P.83-95.

* Waugh D.W., Plumb R.A. Contour Advection with Surgery: A Technique for Investigating Finescale Structure in Tracer Transports // Journal of the Atmospheric Sciences. 1994, Vol.51, N.4, P.530-540.

* Dritschel D.G., Ambaum M.H.P. A Contour-Advective Semi-Lagrangian Numerical Algorithm for simulating fine-scale conservative dynamical fields // Quarterly Journal of the Royal Meteorological Society. 1997, Vol.123, N.540, P.1097-1130.

* Dritschel D.G., Ambaum M.H.P. The Diabatic Contour Advective Semi-Lagrangian model // Monthly Weather Review. 2006, Vol.134, N.9, P.2503-2514.

* Fontane J., Dritschel D.G. The HyperCASL algorithm: A new approach to the numerical simulation of geophysical flows // Journal of Computational Physics. 2009, Vol.228, N.17, P.6411-6425.

* Baranov A.A., Permyakov M.S. An Accelerated Topology Change Algorithm for the Contour Advection Method // Numerical Methods and Programming. 2013, Vol.14, P.75--87.

* Baranov A.A., Permyakov M.S. Analysis of accuracy and computational efficiency of the Contour Advection Method for the barotropic vorticity equation // Numerical Methods and Programming. 2014, Vol.15, P.337-350.

* Baranov A.A., Permyakov M.S. A Contour-Advective semi-Lagrangian numerical algorithm for the problem of interaction between a vortex and an isolated topographic feature on a β-plane // Numerical Methods and Programming. 2014, Vol.15, P.621-630.
