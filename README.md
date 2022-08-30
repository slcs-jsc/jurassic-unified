# JURASSIC-unified
<img align="middle" src="docu/images/projects.png"  width="550" height="275">

## Short description

* JURASSIC-unified has capabilities to simulate scattering of infrared radiation, but also benefits from GPU acceleration.
* Acceleration of simulation of multiple scattering for radiation transport was done by combining the [JURASSIC-scatter](https://github.com/slcs-jsc/jurassic-scatter) and [JURASSIC-GPU](https://github.com/slcs-jsc/jurassic-gpu) projects. 
The radiative transfer simulation for rays from the lowest recursion level from the scattering simulation is similar to the case without scattering, which was accelerated in the JURASSIC-GPU project.
Therefore, in JURASSIC-unified, radiation for rays from higher recursion levels is calculated on the CPU as in JURASSIC-scatter, and for those from the lowest level on the GPU, leveraging the JURASSIC-GPU implementation.

<img align="middle" src="docu/images/execute.png"  width="500" height="270">

* JURASSIC-unified contains no duplicate code, making the code easy to maintain.
* Because it performs the significant part of the simulation on GPUs, this implementation is around <b>24×</b> faster then the JURASSIC-scatter implementation.

## JURASSIC-unified as a library
* todo..
<img align="middle" src="docu/images/library.png"  width="480" height="480">

## Documentation
* todo..
