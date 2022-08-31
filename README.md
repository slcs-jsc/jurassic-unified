# JURASSIC-unified

<img align="middle" src="docu/images/projects.png"  width="550" height="275">

## Short description

* JURASSIC-unified has capabilities to simulate scattering of infrared radiation, but also benefits from GPU acceleration.
* Acceleration of simulation of multiple scattering for radiation transport was done by combining the [JURASSIC-scatter](https://github.com/slcs-jsc/jurassic-scatter) and [JURASSIC-GPU](https://github.com/slcs-jsc/jurassic-gpu) projects. 
The radiative transfer simulation for rays from the lowest recursion level from the scattering simulation is similar to the case without scattering, which was accelerated in the JURASSIC-GPU project.
Therefore, in JURASSIC-unified, radiation for rays from higher recursion levels is calculated on the CPU as in JURASSIC-scatter, and for those from the lowest level on the GPU, leveraging the JURASSIC-GPU implementation.

<img align="middle" src="docu/images/execute.png"  width="500" height="270">

* JURASSIC-unified contains no duplicate code, making the code easy to maintain.
* Because it performs the significant part of the simulation on GPUs, this implementation is around **24×** faster then the JURASSIC-scatter implementation.

## JURASSIC-unified as a library

* JURASSIC-unified can be used from the [JURASSIC](https://github.com/slcs-jsc/jurassic) reference project as a library to speed up the simulation of radiation with neglected scattering.
* The `generate_library.sh` script, in which the code is adapted, compiled and wrapped into the C static library `libjurassic unified.a`, must be run first.
This C static library has to be linked and the files from the `include` folder have to be included, when compiling the JURASSIC reference code.
<img align="middle" src="docu/images/library.png"  width="480" height="480">

### Library functions

#### initialization

```c
jur_ctl_t *jur_unified_init(int argc, char *argv[]);
```
Must be called at the beginning of the progrram. This function reads the control parameters and stores them into a static variable. Similarly,
the emissivity look-up table is initialized in this function.

#### formod

```c
void jur_unified_formod_multiple_packages(atm_t const *atm, obs_t *obs, int num_of_obs_packages, int32_t const *atm_id);
```
Since the JURASSIC reference model ignores scattering, this function does not have `aero` among its parameters. Unlike the function from JURASSIC reference, this formod function supports simultaneous calculation of the radiances for possibly multiple observation packages, for rays with possibly different atmospheres.


# Documentation
* To see more details about the project as well as how to use the library to speed up the JURASSIC reference code, please take a look at Stjepan Požgaj's [master's thesis](docu/stjepan_pozgaj_thesis.pdf).
