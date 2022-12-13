# Simulation of D resonance decays
Repository with code to perform simulations of D resonance decays with [PYTHIA8](https://pythia.org/) and [EvtGen](https://evtgen.hepforge.org/)

## [EvtGen](https://evtgen.hepforge.org/) installation
To install [EvtGen](https://evtgen.hepforge.org/) you can either follow the instructions in the website https://evtgen.hepforge.org/, or use [aliBuild](https://alisw.github.io/alibuild/):
```shell
aliBuild build EVTGEN --defaults o2
```

## Load EVTGEN environment
```shell
alienv enter EVTGEN/latest
```

## Run script
To run the script use:

```shell
root -l -b simulateDresoDecays.cc"(nEvents, pdgCode, decayer, seed, trigger)"
```

where:

- `nEvents = int` is the number of decays to be simulated
- `pdgCode = int` is the pdg code of the resonance (either `435` or `10433`)
- `decayer = int` is the decayer to be used (either `kPythia8` or `kEvtGen`)
- `seed = int` is the seed to be used for reproducibility
- `trigger = string` is the riigger class for the transverse momentum shape (either `\"MB\"` or `\"HB\"`)
