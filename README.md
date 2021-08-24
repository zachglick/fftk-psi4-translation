# fftk-psi4-translation
Set of Gaussian input/output files generated/parsed by FFTK and their corresponding Psi4 translations

## Organization
The directory `PRLD` contains the original Gaussian inputs/outputs from FFTK.

The directory `psi4` contains psi4 inputs/outputs that compute the same quantities as in `PRLD`

## Dependencies
To install the psi4 program, follow the instructions [here](https://psicode.org/installs/v14/).

Depending on the version of psi4 installed, you may have to separately install the psi4 RESP plugin as explained [here](https://github.com/cdsgroup/resp)

Lastly, these calculations require the developer version of the optking geometry optimizer. In the future, this geometry optimizer will be bundled with psi4. For now, FFTK users must separately download the optimizer [here](https://github.com/psi-rking/optking).

## Calculation details

### psi4/charge-opt-resp

The file `prld-resp.py` contains an example of how to calculate RESP charges with psi4. This file can be executed from the command line as:
``` >> python prld-resp.py ```
The RESP charges are printed to the console. This output may be captured (as is shown in `1_default_results.out`) and parsed.

### psi4/geom-opt

The file `PRLD-GeomOpt.dat` contains an example of how to perform a geometry optimization with psi4. This file can be executed from the command line as:
``` >> psi4 -nX PRLD-GeomOpt.dat ``` 
(where X in the number of CPU cores). The output of the geometry optimization (including the geometry and energy of the optimized structure) is automatically written to `PRLD-GeomOpt.dat`







