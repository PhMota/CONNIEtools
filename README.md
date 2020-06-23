# CONNIEtools
Analysis tools for the CONNIE colaboration

This is a set of tool for analysing the CONNIE 

### Simulation

The `simulate` tool is responsable for generating simulated images. The implemented functionalities and defaults are described by the `--help` as
````sh
$ ./simulation image --help
usage: simulation image [-h] [-g CHARGE_GAIN] [-rn READOUT_NOISE]
                        [-dc DARK_CURRENT] [-exp EXPOSE_HOURS]
                        [-os HORIZONTAL_OVERSCAN] [-vos VERTICAL_OVERSCAN]
                        [--ccd-shape CCD_SHAPE CCD_SHAPE]
                        [--rebin REBIN REBIN] [--image-type IMAGE_TYPE]
                        [--image-mode IMAGE_MODE]
                        [--depth-range DEPTH_RANGE DEPTH_RANGE]
                        [--diffusion-function DIFFUSION_FUNCTION]
                        [--charge-efficiency-function CHARGE_EFFICIENCY_FUNCTION]
                        [-N NUMBER_OF_CHARGES]
                        [--charge-range CHARGE_RANGE CHARGE_RANGE]
                        [--number-of-Cu-charges NUMBER_OF_CU_CHARGES]
                        [--number-of-Cu2-charges NUMBER_OF_CU2_CHARGES]
                        [--number-of-Si-charges NUMBER_OF_SI_CHARGES]
                        [--vertical-modulation-function VERTICAL_MODULATION_FUNCTION]
                        [--horizontal-modulation-function HORIZONTAL_MODULATION_FUNCTION]
                        [--default-vertical-modulation DEFAULT_VERTICAL_MODULATION]
                        [--default-horizontal-modulation DEFAULT_HORIZONTAL_MODULATION]
                        [--default-modulation DEFAULT_MODULATION] [--no-fits]
                        [--pdf] [--spectrum] [--verbose VERBOSE] [--csv]
                        name

positional arguments:
  name                  factor to convert charges into ADU

optional arguments:
  -h, --help            show this help message and exit

parameter options:
  -g CHARGE_GAIN, --charge-gain CHARGE_GAIN
                        factor to convert charges into ADU (default: 7.25)
  -rn READOUT_NOISE, --readout-noise READOUT_NOISE
                        sigma of the normal noise distribution in ADU
                        (default: 0)
  -dc DARK_CURRENT, --dark-current DARK_CURRENT
                        lambda of Poisson distribution in 1/(e-·h) (default:
                        0)
  -exp EXPOSE_HOURS, --expose-hours EXPOSE_HOURS
                        exposed hours (default: 1)

geometry options: 
``` 
