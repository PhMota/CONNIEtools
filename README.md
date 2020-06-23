# CONNIEtools
Analysis tools for the CONNIE colaboration

This is a set of tool for analysing the CONNIE 

### Simulation

The `simulate` tool is responsable for generating simulated images. The implemented functionalities and defaults are described by the `--help` as
```sh
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
  name                  basename for simulation output

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
  -os HORIZONTAL_OVERSCAN, --horizontal-overscan HORIZONTAL_OVERSCAN
                        size of the horizontal overscan in pixels (default:
                        150)
  -vos VERTICAL_OVERSCAN, --vertical-overscan VERTICAL_OVERSCAN
                        size of the vertical overscan in pixels (default: 90)
  --ccd-shape CCD_SHAPE CCD_SHAPE
                        shape of the image as 2d pixels (default: (4130,
                        4120))
  --rebin REBIN REBIN   2d rebinning strides (default: [1, 1])
  --image-type IMAGE_TYPE
                        image type (default: int)
  --image-mode IMAGE_MODE
                        set to "1" to use official 1x1 image geomtry or "5" to
                        1x5

depth options:
  --depth-range DEPTH_RANGE DEPTH_RANGE
                        range into which to randomly generate depths (default:
                        [0, 670])
  --diffusion-function DIFFUSION_FUNCTION
                        function to map z-depth into xy-sigma (default:
                        sqrt(-258.817238*log1p(-0.000982*z))/15 if z < 670
                        else 0)
  --charge-efficiency-function CHARGE_EFFICIENCY_FUNCTION
                        function for charge efficiency dependent of z-depth
                        (default: 1. if z < 670 else .9)

charge options:   
  -N NUMBER_OF_CHARGES, --number-of-charges NUMBER_OF_CHARGES
                        number of charges to be randomly generated (default:
                        0)
  --charge-range CHARGE_RANGE CHARGE_RANGE
                        range into which to randomly generate charges
                        (default: [5, 200])
  --number-of-Cu-charges NUMBER_OF_CU_CHARGES
                        number of charges to be randomly generated at the
                        Copper fluorescence energy 8.046keV (default: 0)
  --number-of-Cu2-charges NUMBER_OF_CU2_CHARGES
                        number of charges to be randomly generated at the
                        secundary Copper fluorescence energy 8.904keV
                        (default: 0)
  --number-of-Si-charges NUMBER_OF_SI_CHARGES
                        number of charges to be randomly generated at the
                        Silicon fluorescence energy 1.740keV (default: 0)

modulation options:
  --vertical-modulation-function VERTICAL_MODULATION_FUNCTION
                        function to modulate the vertical axis (default: 0)
  --horizontal-modulation-function HORIZONTAL_MODULATION_FUNCTION
                        function to modulate the horizontal axis (default: 0)
  --default-vertical-modulation DEFAULT_VERTICAL_MODULATION
                        set vertical modulation to "50*cos(y/2/pi/20)"
  --default-horizontal-modulation DEFAULT_HORIZONTAL_MODULATION
                        set horizontal modulation to "-1e3*(x/1000. - 1)**2 if
                        x < 1000. else 0"
  --default-modulation DEFAULT_MODULATION
                        set modulations to "-1e3*(x/1000. - 1)**2 if x < 1000.
                        else 0" and "50*cos(y/2/pi/20)"

output options:   
  --no-fits             suppress fits output
  --pdf                 generate pdf output
  --spectrum            generate energy spectrum
  --verbose VERBOSE     verbose level
  --csv                 generate csv output
```



The only necessary option is the `name` which specifies the basename for the output files, the rest of the options will fall to their defaults.
