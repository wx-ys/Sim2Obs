## Introduction
sim2obs is a python package for preparing simulation data for [SKIRT](https://skirt.ugent.be/root/_home.html) and processing the output data.

## Under Construction

Code Structure:
simtobs.py: a class that collects all modules and functions.
* grabsim: file IO for processing different simulation data.
* runskirt: according to the parameters file, provide the .ski template and set the .ski file parameters.
* ...

Current progress:
* grabsim, one file IO: (Read_AnastrisTNG), roughly completed
* runskirt, roughly completed
* mock, filter, visualization ...  future work


## Usage


```python
from sim2obs import simtobs
simio = simtobs.Sim2Obs('SIM2OBS.ini') 
simio.pre_init()
simio.pre_skirt()

```

See [examples](example.ipynb) for more details,

## Maintainers

[@wx-ys](https://github.com/wx-ys).


## License

[MIT](LICENSE) © Shuai Lu

## Acknowledgments
* [SKIRT9](https://github.com/SKIRT/SKIRT9)
* [PTS9](https://github.com/SKIRT/PTS9)
* [https://github.com/cbottrell/SKIRT](https://github.com/cbottrell/SKIRT)
* [galaxyEmulator](https://github.com/xczhou-astro/galaxyEmulator)