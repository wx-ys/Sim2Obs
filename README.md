## Introduction
sim2obs is a Python package designed to prepare simulation data for use with [SKIRT](https://skirt.ugent.be/root/_home.html) and to process the resulting output data.

Note: sim2obs will be included as a submodule of galyst.

## Under Construction

**Code Structure:**
- `simtobs.py`: Main class that integrates all modules and functions.
    - `grabsim`: Handles file I/O for various simulation data formats.
    - `runskirt`: Generates .ski templates and configures .ski file parameters based on the parameter file.
    - ...

**Current Progress:**
- `grabsim`, file I/O (Read_AnastrisTNG): basic functionality completed
- `runskirt`: basic functionality completed
- `mock`, `filter`, `visualization`: planned for future development

## Usage

```python
from sim2obs import simtobs
simio = simtobs.Sim2Obs('SIM2OBS.ini') 
simio.pre_init()
simio.pre_skirt()
```

See [examples](example.ipynb) for more details.

## Maintainers

[@wx-ys](https://github.com/wx-ys)

## License

[MIT](LICENSE) © Shuai Lu

## Acknowledgments

- [SKIRT9](https://github.com/SKIRT/SKIRT9)
- [PTS9](https://github.com/SKIRT/PTS9)
- [cbottrell/SKIRT](https://github.com/cbottrell/SKIRT)
- [galaxyEmulator](https://github.com/xczhou-astro/galaxyEmulator)
