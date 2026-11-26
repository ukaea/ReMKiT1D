Title: Running
Author: Stefan Mijin
Date: 2026

If the build process is successful, the executable will be in `build/src/executables/ReMKiT1D/ReMKiT1D` (or `<prefix>/bin/ReMKiT1D` if installed via `cmake --install`).

ReMKiT1D requires a `config.json` file as input. Assuming a valid config file is available in the current working directory, the code can be run using, for example:

```bash
mpirun -np [num_procs] ./build/src/executables/ReMKiT1D/ReMKiT1D
```

where `[num_procs]` should be set to the desired number of MPI processes.

### Command-line options

- Specify an alternative path to the configuration file:

```bash
mpirun -np [num_procs] ./ReMKiT1D -with_config_path=/path/to/config.json
```

- Print the ReMKiT1D version string and exit:

```bash
./ReMKiT1D --version
# or
./ReMKiT1D -v
```

For more details on how to configure runs and analyse outputs with the Python support modules, see [ReMKiT1D-Python](https://github.com/ukaea/ReMKiT1D-Python).
