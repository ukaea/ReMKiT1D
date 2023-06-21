Title: Running
Author: Stefan Mijin
Date: 16.05.2023.

If the build process is successful, the executable will be in build/src/executables/ReMKiT1D.
ReMKiT1D requires a config.json file as input. Assuming a valid config file is available, the code can be run from the executable directory using, for example:
```
mpirun -np [num_procs] ./ReMKiT1D
```
where [num_procs] should be set to the desired number of processes. An alternative path to the config file can be specified using the command line option 

```
-with_config_path=/path/to/config/file.json
```

For more details on how to use the code together with the Python support modules, see [Python repos](ADD LINK)