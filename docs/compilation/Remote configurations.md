# Remote configurations for HPC@UniTN
The following wrapping scripts are broken:
>/apps/mpich-3.2.1--gcc-9.1.0/bin/mpicc

>/apps/mpich-3.2.1--gcc-9.1.0/bin/mpicxx 

To resolve this issue, scripts have been copied* and at the beginning of each script the following lines have been added:
```bash
module load gcc91
```
*due to write permission limitations, the modified scripts are stored in the user's home directory, e.g., `/home/username/bin/mpicc` and `/home/username/bin/mpicxx`.

# Useful tips
Print std versions available with the compiler:
```bash
module load gcc91
gcc-9.1.0 -v --help 2> /dev/null | sed -n '/^ *-std=\([^<][^ ]\+\).*/ {s//\1/p}'
```