# Dragonfish

A **functional** metagenomic and metatranscriptomic profiler and toolkit
that leverages modern methodological advances in the field of computational
and algorithmic genomics and relevant published software. Dragonfish uses
[Pufferfish + Puffaligner + Cedar](https://github.com/COMBINE-lab/pufferfish)
at it's core, written by the [COMBINE-lab](https://github.com/COMBINE-lab).
These underlying tools can efficiently store, index, and represent large
collections of genomes and transcriptomes, and perform fast and accurate
alignment of short reads from whole shotgun metagenomic and metatranscriptomic
sequencing experiments to an index. Dragonfish accurately quantifies the
abundance of mapped reads to functional features across genomes and
transcriptomes in the index, while also providing genus-, species-, and
strain-level taxonomic abundance contribution to these functional features.
We believe that by using modern computational genomics methods, and a more
accurate alignment-based approach, that this toolkit will offer significant
advantages to existing functional meta-omic profiling tools that are
available.

An exclusively taxonomic profiler leveraging the above COMBINE-lab tools
exists called [AGAMEMNON](https://github.com/ivlachos/agamemnon), but
Dragonfish's design goals and functionality are complementary and different
that AGAMEMNON's, since we are focused primarily on functional profiling.

## Installation

### Mambaforge / Miniforge3

Install and set up
[Mambaforge](https://github.com/conda-forge/miniforge#mambaforge) or
[Miniforge3](https://github.com/conda-forge/miniforge#miniforge3)

### Project source

To obtain the source of the project and create a conda environment
with the tools needed to run the project, execute the following below (if
using Miniforge3 replace `mamba` with `conda`):

```bash
git clone --recurse-submodules https://github.com/hermidalc/dragonfish.git
cd dragonfish
```


### Pufferfish, Puffaligner, and Cedar

I could not get pufferfish to build within a `conda` envronment using the
required dependencies already installed into that environment. This appears
to be due to the specifics of how conda-forge built their C++ and related
packages. Therefore, currently you have to install pufferfish using
dependencies from your system-wide package manager, e.g. for RHEL/Fedora
Linux `dnf` or Ubuntu `apt`, and with any conda environment deactivated.
In RHEL/Fedora install the following:


```bash
sudo dnf install \
autoconf \
boost-devel \
ca-certificates \
cmake \
curl \
bzip2-devel \
gcc \
gcc-c++ \
git \
jemalloc-devel \
make \
time \
wget \
xz-devel \
zlib-devel
```

Install Pufferfish, Puffaligner, and Cedar:

```bash
mkdir external/pufferfish/build
cd external/pufferfish/build
cmake ../
make
```

Create and activate the base Dragonfish `conda` environment (provides
`snakemake`):

```
mamba env create -f envs/dragonfish.yaml
mamba activate dragonfish
```

Run Dragonfish:

```
snakemake --use-conda --printshellcmds --cores all
```