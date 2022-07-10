# Dragonfish

A novel **functional** metagenomic and metatranscriptomic profiler and toolkit
that leverages modern methdological advances in the field of computational
and algorithmic genomics and relevant published software. It can efficiently
store, index, and represent large collections of genomes and transcriptomes,
perform fast and accurate alignment of short reads from whole shotgun
metagenomic or metatranscriptomic sequencing experiments to an index, and most
importantly, accurately quantify the abundance of mapped reads to functional
groups across genomes and transcriptomes in the index, while also providing
genus-, species-, and strain-level abundance contribution to these functional
groups. We believe that by using modern methods, this toolkit will offer
significant advantages to existing functional profiling tools that are
available, which are based on older methods.

Dragonfish uses [Pufferfish + Puffaligner + Cedar](https://github.com/COMBINE-lab/pufferfish)
at it's core, written by the [COMBINE-lab](https://github.com/COMBINE-lab). A
purely taxonomic profiler called [AGAMEMNON](https://github.com/ivlachos/agamemnon)
exists that was built using the same underlying COMBINE-lab toolkit, but our
design goals and functionality are both complementary and a superset of
AGAMEMNON's, since we are quantifying abundances by functional feature,
taxonomic level contribution to functional feature abundance, and taxonomic
feature abundances.

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

Install the `conda` environment. Only after installing the Pufferfish step
below you can a `activate` it:

```
mamba env create -f envs/dragonfish-mkl.yaml
```

### Pufferfish, Puffaligner, Cedar

I could not get pufferfish to build within a `conda` envronment using the
required dependencies already installed into that environment. This appears
to be due to specifics of to how conda-forge built their C++ and related
packages. Therefore, currently you have to install pufferfish using
dependencies from your system-wide package manager, e.g. for RHEL/Fedora
Linux `dnf` or Ubuntu `apt`, and with the `mamba deactivate` because will
not work inside the `conda` environment. In RHEL/Fedora install the
following:


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
