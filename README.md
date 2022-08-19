# Dragonfish

A new **functional** metagenomic and metatranscriptomic profiler and toolkit
that leverages recent modern methodological advances in the field of
computational and algorithmic genomics and relevant leading software tools.
Dragonfish uses
[Pufferfish, Puffaligner, and Cedar](https://github.com/COMBINE-lab/pufferfish)
at it's core, written by the [COMBINE-lab](https://github.com/COMBINE-lab), and
therefore a colored and compacted de-Bruijn Graph (ccdBG) indexing data
structure. The underlying data structure and software tools can efficiently
represent, store, and index large collections of genomes and transcriptomes,
and perform fast and accurate alignment of short reads from whole shotgun
metagenomic and metatranscriptomic sequencing experiments to its index.
Dragonfish accurately quantifies the abundance of mapped reads to functional
features across genomes and transcriptomes in the index, while also providing
strain-, species-, and genus-level taxonomic abundance contribution to
functional features.

We believe that utilizing using modern computational genomics methods, that now
make it feasible to efficiently index tens of thousands of full reference
sequences, rapidly align reads against an index of full references, and perform
accurate quantification of reads from these alignments, that this profiler will
offer significant advantages to existing functional meta-omic profiling tools
that are available.

An exclusively taxonomic profiler leveraging the above COMBINE-lab tools
exists called [AGAMEMNON](https://github.com/ivlachos/agamemnon), but
Dragonfish's design goals and functionality are broader and complementary to
that of AGAMEMNON's, since we are focused on functional profiling across
taxonomy and taxonomic contributions to community function.

## Installation

### Mambaforge / Miniforge3

Install and set up
[Mambaforge](https://github.com/conda-forge/miniforge#mambaforge) or
[Miniforge3](https://github.com/conda-forge/miniforge#miniforge3)

### Project source

Obtain the latest released project source code:

```bash
git clone --recurse-submodules https://github.com/hermidalc/dragonfish.git
cd dragonfish
```


### Pufferfish, Puffaligner, and Cedar

I could not get pufferfish to build within a conda envrionment using the
required dependencies already installed into that environment. This appears
to be due to the specifics of how conda-forge built their C++ and related
packages. Therefore, currently you have to install pufferfish using
dependencies from your system-wide package manager, e.g. for RHEL/Fedora
Linux `dnf` or Ubuntu `apt`, and _**with any conda environment
deactivated**_. In RHEL/Fedora install the following dependencies:


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

Create and activate the base Dragonfish conda environment, which only
provides snakemake. All the rest of the dependencies are automatically
provided via snakemake and conda when running Dragonfish. If using Miniforge3
replace `mamba` with `conda`:

```
mamba env create -f envs/dragonfish.yaml
mamba activate dragonfish
```

Run Dragonfish:

```
snakemake --use-conda --printshellcmds --cores all
```