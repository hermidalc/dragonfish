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
sequences and millions of their respective CDSs, rapidly align reads against an
index of full references and CDSs, and accurately quantify reads from these
alignments, that this profiler will offer significant advantages to existing
functional meta-omic profiling tools that are available.

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

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
```

After installation, close the terminal and open a new one to get the updated
environment, then:

```bash
conda config --set auto_activate_base false
conda config --set channel_priority strict
mamba update --all
```

### Project source

Obtain the latest released project source code:

```bash
git clone --recurse-submodules https://github.com/hermidalc/dragonfish.git
cd dragonfish
```


### Pufferfish, Puffaligner, and Cedar

I could not get Pufferfish to build within a conda envrionment using the
required dependencies already installed into that environment. The main
issue appears to be with the SeqLib dependency failing to see its required
conda dependency C files or having other compilation issues. Therefore,
you currently have to build Pufferfish using dependencies from your
system-wide package manager (and therefore sudo/root permissions), e.g.
for RHEL/Fedora Linux `dnf` or Ubuntu `apt`. SeqLib successfully builds
this way. In RHEL/Fedora install the following dependencies:

```bash
sudo dnf install \
cmake \
curl \
bzip2-devel \
gcc \
gcc-c++ \
git \
jsoncpp-devel \
jemalloc-devel \
wget \
xz-devel \
zlib-devel
```

Build Pufferfish, Puffaligner, and Cedar:

```bash
mkdir external/pufferfish/build
cd external/pufferfish/build
cmake ../
make
```

Create and activate the base Dragonfish conda environment, which only
provides snakemake. All the rest of the dependencies are automatically
provided via snakemake and conda when running Dragonfish. If using
Miniforge3 replace `mamba` with `conda`:

```
mamba env create -f envs/dragonfish.yaml
mamba activate dragonfish
```

### Running Dragonfish

```
snakemake --use-conda --printshellcmds --cores all
```
