# Dragonfish

Dragonsfish is a functional and taxonomic meta-omic profiler that fully
leverages more recent and signficant advances in algorithms and data structures
for indexing and querying large-scale genomics data. Dragonfish uses
[Pufferfish](https://github.com/COMBINE-lab/pufferfish) at its core, and
therefore a compacted colored de-Bruijn Graph (ccdBG) representation and index.
Pufferfish can efficiently index large collections of genome and transcriptome
sequences and perform fast and accurate alignment of sequencing reads from
genomic and transcriptomic studies, including metagenomic and metatranscriptomic
studies.

Dragonfish efficiently and accurately quantifies the abundance of aligned reads
to annotated functional features shared across all the genomes and
transcriptomes in the index, while also providing strain-, species-, and
genus-level taxonomic abundance contributions to functional features. Dragonfish
also jointly quantifies taxonomic abundances from mappded reads using
[Cedar](https://github.com/COMBINE-lab/pufferfish).

By utilizing modern computational methods for indexing and querying large-scale
genomic data, we now make it feasible to efficiently index on the order of ~100k
full genomic or transcriptomic references, rapidly align reads against an index
of such references, efficiently and accurately quantify their abundances, and map
and summarize abundances to any annotated functional database features as well
as to taxonomy. By accomplishing this, Dragonfish offers significant advantages
to existing functional meta-omic profiling tools that are available.

## Installation

### Mambaforge

Install and set up
[Mambaforge](https://github.com/conda-forge/miniforge#mambaforge), e.g. on
Linux:

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

Obtain the latest Dragonfish project source code:

```bash
git clone https://github.com/hermidalc/dragonfish.git
cd dragonfish
```

## Run Dragonfish

Create and activate the base Dragonfish conda environment, which only
provides snakemake. All the rest of the dependencies are automatically
provided via snakemake and conda when running Dragonfish.

```
mamba env create -f envs/dragonfish.yaml
mamba activate dragonfish
```

Run the full pipeline:

```
snakemake --use-conda --printshellcmds --cores all --scheduler greedy --resources gencode_download_jobs=2
```

## Additional Notes

### Building Pufferfish

**I already provide statically built Pufferfish binaries for Linux x86-64.**
Though here are my notes for trying to build it for other platform and
architectures.

I could not get Pufferfish to build within a conda envrionment using the
required conda dependencies already installed into that environment. The main
issue appears to be with the [SeqLib](https://github.com/walaj/SeqLib)
dependency failing to see its required conda dependency C files or having other
compilation issues.

Therefore, if you want to make you own Pufferfish build, you have to use
dependencies from your system-wide package manager (and therefore potentially
need sudo/root permissions), e.g. for RHEL/Fedora Linux `dnf` or Ubuntu `apt`.
SeqLib successfully builds this way.

### Static build

Get the Dragonfish project Pufferfish submodule:

```bash
git submodule update --init --recursive
```

A static build will create fully self-contained binaries for the OS platform
and CPU architecture you are using to make the build. These binaries can then
simply be copied onto another of the same platform and architecture and will
work. On RHEL/Fedora, install the following dependencies using your package
manager:

```bash
sudo dnf install \
cmake \
curl \
bzip2-static \
gcc \
gcc-c++ \
git \
glibc-static \
hwloc-devel \
jsoncpp-devel \
wget \
xz-static \
zlib-static
```

Download and build [jemalloc](https://github.com/jemalloc/jemalloc):

```bash
cd external
https://github.com/jemalloc/jemalloc/archive/refs/tags/5.3.0.zip
unzip 5.3.0.zip
cd jemalloc-5.3.0

./autogen.sh
make
make install
```

You will find the static library under `lib/libjemalloc.a` and header under
`include/jemalloc/jemalloc.h`. Now go to the Dragonfish project directory and
create the Pufferfish build directory:

```
mkdir external/pufferfish/build
cd external/pufferfish/build
```

If your system has CMake >=3.24, then you can build with the following command:

```bash
cmake \
-DZLIB_USE_STATIC_LIBS \
-DBZip2_USE_STATIC_LIBS \
-DLZMA_USE_STATIC_LIBS \
-DJEMALLOC_LIBRARY=<path to jemalloc>/lib/libjemalloc.a \
-DJEMALLOC_INCLUDE_DIR=<path to jemalloc>/include \
../
```

Otherwise you have to do a temporary bit of hacking to prevent the build from
including your system compression library dynamic files. In RHEL/Fedora run
the following commands:

```bash
sudo mv /usr/lib64/libz.so /usr/lib64/libz.so.ORIG
sudo mv /usr/lib64/libbz2.so /usr/lib64/libbz2.so.ORIG
sudo mv /usr/lib64/liblzma.so /usr/lib64/liblzma.so.ORIG

cmake \
-DJEMALLOC_LIBRARY=<path to jemalloc>/lib/libjemalloc.a \
-DJEMALLOC_INCLUDE_DIR=<path to jemalloc>/include \
../

sudo mv /usr/lib64/liblzma.so.ORIG /usr/lib64/liblzma.so
sudo mv /usr/lib64/libbz2.so.ORIG /usr/lib64/libbz2.so
sudo mv /usr/lib64/libz.so.ORIG /usr/lib64/libz.so
```

Now build with:

```
make
```

After the build successfully completes you now have the statically compiled
Pufferfish binaries under the build `src` directory. Check that they are
indeed static binaries, for example on Linux run `ldd` and you should not see
any dynamic links to any `.so` shared objects or linux kernel dynamic linked
files:

```bash
$ ldd src/*
	not a dynamic executable
src/cedar:
	not a dynamic executable
src/filtersam:
	not a dynamic executable
src/getLineage:
	not a dynamic executable
src/pufferfish:
	not a dynamic executable
```

Again, these are fully self-contained and can simply be copied and used on any
computer with the same OS platform and CPU architecture.

### Dynamic build

In RHEL/Fedora, install the following dependencies:

```bash
sudo dnf install \
cmake \
curl \
bzip2-devel \
gcc \
gcc-c++ \
git \
glibc-devel \
hwloc-devel \
jsoncpp-devel \
jemalloc-devel \
wget \
xz-devel \
zlib-devel
```

Instead of using the Dragonfish forked Pufferfish submodule (where I made
changes to CMakeLists.txt for static build), get the COMBINE-lab official
repository:

```
git clone git@github.com:COMBINE-lab/pufferfish.git
```

Build Pufferfish:

```bash
mkdir pufferfish/build
cd pufferfish/build
cmake ../
make
```
