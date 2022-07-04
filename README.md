# Dragonfish

A novel **functional** metagenomic and metatranscriptomic profiler and framework
that leverages more recent methdological advances in the field of computational
and algorithmic genomics and relevant published software, which can efficiently
store, index, and represent large collections of genomes and transcriptomes,
perform fast and accurate alignment of short reads from whole shotgun
metagenomic or metatranscriptomic sequencing experiments to such an index, and
most importantly, accurately quantify the abundance of aligned reads to
**functional** groups across genomes and transcriptomes in the index, while also
providing species- and strain-level contribution to these functional abundances.

Dragonfish uses [Pufferfish + Puffaligner + Cedar](https://github.com/COMBINE-lab/pufferfish)
at it's core, written by the COMBINE-lab. It is inspired by
[AGAMEMNON](https://github.com/ivlachos/agamemnon), but it's design and
functionality goes significantly further.

## Installation

```bash
git clone --recurse-submodules https://github.com/hermidalc/dragonfish.git
cd dragonfish
```

### Pufferfish, Puffaligner, Cedar

I could not get pufferfish to build within a conda envronment using the
required dependencies already installed into that environment. Therefore,
currently you have to install pufferfish using dependencies from your
system-wide package manager, e.g. for RHEL/Fedora Linux `dnf` the following:


```bash
$ sudo dnf install \
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
