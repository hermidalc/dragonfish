# Functional metagenomic/metatranscriptomic profiler

A novel functional metagenomic and metatranscriptomic profiler and framework
that leverages more recent methdological advances in the field of computational
and algorithmic genomics and relevant published software, which can efficiently
store, index, and represent large collections of genomes and transcriptomes,
perform fast and accurate alignment of short reads from whole shotgun
metagenomic or metatranscriptomic sequencing experiments to such an index, and
accurately quantify the abundance of aligned reads to functional categories
while also providing species- and strain-level contribution to functional
abundances.

## Installation

### Pufferfish

I could not get pufferfish to build within a conda envronment using the
required dependencies already installed into that environment. Therefore,
currently you have to install pufferfish using dependencies from your
system-wide package manager, e.g. for RHEL/Fedora Linux `dnf` the following:


```bash
$ sudo dnf install autoconf boost-devel ca-certificates cmake curl bzip2-devel gcc gcc-c++ git jemalloc-devel make time wget xz-devel zlib-devel
```
