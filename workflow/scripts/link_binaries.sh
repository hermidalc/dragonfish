#!/bin/bash

cd $(git rev-parse --show-toplevel)
mkdir bin
cd bin
ln -s ../external/pufferfish/build/src/pufferfish
ln -s ../external/pufferfish/build/src/cedar
cd ..
