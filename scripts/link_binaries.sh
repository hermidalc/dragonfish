#!/bin/bash

cd $(git rev-parse --show-toplevel)
mkdir bin
cd bin
ln -s ../pufferfish/build/src/pufferfish
ln -s ../pufferfish/build/src/cedar
cd ..