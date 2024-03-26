#!/bin/bash

set -e
set -v

for n in {00..52}
  do
    ln -s /data/anna/animal_metagenome/long-mg-dog/01_assembly/D0${n}/assembly_info.txt D0${n}_assembly_info.txt
  done
