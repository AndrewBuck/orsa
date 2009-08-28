#!/bin/bash
for arg in "$@" 
do
  grep -i step $arg | grep -v nan | grep -v prob > "$arg".dat
done
