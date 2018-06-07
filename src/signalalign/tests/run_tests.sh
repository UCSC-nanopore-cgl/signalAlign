#!/usr/bin/env bash
# This simply executes every python file in the script's current directory
# pytest was having trouble with the cython wrappers

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


for entry in "$DIR"/*.py
do
  echo `python $entry`
done
