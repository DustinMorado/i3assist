#!/bin/bash

# Example of how to fill the maps directory from the guide.

echo
echo "########################################################################"
echo
echo "Here is a multi-line shell command for the Bash Shell that fills in the"
echo "maps directory inside an I3 project directory."
echo
echo "It assumes you are in the 'maps' directory and that the particle volumes"
echo "are in the directory above your I3 project directory. It also assumes"
echo "that your particles are in the 'mrc' format, make the obvious changes if"
echo "your particles are in 'img' or 'em' format."
echo
echo "------------------------------------------------------------------------"
echo
echo "The command is:"
echo 'i=1'
echo 'for j in ../../*.mrc'
echo 'do'
echo '  ln -sv ${i} p$(printf "%05d" ${i}).mrc'
echo '  i=$((i+1))'
echo 'done'
echo 'unset i'
echo
echo "########################################################################"
echo
