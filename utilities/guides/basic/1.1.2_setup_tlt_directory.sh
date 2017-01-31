#!/bin/bash

# Example of how to fill the tilt angles directory using a general template.

echo
echo "########################################################################"
echo
echo "Here is a multi-line shell command for the Bash Shell that fills in the"
echo "tlt directory inside an I3 project directory."
echo
echo "It assumes you are in the 'tlt' directory and that you have already"
echo "filled in the 'maps' directory in step 1.1.1, and that you have copied"
echo "a 'template.tlt' file into your I3 project directory. It also assumes"
echo "that your particles are in the 'mrc' format, make the obvious changes if"
echo "your particles are in 'img' or 'em' format."
echo
echo "------------------------------------------------------------------------"
echo
echo "The command is:"
echo 'for i in ../maps/*.mrc'
echo 'do'
echo '  ln -sv ../template.tlt $(basename ${i} .mrc).tlt'
echo 'done'
echo
echo "########################################################################"
echo
