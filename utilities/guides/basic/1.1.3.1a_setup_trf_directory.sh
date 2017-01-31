
#!/bin/bash

# Example of how to fill the trf directory from the guide.

echo
echo "########################################################################"
echo
echo "Here is a multi-line shell command for the Bash Shell that fills in the"
echo "trf directory inside an I3 project directory."
echo
echo "It assumes you are in the 'trf' directory and that the particle"
echo "positions are in the directory above your I3 project directory. It"
echo "assumes that positions are tomogram coordinates describing the center"
echo "point of each particle (i.e. one position line per particle), and that"
echo "you have no orientation information for the particles. It also assumes"
echo "that your particles are in the directory above your I3 project directory"
echo "and that your particles are in the 'mrc' format, make the obvious"
echo "changes if your particles are in 'img' or 'em' format."
echo
echo "------------------------------------------------------------------------"
echo
echo "The command is:"
echo 'i=1'
echo 'for j in ../../*.mrc'
echo 'do'
echo '  eval $(i3stat -sh -o ${j})'
echo '  tx=$(((ox + nx) / 2))'
echo '  ty=$(((oy + ny) / 2))'
echo '  tz=$(((oz + nz) / 2))'
echo '  fmt_i=$(printf "%05d" ${i})'
echo '  echo -n "p${fmt_i} ${tx} ${ty} ${tz} " > p${fmt_i}.trf'
echo '  echo -n "0.0 0.0 0.0 " >> p${fmt_i}.trf'
echo '  echo "1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0" >> p${fmt_i}.trf'
echo '  i=$((i+1))'
echo 'done'
echo 'unset i'
echo
echo "########################################################################"
echo
