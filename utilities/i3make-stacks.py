from __future__ import print_function, division
import glob
import os
import os.path
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--suffix', '-s', help='map suffix(default=mrc)',
                    metavar='SUFFIX', default='mrc')
parser.add_argument('--directory', '-d', metavar='DIRECTORY',
                    help='directory with maps(default=$PWD)',
                    default=os.getcwd())
parser.add_argument('--particle_prefix', '-p', metavar='PRFX',
                    help='New stack maps names(default=p)',
                    default='p')
parser.add_argument('--script', '-f', metavar='SCRIPT',
                    help='i3concat script name(default=make_stacks.sh)',
                    default='make_stacks.sh')

args = parser.parse_args()

glob_regex = os.path.join(args.directory, "*." + args.suffix)
ptcl_list = sorted(glob.glob(glob_regex))
test_ptcl = ptcl_list[0]

i3stat_output = subprocess.check_output(["i3stat", "-sh", "-o", test_ptcl])
eval(i3stat_output)

test_boxsize_x = nx
test_boxsize_y = ny
test_boxsize_z = nz

test_center_x = (ox + nx) // 2
test_center_y = (oy + ny) // 2
test_center_z = (oz + nz) // 2

for ptcl in ptcl_list[1:]:
    i3stat_output = subprocess.check_output(
        ["i3stat", "-sh", "-o", ptcl])
    eval(i3stat_output)
    try:
        assert(test_boxsize_x == nx)
    except AssertionError:
        print("ERROR: " + ptcl + " has a different X size.")
        raise

    try:
        assert(test_boxsize_y == ny)
    except AssertionError:
        print("ERROR: " + ptcl + " has a different Y size.")
        raise

    try:
        assert(test_boxsize_z == nz)
    except AssertionError:
        print("ERROR: " + ptcl + " has a different Z size.")
        raise

    try:
        assert(test_center_x == ox)
    except AssertionError:
        print("ERROR: " + ptcl + " has a different X origin.")
        raise

    try:
        assert(test_center_y == oy)
    except AssertionError:
        print("ERROR: " + ptcl + " has a different Y origin.")
        raise

    try:
        assert(test_center_z == oz)
    except AssertionError:
        print("ERROR: " + ptcl + " has a different Z origin.")
        raise

trfm_list = [ (os.path.splitext(x)[0] + '.trf') for x in ptcl_list ]
stacks_num = 1000  # We set this at the max of 1k for best speed
avg_ptcl_per_stack = len(ptcl_list) // stacks_num
rem_ptcl_lst_stack = len(ptcl_list) % stacks_num

with open(args.script, 'w') as script:
    script.write('#!/bin/bash\n')
    for stack_idx in range(stacks_num):
        stack_name = '{}{:04d}'.format(args.particle_prefix, stack_idx + 1)
        ptcl_sublist_start = stack_idx * avg_ptcl_per_stack + min(
            stack_idx, rem_ptcl_lst_stack)
        ptcl_sublist_end = (stack_idx + 1) * avg_ptcl_per_stack + min(
            stack_idx + 1, rem_ptcl_lst_stack)
        ptcl_sublist = ptcl_list[ptcl_sublist_start:ptcl_sublist_end]
        trfm_sublist = trfm_list[ptcl_sublist_start:ptcl_sublist_end]

        script.write('i3concat -dim 3 {} {}\n'.format(
            ' '.join(ptcl_sublist), stack_name + '.img'))
        with open(stack_name + '.trf', 'w') as trf, \
             open(stack_name + '.track', 'w') as ref:
            for trf_idx in range(len(trfm_sublist)):
                contents = None
                with open(trfm_sublist[trf_idx]) as ptcl_trf:
                    contents = ptcl_trf.read().split()
                init_subset = contents[0]
                new_subset = stack_name
                center_x = test_center_x
                center_y = test_center_y
                center_z = test_center_z + (trf_idx * test_boxsize_z)
                trf.write('{} {:d} {:d} {:d} {}\n'.format(
                    new_subset, center_x, center_y, center_z,
                    ' '.join(contents[4:])))
                ref.write('{} --> {} {:d} {:d} {:d} {}\n'.format(
                    init_subset, new_subset, center_x, center_y, center_z,
                    ' '.join(contents[4:])))
