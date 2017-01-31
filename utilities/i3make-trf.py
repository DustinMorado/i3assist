#!/usr/bin/env python
#. -*- coding: utf-8 -*-
"""Take a pos file with two point picked and calculate the orientation.

Written By Dustin Reed Morado
Last updated 05.11.2016
"""
from __future__ import print_function, division
import argparse
import numpy
import i3assist

parser = argparse.ArgumentParser(description=("Calculate orientation from two "
                                              "selected points in a pos file"),
                                 epilog="Writted by Dustin Morado 24.01.2017")

parser.add_argument("input_file", help="Input pos file", metavar="IN.pos")
parser.add_argument("subset", help="Field to use for the subset entry",
                    metavar="IN")
parser.add_argument("--reverse", "-r", action="store_true",
                    help="Reverse orientation from point two to point one.")
parser.add_argument("--position", "-p", type=float, default=0.5,
                    help=("Value from 0 to 1 specifying fraction of distance "
                          "from start point (Default=0.5 i.e. midpoint)."))

args = parser.parse_args()

if args.position < 0.0:
    print("WARNING: --position specified is less than 0.0 defaulting to 0.5")
    args.position = 0.5
elif args.position > 1.0:
    print("WARNING: --position specified is more than 1.0 defaulting to 0.5")
    args.position = 0.5
else:
    pass

with open(args.input_file) as position_file:
    is_first_position = True

    for line in position_file:
        if is_first_position is True:
            position_1 = numpy.array(
                [float(x) for x in line.split()], numpy.float_)
            is_first_position = False
        else:
            position_2 = numpy.array(
                [float(x) for x in line.split()], numpy.float_)

            if args.reverse:
                position_difference = position_1 - position_2
            else:
                position_difference = position_2 - position_1

            center = ((position_1 * (1.0 - args.position))
                      + (position_2 * args.position))

            xy_projection = numpy.sqrt(
                position_difference[0] ** 2 + position_difference[1] ** 2)
            phi = ((numpy.pi / 2.0)
                   + numpy.arctan2(position_difference[1],
                                   position_difference[0]))
            theta = ((numpy.pi / 2.0)
                     - numpy.arctan2(position_difference[2], xy_projection))
            trf = "{:s} ".format(args.subset)
            trf += "{:d} {:d} {:d} ".format(*center.astype(numpy.int_))
            trf += "{: 10.8f} {: 10.8f} {: 10.8f} ".format(0., 0., 0.)
            euler = i3assist.Euler(phi=phi, theta=theta, psi=0.0, unit="rad")
            matrix = euler.to_matrix()
            trf += matrix.trf_string()
            print(trf)
