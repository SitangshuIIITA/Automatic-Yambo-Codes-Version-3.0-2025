#!/usr/bin/env python3

# Copyright (C) 2017-2025 elphmod Developers
# This program is free software under the terms of the GNU GPLv3 or later.

import elphmod
import matplotlib.pyplot as plt
import numpy as np
import sys

comm = elphmod.MPI.comm

path = 'KGM'
q, x, corners = elphmod.bravais.path(path, ibrav=4, N=80, moveG=0.01)

if len(sys.argv) > 1 and sys.argv[1] == '--prepare-q':
    if comm.rank == 0:
        q /= 2 * np.pi
        weight = 1 / len(q)

        with open('q.dat', 'w') as filqf:
            filqf.write('%d crystal\n' % len(q))

            for q1, q2, q3 in q:
                filqf.write('%12.10f %12.10f %12.10f %12.10f\n'
                    % (q1, q2, q3, weight))

 
