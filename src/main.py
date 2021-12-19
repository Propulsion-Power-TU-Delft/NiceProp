#! /usr/bin/env python
# -*- coding: utf-8 -*-

#############################################################################
# NiceProp - Interactively learning NICFD
# Authors: ir. A. Giuffre', Dr. ir. M. Pini
# Content: Program main
# 2021 - TU Delft - All rights reserved
#############################################################################


import time
import warnings
from GUI import *


if __name__ == '__main__':
    start_time = time.time()
    warnings.filterwarnings("ignore")

    root = tk.Tk()
    root.title('Welcome to NiceProp!')
    App(root).grid()
    root.mainloop()

    print("\n Elapsed time: %10.1f seconds" % (time.time() - start_time))
