#!/bin/bash
#To run an interactive shell session on a compute node use the -Is flag:
bsub -Is -q interactive bash

#To run an interactive session that can open windows on your X11 server, add the -XF flag:

# bsub -Is -XF -q interactive  bash
