# Master's Thesis: Topology-aware communication node placement for power grid restoration

This repository is dedicated to the Master's thesis project under the Chair of Computer Networks & Communications at the University of Passau. This project is implemented by Lovesh Bishnoi under the supervision of Prof. Dr.-Ing. Hermann de Meer and Anna Volkova.

## How to read this repository?

- CIST.py contains the working of the CIST algorithm.
- PSA_CIST.py contains the working of the PSA-CIST algorithm.
- steiner.py file contains the Steiner point calculation function which is used by CIST and PSA-CIST algorithms. Don't run it directly.
- BS input files Directory/Folder contains 10 different Base Station (BS) Topologies with a varying number of BS islands ranging from 3 to 12.
- PS input files Directory/Folder contains 3 different Power Systems (PSs) Topologies (Real-Life, Distributed, and Worst) along with 1 Synthetic PS Topology (which is specially created for BS Topology of 5 islands).

## How to implement the algorithms?

Implementation is performed on Python 3.9 with PyCharm 2021.3.3 (Professional Edition).

- Download all of the python packages
- Run either CIST.py or PSA_CIST.py
- Load one BS input file from the directory using the "Load BS file" button
- Load one PS input file from the directory using the "Load PS file" button
- Press the "Connect" button to restore the connectivity.
- Press the "Reset Screen" button to clear all of the data on the Tkinter canvas window.
- Press the "Zoom/Pan" button to zoom on canvas window.
