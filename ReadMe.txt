
RJMCMC README

=============

INTRODUCTION

RJMCMC is the program for dynamic clustering of time-course data in which the cluster structure itself can be changed along time. Clusters can be born from the others or dead during the time course.

BUILD

To build the code, run build.sh under the root of the source. Currently only linux is supported (64bit recommended) and the output will be rjmcmc under the root directory.

We are planning to include a makefile for flexible build/install purposes.

DATA FORMAT

Please refer to Data_sample.txt under root for proper data format and usage of the program.

REFERENCE

This is the source code for time-variant clustering published in the following paper:

Wei Huang, Xiaoyi Cao, Fernando H. Biase, Pengfei Yu, Sheng Zhong. Time-variant clustering model for understanding cell fate decisions. PNAS, 111(44):E4797-E4806.

