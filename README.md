# Observing Ocean Wave Conditions Capstone Project

## Overview: 

The Glider Project is an ongoing research project conducted by Kipp Shearman, Pat Welch, and Kai-Fu Chang. The paramount goal of the glider project is to autonomously measure wave patterns and conditions via Autonomous Underwater Vehicle (AUV) gliders, then use this data to build polar plot spectrums for predictive and historical models. Data is already being collected via buoys, but there is no good way to measure wave data in the deep ocean. By incorporating the same technology on these buoys to ocean gliders, we can increase the number of data collection points in the ocean. This data can then be used by many industries. For example, commercial shipping, ship designers, and coastal engineers. 

### Current Conditions  

The state of wave measurements is mainly reliant on fixed moorings. Fixed moorings have the advantage of accessibility because they are close to shore. However, they suffer from several limitations that the Glider Project hopes to address. Namely, fixed moorings are expensive to install, prone to damage, and limited in their radius of measurement. 

### Implementation:  

The project involves 2 implementations of various calculations performed on acceleration data from ocean waves. We have an implementation in Python for onshore calculations and an implementation in C that will be loaded onto an AUV glider module.  

The Python version of our implementation is designed to process all the raw data from the gliders SD card after recovery.  

The C version of our implementation is designed to run on the glider during deployment, process the raw data for just the current time window, and stage the results to be sent over a satellite link pre-equipped on the glider.  
