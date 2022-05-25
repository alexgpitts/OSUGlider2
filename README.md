# Observing Ocean Wave Conditions Capstone Project (CS)

## Overview: 

The Glider Project is an ongoing research project conducted by Kipp Shearman, Pat Welch, and Kai-Fu Chang. The paramount goal of the glider project is to autonomously measure wave patterns and conditions via Autonomous Underwater Vehicle (AUV) gliders, then use this data to build polar plot spectrums for predictive and historical models. Data is already being collected via buoys, but there is no good way to measure wave data in the deep ocean. By incorporating the same technology on these buoys to ocean gliders, we can increase the number of data collection points in the ocean. This data can then be used by many industries. For example, commercial shipping, ship designers, and coastal engineers. 

### Current Conditions  

The state of wave measurements is mainly reliant on fixed moorings. Fixed moorings have the advantage of accessibility because they are close to shore. However, they suffer from several limitations that the Glider Project hopes to address. Namely, fixed moorings are expensive to install, prone to damage, and limited in their radius of measurement. 

### Implementation:  

The glider is already deployed for various other projects. Every few hours the glider surfaces to send its findings back to shore via satellite. The ECE team we are collaborating with implemented a module to be installed on these ocean gliders that will extract acceleration data while the glider floats on the surface. The raw acceleration data is recorded with the accelerometer implemented by the ECE team and saved to an SD card using their custom PCB. 

The project involves 2 implementations of various calculations performed on acceleration data from ocean waves. We have an implementation in Python for onshore calculations and an implementation in C that will be loaded onto an AUV glider module.   

The values we will calculate in our implementation include significant wave height, average wave period, peak wave period, mean zero up cross wave period, peak wave direction, peak wave power spectral density, and more. These calculations will be saved into pre-allocated memory on the glider at the end of each data collection period (about 20 minutes), queuing them to be sent over a satellite link the next time the glider surfaces. In addition, the raw acceleration data in the x, y, and z directions will be written to an SD card to be extracted later for further calculations. Additionally, we created an onshore program in python for more extensive calculations. Using the raw acceleration data from the SD card on the glider, this python program will calculate all the same information as the glider code, using various calculation methods, including the Welch method and banding using different windowing techniques. 
