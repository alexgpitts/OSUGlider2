# Observing Ocean Wave Conditions Capstone Project

This is a capstone project for Oregon State University. The project involves 2 implementations of various calculations performed on acceleration data from ocean waves. We have an implementation in Python for onshore calculations and an implementation in C that will be loaded onto an AUV glider module. 

The Python version of our implementation is designed to process all the raw data from the gliders SD card after recovery. 

The C version of our implementation is designed to run on the glider during deployment, process the raw data for just the current time window, and stage the results to be sent over a satellite link pre-equipped on the glider. 
