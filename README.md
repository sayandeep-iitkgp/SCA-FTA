# SCA-FTA
This is the proof-of-concept code for the SCA-FTA attack. 
This version simulates the faults on different masked and 
unmasked S-Boxes (PRESENT and chi3, to be precise) with differnt
types of error-detection/correction (particularly, the vulnerable ones, 
such as detection on unshared value, detection/majority-voting correction 
on shared values).

The code is written in SAGEMATH and tested for the version 7.4. 

Run:
./SCA_FTA_simulation.sage
