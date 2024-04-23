# Spike_Out
This code removes erroneous diffraction spikes that crosses your spectrum, which can occur in binary or triplet systems. This was developed for HST WFC3/G280 UVIS data, but should be applicable to any dataset.

Most of the script is automated, but you do need to take certain steps to make sure you are running things correctly. Be sure to open the UVIS_diffraction_spike_correction_commented.py file to learn more and make the necessary adjustments.

Overall steps: 

1. Model center of diffraction spike as a line going across frame.
2. Along each row, spike also has a horizontal shape. Determine the width of the diffraction spike in the spectral direction (7 pixels for UVIS).
3. Take the median of the entire diffraction spike shape along a given interval length in the spatial direction (I have the medial per 40 pixels here). 
4. In each interval, empirically model the diffraction spike as the median diffraction spike. 
5. Subtract the background flux from the median diffraction spike. 
6. Subtract the median diffraction spike from each exposure. 
7. Repeat as needed with as many erroneous diffraction spikes as you have! 
