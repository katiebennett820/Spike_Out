'''

IF YOU ARE LOOKING AT A BINARY (OR TRIPLET) SYSTEM:
This code empirically estimates the shape of and removes an unwanted diffraction spike from a second (or third) spectrum on detector 
that crosses your desired spectrum. 
This was designed for HST WFC3/G280 (UVIS) data, but should be applicable to most datasets. 
Most of the script is automated, but you do need to take certain steps to make sure you are running things correctly.
There are also a lot of plots that you can uncomment to check your progress (see 'UNCOMMENT TO PLOT HERE' signs).
The '#sys.exit('Notebook stopped')' lines will stop the rest of the script running so you can check to make sure things are working.  

Overall steps: 

1. Model center of diffraction spike as a line going across frame.
2. Along each row, spike also has a horizontal shape. Determine the width of the diffraction spike in the spectral direction (7 pixels for UVIS).
3. Take the median of the entire diffraction spike shape along a given interval length in the spatial direction (I have the medial per 40 pixels here). 
4. In each interval, empirically model the diffraction spike as the median diffraction spike. 
5. Subtract the background flux from the median diffraction spike. 
6. Subtract the median diffraction spike from each exposure. 


Manual changes you need to input through notebook: 

1. Section 1: Set path to where fits files are kept - see note at beginning of Section 1 for details on how to structure directories
2. Section 1: In last line, set "kind" to whatever fits file type you are using (ima, flt, etc.)
3. Section 2: Enter size of x-axis and y-axis
4. Section 3: Set appropriate estimates for slope (m) and y-intercept (b) of line used to model spike. 
5. Section 3: Uncomment plot to plot image with estimated line - manually iterate this until line overlaps with center of spike. 
    **Be careful here - at some estimated slopes, when you round x and y to integers (which is necessary b/c you use these as your indices), 
    **the points overlap, which will mess up the analysis. Avoid this! 
6. Section 3: You will need to input med_range, the interval length (in the spatial direction) over which you want to calculate median. 
7. Section 3: You can also adjust x_range (the horizontal range around which you will search for spike's shape) and lim (the half width of the spike shape)



written by K. Bennett 1/17/24


'''

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import glob
import sys, os, shutil
import pathlib
from pathlib import Path 

#------SECTITON 1: Gather Data------------------------

#Note that data_path requires fits files to be one directory lower than last directory listed 
#This is set up this way because this is how the data are downloaded from the MAST website
#In this case, files are within folders stored in "HST" directory 
data_path='/Users/katiebennett/Documents/LTT1445Ab/UVIS_Zafar/09_27_20/MAST_2023-11-03T1603/HST' #path to data
save_path='/Users/katiebennett/Documents/LTT1445Ab/UVIS_Zafar/09_27_20_diff_corr/HST_flt_files_final/' #path to save corrected data to

#Pull all raw data files (in this case, .flt files from HST/WFC3 G280 data)
#This function is written courtesy of FIREFLy.py written by Rustamkulov and Sing 
def parse_files(path, kind):
    file_list = []
    paths = sorted(Path(path).iterdir(), key = os.path.getatime)
    paths = [f for f in paths if not str(f).startswith('.')]    
    for p in (paths):
        if '.DS' not in str(p):
            s = sorted(Path(p).iterdir(), key=os.path.getmtime)
            s = [str(q) for q in s]
            for f in s:
                if kind in f and '.DS' not in f:
                    file_list.append(f)    
    return sorted(file_list)

#If not using flt.fits files, change line below to appropriate file type
flts=parse_files(data_path, 'flt.fits') #Get list of all fits files 



#------Section 2: Model center of diffraction spike with a line------------

y_axis=590 #y axis size
x_axis=2250 #x axis size

#Rough region around diffraction spike
rough_lo_end=300
rough_hi_end=420

#Run script through all the files 
for file in flts: 
    data_info=fits.open(file)
    #Gathering pertinent information from header here 
    date=data_info[0].header['DATE-OBS'] 
    time=data_info[0].header['TIME-OBS']
    start_time=data_info[0].header['EXPSTART']
    rootname=data_info[0].header['ROOTNAME']

    img=data_info['SCI'].data #Get image data
    img=img.T #transpose image so it reads (x,y) 
    img=np.array(img) #Need to turn image into a numpy array in order to manipulate it 
    img_corr=np.array(img) #Save a copy, which will become "corrected image"

#------Section 3: Model diffraction spike---------------------

    #Line 1 (brighter of the two, to the right)
    #Need to manually adjust m and b - uncomment plot below to see if line overlaps with desired diffraction spike
    m= -1.06 #0.942
    b= 1161.99 #1203 #-485  
    x=np.arange(0, x_axis) #Set appropriate x range 
    y=m*x+b

    for i in range(len(y)): #Need x and y to be integers 
        y[i]=round(y[i])
 
    #Only need x-range that lies within detector, so use top and bottom indices here to cut off rest of line equation that does not fall on detector
    top=np.where(y==y_axis)
    bottom=np.where(y==0)
    x_top, x_bottom=x[top], x[bottom]
    x=np.ndarray.flatten(np.linspace(x_bottom, x_top, y_axis))
    y=m*x+b 
    
    for i in range(len(y)): #Again, round new x and y values to integers
        x[i]=round(x[i])
        y[i]=round(y[i])

    '''
    #UNCOMMENT TO PLOT HERE 
    #Plot this 
    plt.plot(x, y, 'o', ms=1, color='aqua')
    plt.imshow(img.T, norm='log', origin='lower', cmap = 'Greys_r')
    plt.show()
    
    sys.exit('Notebook stopped')
    '''
    
    #Now, we look at things row by row to see what the shape of the spike looks like. 
    #We will use the median flux per 5 pixels (in the spatial direction) to model and remove flux 
    med_range=40 #Can change over what pixel range you want flux to be calculated
    med_range_width = med_range-1 #value needed to calculated upper end of each interval range 
    for i in np.arange(0, y_axis, med_range): 
        #This first part is just getting the correct indices for the data image and y line we have created
        #Set indices here - each interval goes from (lo_index, hi_index) along shape of spike
        lo_index=i 
        hi_index=i+med_range_width
        x_range=30 #set x-range - this is the horizontal distance from the spike center (y) around which we search for flux from diffraction spike 
        lim=7 #Horizontal number of pixels around center of spike that are part of diffraction pattern
        
        #Get equivalent indices for y line (Find indices where y has values in desired range)
        index=np.where(np.logical_and(y>=lo_index, y<=hi_index)) 
        index=np.ndarray.flatten(np.asarray(index)) 
        up_lim=index[0] #Set lower index limits f
        lo_lim=index[-1] #Set upper index limits      
        
        '''
        #UNCOMMENT TO PLOT HERE

        #Save value of center of spike along given interval 
        val=[]
        for i in index:
            val.append(img[round(x[i]),round(y[i])])

        #This plots flux value along center of spike for given interval along spatial direction as defined above - 
        #Useful if you set med_range=y_axis to see how flux value changes along whole detector. Hopefully, it does not 
        #vary much except where it crosses the detector 
        plt.plot(y[up_lim:lo_lim+1], val)
        plt.xlabel('Y Index')
        plt.ylabel('Flux (e/s)')
        plt.show()

        sys.exit('Notebook stopped')
        '''

        #Figure out shape of spike along a single row 
        val_range_arr=np.zeros((x_range*2, len(index))) #Make empty array to save flux values to along each row within given interval 
        q=0 
        for i in index:
            #Create x min and max along each row (index) - will depend on size you have made your range
            #This is just for visualization, so it's okay if the x_range value is not exact 
            x_lo=int(x[i]-x_range) 
            x_hi=int(x[i]+x_range)

            '''
            #UNCOMMENT TO PLOT HERE
            #This shows the raw data overplotted with our region of focus 
            plt.imshow(img.T, norm='log', origin='lower', cmap = 'Greys_r')
            plt.plot(x, y, color='aqua')
            plt.plot(x+x_range, y, color='aqua')
            plt.plot(x-x_range, y, color='aqua')
            plt.ylim(0,590)
            plt.show()

            sys.exit('Notebook stopped')
            '''

            #Along each row (index), go through each column and save value of image at that pixel
            val_range=[]
            for j in range(x_lo, x_hi):
                val_range.append(img[j, i])
            val_range_arr[:,q]=val_range #Also save values in this array, which we will use to find the median
            q+=1
            
            #HERE IS WHERE THE MAGIC HAPPENS - calculation of median diffraction spike 
            val_range_med=np.median(val_range_arr, axis=1) #find median of flux along rows per given interval!   
            
            '''
            #UNCOMMENT TO PLOT HERE
            #Plot flux per row in given interval along with median flux - good visual tool to make sure everything is working 
            plt.plot(range(x_range*2), val_range, alpha=0.2, color='blue') #Plot flux per row               
        plt.plot(range(x_range*2), val_range_med, alpha=1, color='black') #Plot median flux per frame 
        plt.axvline(x=x_range+lim, ls='--', color='black')
        plt.axvline(x=x_range-lim, ls='--', color='black')
        #plt.ylim((0, 100))
        plt.xlabel('X Range (middle is center of spike)')
        plt.ylabel('Flux (e/s)')
        plt.show()
        '''

        #----------SECTION 4: Subtract Median Diffraction Spike----------------

        #It's time! Subtract median flux from all indices

        #Set x range affected by spike (based on 'lim', which you have set earlier in this notebook)
        spike_x_lo=x_range-lim 
        spike_x_hi=x_range+lim 

        #Subtract background flux from value of diffraction spike
        bgd_vals=np.concatenate((val_range_med[15:spike_x_lo], val_range_med[spike_x_hi:45]))
        val_range_med_bgd=np.median(bgd_vals) #background flux
        val_range_med=val_range_med-val_range_med_bgd

        spike_val=val_range_med[spike_x_lo:spike_x_hi] #only going to subtract flux in this range

        #Go through each row to subtract median flux:
        for i in index:
            #Get equivalent x,y integers and the x min and max per row
            x_center, y_center = int(x[i]), int(y[i]) 
            x_min = x_center - lim #-1
            x_max= x_center + lim  #-1
            y_arr=np.repeat(y_center, len(range(x_min,x_max))) #make array of y_center values 
            
            #HERE IS WHERE THE MAGIC HAPPENS - subtract median flux from image 
            if rough_lo_end<y_center<rough_hi_end:
                img_corr[x_min:x_max,y_center]=img[x_min:x_max,y_center]-spike_val #Subtract median spike from image 

            #This is portion we are interested in that is NOT getting subtracted (the region between x_range and lim)
            #Use this to estimate what we want flux to really be 
            x_lo=x_center-x_range #Lower limit of x_range
            x_hi=x_min #Lower junction between x_range and lim
            x_lo2=x_max #Upper junction between x_range and lim
            x_hi2=x_center+x_range #Upper limit of x_range 
            '''
            
            #Unfortunately, not that simple! If you have estimated the median flux well, the resulting corrected values
            #will plunge to zero, but we want them to reflect the flux in the surrounding pixels 

            #Calculating average values in x-range region that did not get subtracted. Do each side indivually, then average together
            avg1=np.mean(img_corr[x_lo:x_hi, y_center])
            avg2=np.mean(img_corr[x_lo2:x_hi2, y_center])
            avg=np.mean((avg1, avg2))

            #Calculate average scatter in x_range flux in order to figure out what corrected values to reject
            r_sq=[]
            for m in range(x_lo2, x_hi2):
                resid=np.abs(img_corr[m, y_center]-avg2) #Calculate residual on low side
                r_sq.append(resid)
            resid_avg=np.mean(np.array(r_sq)) #Calculate average residual 
            scale_factor=1 #Arbitrarily decided scale factor (you can mess with this )
            min_accepted_value=avg2-resid_avg*scale_factor #Set limit on "acceptable value" by subtracting average residual from mean 

            #Need to offset corrected flux  with line going through flux 
            offset=1 #This is just an indexing offset, to make sure the slope of our line is accurate 

            #Calculate a line the old fashioned way! Get (x1, y1), (x2, y2) then calculate m and b
            x1=x_min-offset
            y1=img_corr[x_min-offset,y_center]
            x2=x_max+offset
            y2=img_corr[x_max+offset,y_center]

            slope=(y2-y1)/(x2-x1)
            xline=np.arange(x_min-1, x_max+2)
            yline=slope*(xline-x1)+y1
            bline=y1-slope*x1
            yline_real=slope*xline+bline #Get into slope-intercept form 

            #Crop line a bit so it is the same length as the number of pixels with corrected values 
            xline=xline[1:-1] 
            yline_real=yline_real[1:-1]
            
            #If value of corrected flux drops below accepted value, apply this correction 
            ex=0
            for i in range(x_min, x_max):
                if img_corr[i, y_center]<min_accepted_value:
                    img_corr[i, y_center]=img_corr[i, y_center]+(yline_real[ex]-img_corr[i, y_center])
                ex+=1
            '''
            '''
            #UNCOMMENT TO PLOT HERE
            #Plot line, original flux, median flux, and corrected flux to see how you did! 
            #up and down are just to investigate the region right around where the spike crosses the spectrum. Adjust accordingly 
            up=350
            down=370
            if y_center>up and y_center<down:
                #This plots flux, median flux, and corrected flux
                plt.title(y_center)
                #plt.plot(x1, y1, 'x', color='fuchsia')
                #plt.plot(x2, y2, 'x', color='fuchsia')
                #plt.plot(xline, yline_real, color='fuchsia')
                plt.plot(range(x_min, x_max), spike_val, color='black', label='Median Diffraction Flux')
                plt.plot(range(x_lo, x_hi2), img[x_lo:x_hi2,y_center], color='red', label='Flux at ith row')
                plt.plot(range(x_lo, x_hi2), img_corr[x_lo:x_hi2,y_center], color='purple', label='Corrected flux at ith row')                    
                plt.plot()
                plt.legend()
                plt.show()
           
    
    #UNCOMMENT TO PLOT HERE
    #Plot original data and corrected data to see how you did! 
    fig, axs = plt.subplots(2,1)
    axs[0].imshow(img.T, norm='log', origin='lower', cmap = 'magma', vmin=0.1, vmax=8.5e4)
    axs[1].imshow(img_corr.T, norm='log', origin='lower', cmap = 'magma', vmin=0.1, vmax=8.5e4)
    plt.show()
    '''
    #sys.exit('stop')
    

    
#------SAME THING AGAIN FOR SECOND DIFFRACTION SPIKE, IF NEEDED

#------Section 3: Model diffraction spike---------------------

    #Line 1 (brighter of the two, to the right)
    #Need to manually adjust m and b - uncomment plot below to see if line overlaps with desired diffraction spike
    m= -1.06 #0.942
    b= 1125.99 #1153 #1203 #-485  
    x=np.arange(0, x_axis) #Set appropriate x range 
    y=m*x+b

    for i in range(len(y)): #Need x and y to be integers 
        y[i]=round(y[i])
 
    #Only need x-range that lies within detector, so use top and bottom indices here to cut off rest of line equation that does not fall on detector
    top=np.where(y==y_axis)
    bottom=np.where(y==0)
    x_top, x_bottom=x[top], x[bottom]
    x=np.ndarray.flatten(np.linspace(x_bottom, x_top, y_axis))
    y=m*x+b 
    
    for i in range(len(y)): #Again, round new x and y values to integers
        x[i]=round(x[i])
        y[i]=round(y[i])

    '''
    #UNCOMMENT TO PLOT HERE 
    #Plot this 
    plt.plot(x, y, 'o', ms=1, color='aqua')
    plt.imshow(img.T, norm='log', origin='lower', cmap = 'Greys_r')
    plt.show()
    
    #sys.exit('Notebook stopped')
    '''
    
    #Now, we look at things row by row to see what the shape of the spike looks like. 
    #We will use the median flux per 5 pixels (in the spatial direction) to model and remove flux 
    med_range=40 #Can change over what pixel range you want flux to be calculated
    med_range_width = med_range-1 #value needed to calculated upper end of each interval range 
    for i in np.arange(0, y_axis, med_range): 
        #This first part is just getting the correct indices for the data image and y line we have created
        #Set indices here - each interval goes from (lo_index, hi_index) along shape of spike
        lo_index=i 
        hi_index=i+med_range_width
        x_range=30 #set x-range - this is the horizontal distance from the spike center (y) around which we search for flux from diffraction spike 
        lim=7 #Horizontal number of pixels around center of spike that are part of diffraction pattern
        
        #Get equivalent indices for y line (Find indices where y has values in desired range)
        index=np.where(np.logical_and(y>=lo_index, y<=hi_index)) 
        index=np.ndarray.flatten(np.asarray(index)) 
        up_lim=index[0] #Set lower index limits f
        lo_lim=index[-1] #Set upper index limits      
        
        '''
        #UNCOMMENT TO PLOT HERE

        #Save value of center of spike along given interval 
        val=[]
        for i in index:
            val.append(img[round(x[i]),round(y[i])])

        #This plots flux value along center of spike for given interval along spatial direction as defined above - 
        #Useful if you set med_range=y_axis to see how flux value changes along whole detector. Hopefully, it does not 
        #vary much except where it crosses the spectra 
        plt.plot(y[up_lim:lo_lim+1], val)
        plt.xlabel('Y Index')
        plt.ylabel('Flux (e/s)')
        plt.show()

        #sys.exit('Notebook stopped')
        '''

        #Figure out shape of spike along a single row 
        val_range_arr=np.zeros((x_range*2, len(index))) #Make empty array to save flux values to along each row within given interval 
        q=0 
        for i in index:
            #Create x min and max along each row (index) - will depend on size you have made your range
            #This is just for visualization, so it's okay if the x_range value is not exact 
            x_lo=int(x[i]-x_range) 
            x_hi=int(x[i]+x_range)

            '''
            #UNCOMMENT TO PLOT HERE
            #This shows the raw data overplotted with our region of focus 
            plt.imshow(img.T, norm='log', origin='lower', cmap = 'Greys_r')
            plt.plot(x, y, color='aqua')
            plt.plot(x+x_range, y, color='aqua')
            plt.plot(x-x_range, y, color='aqua')
            plt.ylim(0,590)
            plt.show()

            #sys.exit('Notebook stopped')
            '''

            #Along each row (index), go through each column and save value of image at that pixel
            val_range=[]
            for j in range(x_lo, x_hi):
                val_range.append(img[j, i])
            val_range_arr[:,q]=val_range #Also save values in this array, which we will use to find the median
            q+=1
            
            #HERE IS WHERE THE MAGIC HAPPENS - calculation of median diffraction spike 
            val_range_med=np.median(val_range_arr, axis=1) #find median of flux along rows per given interval!   
            
            '''
            #UNCOMMENT TO PLOT HERE
            #Plot flux per row in given interval along with median flux - good visual tool to make sure everything is working 
            plt.plot(range(x_range*2), val_range, alpha=0.2, color='blue') #Plot flux per row               
        plt.plot(range(x_range*2), val_range_med, alpha=1, color='black') #Plot median flux per frame 
        plt.axvline(x=x_range+lim, ls='--', color='black')
        plt.axvline(x=x_range-lim, ls='--', color='black')
        #plt.ylim((0, 100))
        plt.xlabel('X Range (middle is center of spike)')
        plt.ylabel('Flux (e/s)')
        plt.show()
        '''
        

        #----------SECTION 4: Subtract Median Diffraction Spike----------------

        #It's time! Subtract median flux from all indices

        #Set x range affected by spike (based on 'lim', which you have set earlier in this notebook)
        spike_x_lo=x_range-lim 
        spike_x_hi=x_range+lim 

        #Subtract background flux from value of diffraction spike
        bgd_vals=np.concatenate((val_range_med[15:spike_x_lo], val_range_med[spike_x_hi:45]))
        val_range_med_bgd=np.median(bgd_vals) #background flux
        val_range_med=val_range_med-val_range_med_bgd

        spike_val=val_range_med[spike_x_lo:spike_x_hi] #only going to subtract flux in this range

        #Go through each row to subtract median flux:
        for i in index:
            #Get equivalent x,y integers and the x min and max per row
            x_center, y_center = int(x[i]), int(y[i]) 
            x_min = x_center - lim #-1
            x_max= x_center + lim #-1
            y_arr=np.repeat(y_center, len(range(x_min,x_max))) #make array of y_center values 
            
            #HERE IS WHERE THE MAGIC HAPPENS - subtract median flux from image 
            #if np.median(img[x_min:x_max,y_center])>np.median(spike_val):
            if rough_lo_end<y_center<rough_hi_end:
                img_corr[x_min:x_max,y_center]=img[x_min:x_max,y_center]-spike_val #Subtract median spike from image 
            #else: 
                #img_corr[x_min:x_max,y_center]=img[x_min:x_max,y_center]

            #This is portion we are interested in that is NOT getting subtracted (the region between x_range and lim)
            #Use this to estimate what we want flux to really be 
            x_lo=x_center-x_range #Lower limit of x_range
            x_hi=x_min #Lower junction between x_range and lim
            x_lo2=x_max #Upper junction between x_range and lim
            x_hi2=x_center+x_range #Upper limit of x_range 
            '''
            #Unfortunately, not that simple! If you have estimated the median flux well, the resulting corrected values
            #will plunge to zero, but we want them to reflect the flux in the surrounding pixels 

            #Calculating average values in x-range region that did not get subtracted. Do each side indivually, then average together
            avg1=np.mean(img_corr[x_lo:x_hi, y_center])
            avg2=np.mean(img_corr[x_lo2:x_hi2, y_center])
            avg=np.mean((avg1, avg2))

            #Calculate average scatter in x_range flux in order to figure out what corrected values to reject
            r_sq=[]
            for m in range(x_lo2, x_hi2):
                resid=np.abs(img_corr[m, y_center]-avg2) #Calculate residual on low side
                r_sq.append(resid)
            resid_avg=np.mean(np.array(r_sq)) #Calculate average residual 
            scale_factor=1 #Arbitrarily decided scale factor (you can mess with this )
            min_accepted_value=avg2-resid_avg*scale_factor #Set limit on "acceptable value" by subtracting average residual from mean 

            #Need to offset corrected flux  with line going through flux 
            offset=1 #This is just an indexing offset, to make sure the slope of our line is accurate 

            #Calculate a line the old fashioned way! Get (x1, y1), (x2, y2) then calculate m and b
            x1=x_min-offset
            y1=img_corr[x_min-offset,y_center]
            x2=x_max+offset
            y2=img_corr[x_max+offset,y_center]

            slope=(y2-y1)/(x2-x1)
            xline=np.arange(x_min-1, x_max+2)
            yline=slope*(xline-x1)+y1
            bline=y1-slope*x1
            yline_real=slope*xline+bline #Get into slope-intercept form 

            #Crop line a bit so it is the same length as the number of pixels with corrected values 
            xline=xline[1:-1] 
            yline_real=yline_real[1:-1]
            
            #If value of corrected flux drops below accepted value, apply this correction 
            ex=0
            for i in range(x_min, x_max):
                if img_corr[i, y_center]<min_accepted_value:
                    img_corr[i, y_center]=img_corr[i, y_center]+(yline_real[ex]-img_corr[i, y_center])
                ex+=1
            '''
            '''
            #UNCOMMENT TO PLOT HERE
            #Plot line, original flux, median flux, and corrected flux to see how you did! 
            #up and down are just to investigate the region right around where the spike crosses the spectrum. Adjust accordingly 
            up=350
            down=370
            if y_center>up and y_center<down:
                    #This plots flux, median flux, and corrected flux
                plt.title(y_center)
                #plt.plot(x1, y1, 'x', color='fuchsia')
                #plt.plot(x2, y2, 'x', color='fuchsia')
                #plt.plot(xline, yline_real, color='fuchsia')
                plt.plot(range(x_min, x_max), spike_val, color='black', label='Median Diffraction Flux')
                plt.plot(range(x_lo, x_hi2), img[x_lo:x_hi2,y_center], color='red', label='Flux at ith row')
                plt.plot(range(x_lo, x_hi2), img_corr[x_lo:x_hi2,y_center], color='purple', label='Corrected flux at ith row')                    
                plt.plot()
                plt.legend()
                plt.show()
            
    
    #UNCOMMENT TO PLOT HERE
    #Plot original data and corrected data to see how you did! 
    fig, axs = plt.subplots(2,1)
    axs[0].imshow(img.T, norm='log', origin='lower', cmap = 'magma', vmin=0.1, vmax=8.5e4)
    axs[1].imshow(img_corr.T, norm='log', origin='lower', cmap = 'magma', vmin=0.1, vmax=8.5e4)
    plt.show()
    '''
    #sys.exit('stop')
    

    #Save new data into fits file
    data_info['SCI'].data=img_corr.T
    
    #Save into new directory 
    os.mkdir(save_path+rootname+'/')
    new_path=save_path+rootname+'/'
    data_info.writeto(new_path+date+'_'+rootname+'_flt.fits', overwrite=True)
