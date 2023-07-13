# Ganz-Baker-Image-Velocimetry-Analysis

A set of Matlab functions for analyzing image velocimetry data on biological motility

---

## Setting up

### 0.) MATLAB

The functions make use of MATLAB's 

- Image Processing Toolbox
- Signal Processing Toolbox

There are three functions called from the Signal Processing toolbox (designfilt, filtfilt called by gutFreqWaveSpeedFinder.m and obtainMotilityParameters.m; xcorr.m called by obtainMotilityParameters.m). Again, if you don't have this, you can probably find something equivalent.

- intersections.m

In the original version of this program, there was one function (polyxpoly.m) required from the mapping toolbox; it called various sub-functions. We can now (Sept. 3, 2018) avoid this, downloading "intersections.m" (by Douglas Schwarz) from [here](https://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections);  this should be put in the MATLAB path. 


### 1.) Download PIVLab

* Download PIVLab from [here](http://www.mathworks.com/matlabcentral/fileexchange/27659-pivlab-time-resolved-particle-image-velocimetry--piv--tool). Remember they have a [website](http://pivlab.blogspot.com) if you'd like documentation on anything.
* Unzip the package into whichever directory you prefer.
* Add the PIVLab folders and subfolders to the [Matlab path](https://www.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html) ("Add with Subfolders" button, then choose the PIVLab main folder).

### 2.) Download this Image Velocimetry Analysis Suite

* Git Clone/Pull requests are recommended, though users may also download a .zip file and unzip it into whichever directory they prefer.
* Add this directory to the Matlab path, similar to the previous step.

### 3.) Have experimental motility images (or download our [sample image set](https://www.dropbox.com/s/zendq38z74jsxh8/ImageVelocimetryAnalysis_TestData.tif?dl=0))

* Experiments should have a directory structure as follows (folder names may be whatever): Experimental Directory -> Specimen folders -> Individual video folders -> .tif video. An example is [here](https://www.dropbox.com/s/uiveank3rndqsr7/GanzBakerIVADirectoryStructure.png?dl=0) (the main experimental directory is the folder with the date).

* The tiff videos can be a stack of tiffs, tiff images in numbered sequence, or even a combination of the two.

---

## Run Analysis

This section will use the sample image set as an example. Feel free to follow along with your own experimental motility image data.

### 1.) Make sure directory structure is correct.

* Begin by creating the right directory structure for your data. If using the sample motility data, prepare it like [step 1 here](https://www.dropbox.com/s/4wrd6fbymf4hs91/Starting%20Program.png?dl=0).

### 2.) Initialize a directory 

* Create a new directory for storing your analysis ([step 2](https://www.dropbox.com/s/4wrd6fbymf4hs91/Starting%20Program.png?dl=0))(or use the same directory as your experiment).

### 3.) Run program

* Open Matlab. Type "analyzeMotility;" (no quotes) to begin the program ([step 3](https://www.dropbox.com/s/4wrd6fbymf4hs91/Starting%20Program.png?dl=0)). If you want, enter the experiment and data directories as string arguments to avoid the next few prompts.

* If you did not add any arguments, it will now prompt you for the experimental directory, then it will prompt you for the directory to contain your analysis.

### 4.) Click yes 

* If this is the first time initializing your analysis directory, you will be prompted to create a new analysis file ([step 4](https://www.dropbox.com/s/4wrd6fbymf4hs91/Starting%20Program.png?dl=0)).

### 5.) Click "[Analyze Selection](https://www.dropbox.com/s/v1yhk16up409sp7/Program.png?dl=0)" (red text)

* Set the variables in the upper left "Variable Panel." These cannot be changed once the analysis has started.

* Click the "Analyze Selection" button to begin analysis on everything with a checkmark ([red text](https://www.dropbox.com/s/v1yhk16up409sp7/Program.png?dl=0)).

* About the checkboxes: the color represents what data you have performed analysis on. The checkmark indicates what analysis you wish to do. 

* Analysis is completed for an entire column before it moves to the next column (so all PIV is completed before the program asks you to outline any masks). 

* The last checkbox column indicates which data is included with the output when "Collect Analysis" is clicked (green text).

* Specimen directory names in the "Image Processing Control Panel" are highlighted orange, video directory names are highlighted cyan.

* The buttons in the "Image Processing Control Panel" next to the directory name will be available after analysis has been completed for that directory. See the [image](https://www.dropbox.com/s/v1yhk16up409sp7/Program.png?dl=0) for a brief description of their purpose.

### 6.) PIV will run.

* [PIV analysis](https://www.dropbox.com/s/o6gb7k0ggd184a8/PIV%20Running.png?dl=0) takes a long time. For a moderately powerful computer, data sets on the order of a few hundred GB take 12-24 hours to complete. Other computers have taken longer (several days).

### 7.) Outline the mask

* For our analysis, we only care about velocimetry data of zebrafish guts, so we mask it an discard data outside of the mask.

* To mask: The program will show an [image](https://www.dropbox.com/s/3a0l0tpndhk2bii/Begin%20Outlining.png?dl=0) of the current data set. Click on it to start drawing a polygon. Each click adds a new vertex.

* [Clicking](https://www.dropbox.com/s/u8jy0gb5cxz75vj/outline.png?dl=0) on the image continues drawing a polygon. Double click to enclose the area.

* You will then draw a [center line](https://www.dropbox.com/s/wj7yxwrt9egx8og/middleline.png?dl=0). The curvature of this line is important: new "local" coordinates will be generated from a smoothed version of the line you drew. That is, the curve you draw will define an anterior-posterior (AP) direction that changes as you move along the curve (the dorsal-ventral direction at any point along the curve is the perpendicular to whatever direction the AP is currently pointing). Double click to end the drawing. Click "Yes" if you like your drawing, "No" if you'd like to start over (just with the current image).

### 8.) Define the first correlation peak

* A set of images will pop up. The [top one](https://www.dropbox.com/s/qiarsnctl660yog/Begin%20XCorr.png?dl=0) (the cross-correlation curve, with a window named "Cross-correlation") will show two choices: the cross-correlation, and the cross-correlation with the median at each time-point subtracted (which removes position-independent motion). In the dialog box, pick which curve to use to find peaks. Then, you'll get a window with only this cross-correlation visualization. [Draw a line](https://www.dropbox.com/s/rko8sn36kfvncdf/Define%20line.png?dl=0), clicking two points along the first non-zero maximum. It doesn't need to be very accurate (the program will find the real maxima around the line you draw); the line you draw is simply helping the program distinguish **which** correlation peak out of the many you want to analyze.

* After double clicking to end the drawing, a [prompt](https://www.dropbox.com/s/rko8sn36kfvncdf/Define%20line.png?dl=0) will ask three questions. You can normally skip the first (it is asking how many seconds plus or minus it should look in y). The second and third are asking for the start and end of the correlation maximum, respectively. As shown in the figure, I tend to choose x=1 for the starting x and an x for which there is still a very well defined peak in the correlation plot (in this case, that's x=15). 

* Note that these values are not directly used in the analysis, so don't fret over them too much. They are simply used to figure out things like wave-speed (the inverse slope of the maxima) and the frequency of motility (the y-intercept).

* After you click ok, a [prompt](https://www.dropbox.com/s/s6v276ldoh3460e/Is%20XCorrGood.png?dl=0) will come up and show you the values the it has found. This is just a sanity check. For instance, the wave speed it found was ~25 seconds, which corresponds to the y-intercept of the correlation maxima, which looks good. Further, things like the R-squared value look good, so I feel confident that the analysis worked correctly. Click "Yes."

* **If you have an ugly plot**, this is likely still ok. It just means there was no noticeable repeating motility. If this is the case, draw a random line, click ok for the bounds, and the prompt will come up just like before. Select "No" this time. The program will ask if you want to fill with NaN's, click yes. It will then ask you to enter a frequency: just enter something around what the other motility data sets have been. It will then attempt to find an amplitude for this data set in a similar frequency range to what you entered.

### 9.) Once finished, collect your data

* In the [program screen](https://www.dropbox.com/s/v1yhk16up409sp7/Program.png?dl=0), click the "Collect Analysis" button (green text) to collect your data (you may now also revisit the buttons that were initially unavailable and do things like creating a movie with PIV vectors overlaying it).

* This button does two things: It sends your data to the [Matlab workspace](https://www.dropbox.com/s/cgj74ary6tspr9f/Variables%20sent%20to%20workspace.png?dl=0) (if you prefer to work in Matlab) and it also saves a [.csv file in the main analysis directory](https://www.dropbox.com/s/m011h0copbfbolw/Variables%20set%20to%20csv.png?dl=0) if you prefer working with any other software (such as excel).

* Output units: Frequencies are in units of 1/minutes. Fourier Transform Peaks, i.e. motility Amplitudes, are in units of microns. WaveSpeed slope: units of microns/second. (Assuming scales are correctly input into the GUI.) 

### 10.) Known bugs / problems

* There may be an error when using the 2019 version of PIVLab. This is probably easy to fix, but we haven't looked into it yet. (RP, April 9, 2019.)