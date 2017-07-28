# Ganz-Baker-Image-Velocimetry-Analysis

A set of Matlab functions for analyzing image velocimetry data on biological motility

---

## Setting up

### 1.) Download PIVLab

* Download PIVLab from [here](http://www.mathworks.com/matlabcentral/fileexchange/27659-pivlab-time-resolved-particle-image-velocimetry--piv--tool). Remember they have a [website](http://pivlab.blogspot.com) if you'd like documentation on anything.
* Unzip the package into whichever directory you prefer.
* Add the PIVLab folders and subfolders to the [Matlab path](https://www.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html) ("Add with Subfolders" button, then choose the PIVLab main folder).

### 2.) Download this Image Velocimetry Analysis Suite

* Git Clone/Pull requests are recommended, though users may also download a .zip file and unzip it into whichever directory they prefer.
* Add this directory to the Matlab path, similar to the previous step.

### 3.) Have experimental motility images (or download our **sample image set**)

* Experiments should have a directory structure as follows (folder names may be whatever): Experimental Directory -> Specimen folders -> Individual video folders -> .tif video. An example is [here](https://www.dropbox.com/s/uiveank3rndqsr7/GanzBakerIVADirectoryStructure.png?dl=0) (the main experimental directory is the folder with the date).

* The tiff videos can be a stack of tiffs, tiff images in numbered sequence, or even a combination of the two.
