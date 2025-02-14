README.txt

Last update: 2025/2/14

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data analysis and visualization of time-lapse recordings were performed using Fiji software.
Fiji is an image processing package—a “batteries-included” distribution of ImageJ2, bundling a lot of plugins which facilitate scientific image analysis. The package can be downloaded and installed from the following webpage, which also has information about installation on the various OSs: Windows, Linux, MacOS. 
http://imagej.net/Fiji
Fiji is released as open source under the GNU General Public License.
Plugins and other components have their own licenses.

Publication
Schindelin, J., Arganda-Carreras, I., Frise, E., Kaynig, V., Longair, M., Pietzsch, T., … Cardona, A. (2012). Fiji: an open-source platform for biological-image analysis. Nature Methods, 9(7), 676–682. doi:10.1038/nmeth.2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

For our dataset: Selection of individual subregions in the imaging field.
Due to the diffuse expression of iAChSnFR sensor (on somata and neuropil), the field of view (FOV) was overlaid with equally sized grids (~50 μm x ~50 μm)for selection and trace extraction. Using this procedure, we divided the FOV into 15x15 grids and obtained a total of 225 individual subregions in each mouse. 

The fluorescence trace of individual subregions over time was calculated by averaging the corresponding pixel values within each specified region. Time-series traces of all processed movies were extracted.

Instructions for use:

The macro "Split Image Into ROIs" splits the imaging field into equally sized grids of nxn tiles, where n is number of grids chosen by the user.

1. For installation and use, copy and paste the macros into "Startup Macros" in Fiji.
Plugins >> Macros >> Startup Macros. Click save.
2. To apply the changes, close and reopen Fiji software.
3. Load your movie file and run the macro 
Plugins >> Macros >> Split Image Into ROIs
4. Enter your desired number of grids

To view the time-series data of all selected ROIs, run the macro "Plot All Z-Axis profiles"
Plugins >> Macros >> Plot All Z-Axis Profiles

The data can then be extracted for further processing and analysis.
