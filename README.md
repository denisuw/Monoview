# Monoview
Matlab-based tool for interactive monoplotting from aerial images

The tool is a simple interface written in Matlab using its graphical user interface (GUI) module, GUIDE. The tool enables the user to click on one main image (Cam1) and to be immediately directed towards other photos where the clicked point is visible (Cam2). Furthermore, by using the monoplotting concept, the user will also be able to inquire the 3D coordinates of the selected point. Other features of the tool include a distance calculator between two points in one image, and a vectorizer which can be used to digitize polygonal features from an image.

The principal inputs for this program are:
1. Aerial photos in .jpeg format placed in the same folder,
2. Initial or calculated exterior orientation (EO) parameters in a .txt file with the Omega Phi Kappa format, directly exportable from Photoscan,
3. Internal orientation (IO) parameters obtained from camera calibration and stored in an .xml file with Photoscan calibration format,
4. A DEM (Digital Elevation Model) of the area in .tif format, alongside its world file .tfw, and
5. The size of the pixel of the camera’s CCD in millimeters.
