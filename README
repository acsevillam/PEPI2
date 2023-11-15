Geant4PEPI 
-------------

 Geant4 simulation of an Edge illumination setup with tungsten anode source, sample and object masks (optional) and pixel based CdTe sensor 650um thick.
 3 geometries can be implemented: 1. double-mask (default), single-mask, conventional (no masks)
 Detector can be operated as ideal photon counter, realistic energy response with 1 energy threshold (default) or 2 thresholds

## REFERENCES
 Feel free to use and modify the code. If you intend to publish, please use these references:
 
 - Brombal, L., et al. "X-ray differential phase-contrast imaging simulations with Geant4." Journal of Physics D: Applied Physics 55.4 (2021): 045102.
 - Brombal, L., et al. "A Geant4 tool for edge-illumination X-ray phase-contrast imaging." Journal of Instrumentation 17.01 (2022): C01043.
 
## Authors
 L. Brombal, N. Poles
 
 INFN - Division of Trieste
 University of Trieste, Department of Physics

## GEOMETRY DEFINITION

standard setup, top view sketch
```
╔═════════════════════════════════════════════════════════════════════════════╗
║                                                                             ║
║                                                       █░║                   ║
║                                                       █░║ <-Pixirad         ║
║                       █░   ╔═╗       	  	 mask-> █░║   sensor          ║
║                       █░   ╚═╝            substrate  	█░║  	     	      ║
║ beam                  █░   ╔═╗                        █░║  		      ║	
║ ======>         	█░   ╚═╝                        █░║  		      ║
║                 	█░   ╔═╗ <-sample             	█░║                   ║
║                mask->	█░   ╚═╝			█░║		      ║	
║	    substrate	 ^                	  	█░║                   ║
║                        sample                         █░║ 	              ║
║                        mask                            ^                    ║
║                                                        detector             ║
║							 mask		      ║
╚═════════════════════════════════════════════════════════════════════════════╝
```
The definition of masks' pitches, apertures, thickness and position is flexible
Sample is composed by trapezoidal prisms with different inclinations, i.e. different refraction
 
## REQUIREMENTS
 
Geant4-10.05 or higher.
Older versions of Geant4 might work but they have not been tested

Open GL Utility Toolkit (GLUT) is needed to take advantage of Open GL visualization
	    
## PHYSICS
 
This simulation includes a custom implementation of the X-ray refraction.

X-ray refraction process is included in the PepiXrayRefraction.cc file and it is assigned to 'Gamma' particles through in the PepiPhysicsXrayRefraction.cc file.

Tungsten anode X-ray tube spectra, generated with the SpekCalc software, are included in the 'spectra' folder and span voltages from 40 to 100 kV at 10 kV interval.

The refractive index of materials are sourced externally. The data for a limited list of light materials is included in the 'data' folder along with their respective energy range (in keV).
 				
## VISUALIZATION
 
  The Visualization Manager is set in the main().
  The initialisation of the drawing is done via the commands :
  /vis/... in the macro vis.mac. In interactive session:
  PreInit or Idle > /control/execute vis.mac
 	
  The default view is a lateral view of the setup.
 	
  The tracks are drawn at the end of event, and erased at the end of run.

## HOW TO START ?
  1 Install Geant4 and configure the environment (source pathto/geant4.sh)
  
  2 Create and move in the build folder 
	% mkdir build
	% cd build
  
  3 Execute
  	% cmake ..
  	% make -j <nothreads>
  
  4a Execute pepi in 'batch mode' with macro file
  	% ./pepi run.mac	

  4b Execute pepi in 'interactive mode' with visualization and type your commands
  	% ./pepi
  	....
	Idle> run/beamOn 100
	....
	Idle> exit
 
 The run.mac file provides a template for generating EI images at different stepping and dithering positions
 
  
## OUTPUT FILES
 
G4Pepi produce raw images in 32-bit little endian format with dimension of 512 X 402
- When the detector type is in '0COL' mode, 1 image is produced (512 X 402 X 1)
- When the detector type is in '1COL' mode, 2 images are produced (512 X 402 X 2), the first being ideal and the second with Threshold 1
- When the detector type is in '2COL' mode, 3 images are produced (512 X 402 X 3), the first being ideal, the second with Threshold 1, the third with Threshold 2

 

for futher details on installation and building of Geant4 application refer to the Book for Application Developers
https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/index.html
