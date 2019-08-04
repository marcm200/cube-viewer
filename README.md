# cube-viewer
Simple viewer for 3D cubes built from stacks of bitmaps.

## Organization of README

1. Quick start
2. Background
3. Command line options
4. Limitations
5. Contact


### 1. Quick start

Compile the code using a C++ compiler to create the executable (assuming it is named cube-viewer.exe).
Run it in a command prompt or by using the DOS batch files `_EXAMPLE_construct.bat` and `_EXAMPLE_view.bat` 
 
`cube-viewer.exe cmd=view file=_filter400.ccb`

This creates a number of files named `_filter400_*.bmp` that contain what an observer sees from different positions in space.

`cube-viewer.exe cmd=construct file=bc1.ccb`

This constructs a cube file named `bc1.ccb` from images `00000NNN.bmp` in the current directory (see binary release for a set of images).

To view this cube use:

`cube-viewer.exe cmd=view file=bc1.ccb observer=1`

and to look into the object be defining some RGB colors as transparent (black is always considered transparent)

`cube-viewer.exe cmd=view file=bc1.ccb transp=0,0,0,100,100,100´`


### 2. Background

Constructing some Lyapunov images in 2D I wanted to look at some parameter changes as a 3D object. The code here is designed to do exactly that: Take images created elsewhere, put them into a 3D cube and look at it from various positions, defining some RGB interval as being transparent so as to maybe look into an object. The observer can be arbitrarily placed into space or predefined 
standard positions can be used (see command line options). Shading is accomplished by dimming the brightness of a resulting pixel by distance to the observer and by angle of a normal vector from a 2x2 pixel grid to the observer's view axis which is also the light source.


### 3. Command line options

The software can be used with two goals: construct a cube from a stack of bitmap images or to view an existing cube from various observer positions.

#### 3.1 Constructing a cube

To construct a cube one needs a stack of 24-bit RGB bitmaps of quadratic size. Those ought to be numbered from 1 to N where N=width of the image using names with leading zeros to 8 digits, e.g. 00000001.bmp to 00000400.bmp for a set of 400x400 pixel images.

The cube is then constructed via (image files must be in the same directory as the executable):

`cube-viewer.exe cmd=construct file=YOURFILENAME.ccb`

It saves the cube into a file called YOURFILENAME.ccb which stores the RGB data of the combined images.

#### 3.2 Viewing a cube

An existing cube can be viewed by using `cmd=view` as a command line value and takes several parameters.

E.g. `cube-viewer.exe cmd=view file=_jupiter400.ccb observer=-1 transp=0,0,0,10,11,12`

`file`being the filename of the cube to be read.

`óbserver=-1` being the observer's position. The -1 denotes that all predefined observer positions should be used one after the other. Predefined are positions above every edge of the cube and looking into the cube's center point as well as looking along every axis to a side plane's middle point and take values from 1 to 14.

It is also possible to explicitly state an observer position by using `óbserver=x,y,z` with three coordinates (whole numbers). The cube always occupies coordinates 0,0,0 (included) to position +width,+width,+width (excluded).

`transp=0,0,0,10,11,12` denotes the RGB range which is considered transparent. Given are two RGB values, here all points whose color lies between 0 <= R <=10 and 0 <= G <= 11 and 0 <= B <= 12 are considered to be transparent. Values are given in lowR,lowG,lowB,highR,highG,highB format. 

<b>Note:</b> Black (0,0,0) is always considered transparent, no matter what value is given here.

The resulting images are saved under the name of the cube with added observer position  coordinates.


### 4. Limitations

- The software comes with no warranty.
- It is tested on cubes where all sides are equally long and at most 1024 pixels.
- It lacks perspective in the sense that an object will not shrink in size for points farther away from the observer.
- The observer always looks to the center of the cube.
- The distance of the observer to the cube does not affect the apparent size of what is seen in the resulting image, the observer is always jumped to just still outside the cube.
- Axis are not checked to be in the three finger orientation.
- The light source is always the observer's initial position (head lamp).
- The cube always occupies the coordinates 0,0,0 (included) to position +width,+width,+width (excluded).


### 5. Contact

Please direct comments to: marcm200@freenet.de

Marc Meidlinger, August 2019
