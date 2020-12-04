

# legoMosaic
This program turns a photo into a lego mosaic. It can build top down mosaic and sides out mosaics, and can use dithering to increase the number of colors
# Usage
## RunLegoMosaic
## options
- -p,  --picture : Reference picture to create mosaic from
- -o, --output : The output image name
- -n, --numHorizontal : Number of tiles horizontally
- -t, --tileSize : Size of the lego tiles in pixels
- -g, --gamma : Gamma exponent
- -s, --sidesOut : Generate lego mosaic with sides faced outwards
- -d, --dither : Use dithering in the mosaic
- -c, --colors : File listing the lego colors to use

#### Dither
Setting the dithering flag utilizes Floyd-Steinberg dithering in order to approximate larger pallette. This creates a more complex mosaic but the output will usually look better, especially at large scales.
![Manual Polygon](http://nathanbain.com/wikiImages/RunLegoMosaic/dither.png)
#### SidesOut
By default, the program will generate a mosaic of lego blocks with the studs facing toward the viewer. By setting the sidesOut flag, it will instead create a mosaic of the sides of blocks facing the viewer, and will even rotate the blocks 90 degrees in parts of the mosaic to improve the mosaic quality. This mosaic will usually look better, especially small mosaics will sharp details.
![Manual Polygon](http://nathanbain.com/wikiImages/RunLegoMosaic/sidesOut.png)
#### Gamma
Using [gamma correction](https://en.wikipedia.org/wiki/Gamma_correction) will alter the color to make it brighter or darker.
![Manual Polygon](http://nathanbain.com/wikiImages/RunLegoMosaic/gamma.png)
#### Color
There are 68 possible colors to choose from, and three sets of colors already specified: black and white, original colors from 1973, and all 68 colors. By using dithering, more colors can be approximated from fewer colors.
![Manual Polygon](http://nathanbain.com/wikiImages/RunLegoMosaic/colors.png)
