#include "LegoMosaic.h"
#include <tclap/CmdLine.h>

using namespace TCLAP;

int main( int argc, char **argv )
{
  try
  {
    CmdLine cmd("Creates an image mosaic.", ' ', "2.0");

    ValueArg<string> colorNameArg( "c", "colors", "Colors to use", false, "legoData/allColors.txt", "string", cmd);

    SwitchArg randomizeArg( "r", "randomize", "Randomize the placement of the blocks", cmd, false );

    SwitchArg ditherArg( "d", "dither", "Use dithering in the mosaic", cmd, false );

    ValueArg<int> styleArg( "s", "style", "0 = studs-out; 1 = sides-out; 2 = more blocks; 3 = halfscale", false, 0, "int", cmd );

    ValueArg<double> gammaArg( "g", "gamma", "Gamma exponent", false, 1.0, "double", cmd );

    ValueArg<int> dilateDistanceArg( "e", "expand", "Expand dark areas", false, 0, "int", cmd );

    ValueArg<int> tileSizeArg( "t", "tileSize", "Size of the lego tiles in pixels", false, 10, "int", cmd );

    ValueArg<int> numberArg( "n", "numHorizontal", "Number of tiles horizontally", true, 20, "int", cmd);

    ValueArg<string> outputArg( "o", "output", "Output name", true, "out", "string", cmd);

    ValueArg<string> pictureArg( "p", "picture", "Reference picture to create mosaic from", true, "in", "string", cmd);

    cmd.parse( argc, argv );

    bool randomize                    = randomizeArg.getValue();
    bool dither                       = ditherArg.getValue();
    int style                         = styleArg.getValue();
    string inputName                  = pictureArg.getValue();
    string outputName                 = outputArg.getValue();
    string colorName                  = colorNameArg.getValue();
    int numHorizontal                 = numberArg.getValue();
    int tileSize                      = tileSizeArg.getValue();
    int dilateDistance                = dilateDistanceArg.getValue();
    double gamma                      = gammaArg.getValue();

    if( VIPS_INIT( argv[0] ) ) return( -1 );

    RunLegoMosaic( inputName, outputName, numHorizontal, tileSize, style, dither, randomize, dilateDistance, gamma, colorName );
  }
  catch (ArgException &e)  // catch any exceptions
  {
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }
  return 0;
}