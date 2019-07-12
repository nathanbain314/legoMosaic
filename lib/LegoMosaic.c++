#include "LegoMosaic.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

double si = 1;


void fit(float &c, int type)
{
  if( type == 0 )
  {
    c = max(c,0.0f);
    c = min(c,100.0f);
  }
  else
  {
    c = max(c,-128.0f);
    c = min(c,127.0f);
  }
}

void rgbToLab( int r, int g, int b, float &l1, float &a1, float &b1 )
{
  vips_col_sRGB2scRGB_8( r, g, b, &l1,&a1,&b1 );
  vips_col_scRGB2XYZ( l1, a1, b1, &l1, &a1, &b1 );
  vips_col_XYZ2Lab( l1, a1, b1, &l1, &a1, &b1 );
}

void changeColorspace( Tree &tree, Point &center, unsigned char * c, float * c2, int start, int end, int width, double gamma, bool gammutMapping, ProgressBar * changingColorspace, bool labSpace, bool quiet )
{
  bool show = !quiet && !(start);

  for( int height = start, index = start*width; height < end; ++height )
  {
    for( int w = 0; w < width; ++w, ++index )
    {
      for( int k = 0; k < 3; ++k )
      {
        int n10 = c[3*index+k];
        
        float v = (float)n10;

        v /= 256.0;

        v = pow( v, gamma );
        c2[3*index+k] = 256 * v;
      }
      
      if( labSpace ) rgbToLab( c2[3*index+0], c2[3*index+1], c2[3*index+2], c2[3*index+0], c2[3*index+1], c2[3*index+2] );      

      if(gammutMapping)
      {
        Point point_query(c2[3*index+0], c2[3*index+1], c2[3*index+2]);

        Segment segment_query(center,point_query);

        if(tree.do_intersect(segment_query))
        {
          Point closest_point = tree.closest_point(point_query);

          for( int k = 0; k < 3; ++k )
          {
            c2[3*index+k] = closest_point[k];
          }
        }
      }
    }

    if( show ) changingColorspace->Increment();
  }
}

void generateImages( vector< vector< vector< unsigned char > > > &brickImages, vector< pair< int, int > > &dimensions )
{
  ifstream shapeData( "legoData/studsOut.dat", ios::binary );

  size_t len;

  int rotImage[14] = 
  {
    1,
    2,
    0,
    2,
    1,
    2,
    2,
    1,
    2,
    0,
    1,
    2,
    1,
    0,
  };

  ProgressBar *generatingImages = new ProgressBar( 14*68, "Loading images" );

  for( int i = 0, index = 0; i < 14; ++i )
  {
    for( int j = 0; j < 68; ++j )
    {
      shapeData.read( (char *)&len, sizeof(size_t) );

      void * data = (void *) malloc (sizeof(unsigned char)*len);

      shapeData.read( (char *)data, len*sizeof(char) );

      VipsBlob *blob = vips_blob_new( NULL, data, len );

      VImage image = VImage::pngload_buffer(blob).resize(si);

      unsigned char * c;

      if( rotImage[i] < 2 )
      {
        c = ( unsigned char * )image.data();
        dimensions[index] = pair< int, int >( image.width(), image.height() );
      }
      else
      {
        c = (unsigned char *)image.rot(VIPS_ANGLE_D90).data();
        dimensions[index] = pair< int, int >( image.height(), image.width() );
      }

      brickImages[index].push_back( vector< unsigned char >(c, c + 3*image.width()*image.height() ) );

      if( rotImage[i] == 1 )
      {
        dimensions[index+1] = pair< int, int >( image.height(), image.width() );

        c = (unsigned char *)image.rot(VIPS_ANGLE_D90).data();
        brickImages[index+1].push_back( vector< unsigned char >(c, c + 3*image.width()*image.height() ) );
      }
      else if( rotImage[i] == 2 )
      {
        dimensions[index+1] = pair< int, int >( image.width(), image.height() );

        c = ( unsigned char * )image.data();
        brickImages[index+1].push_back( vector< unsigned char >(c, c + 3*image.width()*image.height() ) );
      }

      generatingImages->Increment();
    }

    index += rotImage[i] > 0 ? 2 : 1;
  }

  generatingImages->Finish();
}

void buildTopLevel( string outputImage, int start, int end, int outputWidth, int outputHeight, vector< vector< int > > &brickData, vector< vector< vector< unsigned char > > > &brickImages, int maxTileWidth, int maxTileHeight, ProgressBar *topLevel )
{
  bool singleImage = (maxTileWidth > 256);

  for( int tileOffsetY = start; tileOffsetY < end; tileOffsetY += maxTileHeight )
  {
    int tileHeight = min( outputHeight - tileOffsetY, maxTileHeight );

    for( int tileOffsetX = 0; tileOffsetX < outputWidth; tileOffsetX += maxTileWidth )
    {
      if( start == 0 && !singleImage ) topLevel->Increment();

      int tileWidth = min( outputWidth - tileOffsetX, maxTileWidth );

      unsigned char * tileData = ( unsigned char * )calloc (tileHeight*tileWidth*3,sizeof(unsigned char));

      for( int k = 0; k < brickData.size(); ++k )
      {
        if( singleImage ) topLevel->Increment();

        int brick = brickData[k][0];
        int colorIndex = brickData[k][1];
        int bh = brickData[k][2];
        int bw = brickData[k][3];
        int bi = brickData[k][4];
        int bj = brickData[k][5];

        if( bj >= tileOffsetX + tileWidth || bj + bw < tileOffsetX || bi >= tileOffsetY + tileHeight || bi + bh < tileOffsetY ) continue;

        for( int y1 = 0; y1 < bh; ++y1 )
        {
          for( int x1 = 0; x1 < bw; ++x1 )
          {
            int x = x1 + bj - tileOffsetX;
            int y = y1 + bi - tileOffsetY;

            if( !(x >= 0 && x < tileWidth && y >= 0 && y < tileHeight) ) continue;

            int index1 = 3*(y*tileWidth+x);

            int index2 = 3*(y1*bw+x1);

            tileData[index1+0] = brickImages[brick][colorIndex][index2+0];
            tileData[index1+1] = brickImages[brick][colorIndex][index2+1];
            tileData[index1+2] = brickImages[brick][colorIndex][index2+2];
          }
        }
      }

      string tileOutputName = singleImage ? outputImage : string(outputImage).append(to_string(tileOffsetX/256)+"_"+to_string(tileOffsetY/256)+".jpeg");

      if( singleImage )
      {
        topLevel->Finish();

        cout << "Saving image" << endl;

        VImage::new_from_memory( tileData, tileHeight*tileWidth*3, tileWidth, tileHeight, 3, VIPS_FORMAT_UCHAR ).vipssave((char *)tileOutputName.c_str());
      }
      else
      {
        VImage::new_from_memory( tileData, tileHeight*tileWidth*3, tileWidth, tileHeight, 3, VIPS_FORMAT_UCHAR ).jpegsave((char *)tileOutputName.c_str(), VImage::option()->set( "optimize_coding", true )->set( "strip", true ) );
      }

      free( tileData );
    }
  }
}

void buildImage( string outputImage, int outputWidth, int outputHeight, vector< vector< int > > &brickData, vector< vector< vector< unsigned char > > > &brickImages )
{
  // Save image as static image or zoomable image
  if( vips_foreign_find_save( outputImage.c_str() ) != NULL )
  {
    ProgressBar *topLevel = new ProgressBar(brickData.size(), "Generating image");

    buildTopLevel( outputImage, 0, outputHeight, outputWidth, outputHeight, brickData, brickImages, outputWidth, outputHeight, topLevel );
  }
  else
  {
    int level = (int)ceil(log2( max(outputWidth,outputHeight) ) );

    if( outputImage.back() == '/' ) outputImage = outputImage.substr(0, outputImage.size()-1);

    g_mkdir(string(outputImage).append("_files/").c_str(), 0777);
    g_mkdir(string(outputImage).append("_files/").append(to_string(level)).c_str(), 0777);

    int threads = sysconf(_SC_NPROCESSORS_ONLN);

    ProgressBar *topLevel = new ProgressBar(ceil((double)outputWidth/256.0)*ceil((double)outputHeight/((double)threads*256.0)), "Building top level");

    future< void > ret[threads];

    for( int k = 0; k < threads; ++k )
    {
      int start = k*outputHeight/threads;
      int end = (k+1)*outputHeight/threads;

      start = (start / 256 ) * 256;
      if( k+1 < threads ) end = (end / 256 ) * 256;

      ret[k] = async( launch::async, &buildTopLevel, string(outputImage).append("_files/"+to_string(level)+"/"), start, end, outputWidth, outputHeight, ref(brickData), ref(brickImages), 256, 256, topLevel );
    }

    // Wait for threads to finish
    for( int k = 0; k < threads; ++k )
    {
      ret[k].get();
    }
    
    topLevel->Finish();

    ofstream dzi_file;
    dzi_file.open(string(outputImage).append(".dzi").c_str());
    dzi_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
    dzi_file << "<Image xmlns=\"http://schemas.microsoft.com/deepzoom/2008\" Format=\"jpeg\" Overlap=\"0\" TileSize=\"256\">" << endl;
    dzi_file << "    <Size Height=\"" << outputHeight << "\" Width=\"" << outputWidth << "\"/>" << endl;
    dzi_file << "</Image>" << endl;
    dzi_file.close();

    int numLevels = 0;
    for( int o = ceil((double)outputWidth/256.0); o > 1; o = ceil((double)o/2.0) ) numLevels += o;

    ProgressBar *lowerLevels = new ProgressBar(numLevels, "Building lower levels");

    outputWidth = (int)ceil((double)outputWidth/ 128.0 );
    outputHeight = (int)ceil((double)outputHeight/ 128.0 );

    for( ; level > 0; --level )
    {
      outputWidth = (int)ceil((double)outputWidth/ 2.0 );
      outputHeight = (int)ceil((double)outputHeight/ 2.0 );

      string current = string(outputImage).append("_files/"+to_string(level-1)+"/");
      string upper = string(outputImage).append("_files/"+to_string(level)+"/");

      g_mkdir((char *)current.c_str(), 0777);

      for(int i = 1; i < outputWidth; i+=2)
      {
        lowerLevels->Increment();

        for(int j = 1; j < outputHeight; j+=2)
        {
          (VImage::jpegload((char *)string(upper).append(to_string(i-1)+"_"+to_string(j-1)+".jpeg").c_str(), VImage::option()->set( "shrink", 2)).
          join(VImage::jpegload((char *)string(upper).append(to_string(i)+"_"+to_string(j-1)+".jpeg").c_str(), VImage::option()->set( "shrink", 2)),VIPS_DIRECTION_HORIZONTAL)).
          join(VImage::jpegload((char *)string(upper).append(to_string(i-1)+"_"+to_string(j)+".jpeg").c_str(), VImage::option()->set( "shrink", 2)).
          join(VImage::jpegload((char *)string(upper).append(to_string(i)+"_"+to_string(j)+".jpeg").c_str(), VImage::option()->set( "shrink", 2)),VIPS_DIRECTION_HORIZONTAL),VIPS_DIRECTION_VERTICAL).
          jpegsave((char *)string(current).append(to_string(i>>1)+"_"+to_string(j>>1)+".jpeg").c_str() );
        }
      }
      if(outputWidth%2 == 1)
      {
        for(int j = 1; j < outputHeight; j+=2)
        {
          (VImage::jpegload((char *)string(upper).append(to_string(outputWidth-1)+"_"+to_string(j-1)+".jpeg").c_str(), VImage::option()->set( "shrink", 2)).
          join(VImage::jpegload((char *)string(upper).append(to_string(outputWidth-1)+"_"+to_string(j)+".jpeg").c_str(), VImage::option()->set( "shrink", 2)),VIPS_DIRECTION_VERTICAL)).
          jpegsave((char *)string(current).append(to_string(outputWidth>>1)+"_"+to_string(j>>1)+".jpeg").c_str(), VImage::option()->set( "optimize_coding", true )->set( "strip", true ) );
        }
      }
      if(outputHeight%2 == 1)
      {
        for(int j = 1; j < outputWidth; j+=2)
        {
          (VImage::jpegload((char *)string(upper).append(to_string(j-1)+"_"+to_string(outputHeight-1)+".jpeg").c_str(), VImage::option()->set( "shrink", 2)).
          join(VImage::jpegload((char *)string(upper).append(to_string(j)+"_"+to_string(outputHeight-1)+".jpeg").c_str(), VImage::option()->set( "shrink", 2)),VIPS_DIRECTION_HORIZONTAL)).
          jpegsave((char *)string(current).append(to_string(j>>1)+"_"+to_string(outputHeight>>1)+".jpeg").c_str(), VImage::option()->set( "optimize_coding", true )->set( "strip", true ) );
        }
      }
      if(outputWidth%2 == 1 && outputHeight%2 == 1)
      {
        VImage::jpegload((char *)string(upper).append(to_string(outputWidth-1)+"_"+to_string(outputHeight-1)+".jpeg").c_str()).resize(0.5).
        jpegsave((char *)string(current).append(to_string(outputWidth>>1)+"_"+to_string(outputHeight>>1)+".jpeg").c_str(), VImage::option()->set( "optimize_coding", true )->set( "strip", true ) );  
      }
    }

    lowerLevels->Finish();

    // Generate html file to view deep zoom image
    ofstream htmlFile(string(outputImage).append(".html").c_str());
    htmlFile << "<!DOCTYPE html>\n<html>\n<head><script src=\"js/openseadragon.min.js\"></script></head>\n<body>\n<style>\nhtml,\nbody,\n#collage\n{\nposition: fixed;\nleft: 0;\ntop: 0;\nwidth: 100%;\nheight: 100%;\n}\n</style>\n\n<div id=\"collage\"></div>\n\n<script>\nvar viewer = OpenSeadragon({\nid: 'collage',\nprefixUrl: 'icons/',\ntileSources:   \"" + outputImage + ".dzi\",\nminZoomImageRatio: 0,\nmaxZoomImageRatio: 1\n});\n</script>\n</body>\n</html>";
    htmlFile.close();
  }
}

void generateOutput( vector< vector< int > > mosaic, string outputImage )
{
  vector< vector< vector< unsigned char > > > brickImages(25);
  vector< pair< int, int > > dimensions(25);

  generateImages( brickImages, dimensions );

  size_t lastindex = outputImage.find_last_of("."); 
  string outputFile = outputImage.substr(0, lastindex).append(".ldr"); 

  ofstream output( outputFile );

  int width = mosaic[0].size();
  int height = mosaic.size();

  int width2 = (int)ceil(width*round(77*si));
  int height2 = (int)ceil(height*round(77*si));

  vector< vector< int > > brickData;

  for( int brick = 0; brick < numBricks; ++brick )
  {
    int w = brickSizes[brick][0];
    int h = brickSizes[brick][1];

    vector< int > indices( (width-w+1)*(height-h+1) );
    iota( indices.begin(), indices.end(), 0 );

    // Shuffle the points so that patterens do not form
    shuffle( indices.begin(), indices.end(), default_random_engine(0));

    for( int k = 0; k < indices.size(); ++k )
    {
      int i = indices[k]/(width-w+1);
      int j = indices[k]%(width-w+1);

      int squareColor = mosaic[i][j];

      if( squareColor < 0 ) continue;

      bool badColor = false;

      for( int y = i; y < i + h && !badColor; ++y )
      {
        for( int x = j; x < j + w && !badColor; ++x )
        {
          if( mosaic[y][x] != squareColor )
          {
            badColor = true;
          }
        }
      }

      if( !badColor )
      {
        for( int y = i; y < i + h; ++y )
        {
          for( int x = j; x < j + w; ++x )
          {
            mosaic[y][x] = -1;
          }
        }

        int bh = dimensions[brick].second;
        int bw = dimensions[brick].first;

        int bi = round(77*si)*i;
        int bj = round(77*si)*j;

        int colorIndex = 0;

        for( int l = 0; l < 68; ++l )
        {
          if( color_codes[l] == squareColor )
          {
            colorIndex = l;
            break;
          }
        }

        brickData.push_back({brick,colorIndex,bh,bw,bi,bj});

        output << " 1 " << squareColor << " " << j*20+w*10 << " 0 " << i*20+h*10 << brickNames[brick] << endl;
      }
    }
  }

  for( int i = 0; i < height; ++i )
  {
    for( int j = 0; j < width; ++j )
    {
      if( mosaic[i][j] >= 0 )
      {
        output << " 1 " << mosaic[i][j] << " " << j*20+10 << " 0 " << i*20+10 << " 1 0 0 0 1 0 0 0 1 3024.dat\n";

        int bh = dimensions[24].second;
        int bw = dimensions[24].first;

        int bi = round(77*si)*i;
        int bj = round(77*si)*j;

        int colorIndex = 0;

        for( int l = 0; l < 68; ++l )
        {
          if( color_codes[l] == mosaic[i][j] )
          {
            colorIndex = l;
            break;
          }
        }

        brickData.push_back({24,colorIndex,bh,bw,bi,bj});
      }
    }
  }

  output.close();

  buildImage( outputImage, width2, height2, brickData, brickImages );
}

double de00( int r, int g, int b, int rv, int gv, int bv )
{
  float l1,a1,b1,l2,a2,b2;
  // Convert rgb to rgb16
  vips_col_sRGB2scRGB_8( r, g, b, &l1,&a1,&b1 );
  vips_col_sRGB2scRGB_8( rv, gv, bv, &l2,&a2,&b2 );
  // Convert rgb16 to xyz
  vips_col_scRGB2XYZ( l1, a1, b1, &l1, &a1, &b1 );
  vips_col_scRGB2XYZ( l2, a2, b2, &l2, &a2, &b2 );
  // Convert xyz to lab
  vips_col_XYZ2Lab( l1, a1, b1, &l1, &a1, &b1 );
  vips_col_XYZ2Lab( l2, a2, b2, &l2, &a2, &b2 );

  return vips_col_dE00( l1, a1, b1, l2, a2, b2 );
}

void generateLegoMosaic( string inputImage, string outputImage, int numAcross, bool dither, double gamma, vector< int > &colorsToUse )
{
  // Load input image
  VImage image = VImage::thumbnail((char *)inputImage.c_str(),numAcross,VImage::option()->set( "auto_rotate", true ));//.flip(VIPS_DIRECTION_HORIZONTAL);

  // Convert to a three band image
  if( image.bands() == 1 )
  {
    image = image.bandjoin(image).bandjoin(image);
  }
  if( image.bands() == 4 )
  {
    image = image.flatten();
  }

  int width = image.width();
  int height = image.height();

  // Get image data
  unsigned char * inputData = ( unsigned char * )image.data();
  float * c2 = new float[3*width*height];

  vector< vector< vector< float > > > floatData(height,vector< vector< float > >(width,vector< float >(3,0)));

  int threads = sysconf(_SC_NPROCESSORS_ONLN);

  ProgressBar *changingColorspace = new ProgressBar(height/threads, "Changing colorspace");

  vector<Point> points;

  double center1[3] = {0.0,0.0,0.0};

  float colors[68][3];

  int n100 = 0;

  for( int k1 = 0; k1 < colorsToUse.size(); ++k1 )
  {
    int k = colorsToUse[k1];

    ++n100;

    rgbToLab(color_values[k][0],color_values[k][1],color_values[k][2],colors[k][0],colors[k][1],colors[k][2]);

    points.push_back(Point(colors[k][0],colors[k][1],colors[k][2]));

    for( int i1 = 0; i1 < 3; ++i1 )
    {
      center1[i1] += colors[k][i1];
    }
  }

  Point center(center1[0]/(float)n100,center1[1]/(float)n100,center1[2]/(float)n100);

  Polyhedron poly;

  CGAL::convex_hull_3(points.begin(), points.end(), poly);

  Tree tree(faces(poly).first, faces(poly).second, poly);

  future< void > ret2[threads];

  for( int k = 0; k < threads; ++k )
  {
    ret2[k] = async( launch::async, &changeColorspace, ref(tree), ref(center), inputData, c2, k*height/threads, (k+1)*height/threads, width, gamma, true, changingColorspace, true, false );
  }

  // Wait for threads to finish
  for( int k = 0; k < threads; ++k )
  {
    ret2[k].get();
  }

  changingColorspace->Finish();

  for( int i = 0, index = 0; i < height; ++i )
  {
    for( int j = 0; j < width; ++j, ++index )
    {
      for( int k = 0; k < 3; ++k )
      {
        floatData[i][j][k] = c2[3*index+k];
      }
    }
  }

  vector< vector< int > > mosaic(height,vector< int >(width,-1));

  ProgressBar *buildingMosaic = new ProgressBar(height, "Building Mosaic");

  for( int i = 0, index = 0; i < height; ++i )
  {
    for( int j = 0; j < width; ++j, ++index )
    {
      int bestColor = -1;
      double bestDifference = 1000000000;

      for( int k1 = 0; k1 < colorsToUse.size(); ++k1 )
      {
        int k = colorsToUse[k1];

        double difference = vips_col_dE00( colors[k][0], colors[k][1], colors[k][2], floatData[i][j][0], floatData[i][j][1], floatData[i][j][2] );

        if( difference < bestDifference )
        {
          bestDifference = difference;
          bestColor = k;
        }
      }

      if( dither )
      {
        float error[3] = {floatData[i][j][0]-colors[bestColor][0],floatData[i][j][1]-colors[bestColor][1],floatData[i][j][2]-colors[bestColor][2]};
      
        if( j < width - 1 )
        {
          for( int i1 = 0; i1 < 3; ++i1 )
          {
            floatData[i][j+1][i1] += error[i1] * 7.0 / 16.0;
            fit(floatData[i][j+1][i1],i1);
          }
        }

        if( i < height - 1 )
        {
          if( j > 0 )
          {
            for( int i1 = 0; i1 < 3; ++i1 )
            {
              floatData[i+1][j-1][i1] += error[i1] * 3.0 / 16.0;
              fit(floatData[i+1][j-1][i1],i1);
            }
          }

          {
            for( int i1 = 0; i1 < 3; ++i1 )
            {
              floatData[i+1][j][i1] += error[i1] * 5.0 / 16.0;
              fit(floatData[i+1][j][i1],i1);
            }
          }

          if( j < width - 1 )
          {
            for( int i1 = 0; i1 < 3; ++i1 )
            {
              floatData[i+1][j+1][i1] += error[i1] * 1.0 / 16.0;
              fit(floatData[i+1][j+1][i1],i1);
            }
          }
        }
      }

      mosaic[i][j] = color_codes[bestColor];

    }

    buildingMosaic->Increment();
  }

  buildingMosaic->Finish();

  generateOutput( mosaic, outputImage );
}

void generateImages2( vector< vector< vector< unsigned char > > > &brickImages, vector< pair< int, int > > &dimensions )
{
  ifstream shapeData( "legoData/sidesOut.dat", ios::binary );

  size_t len;

  ProgressBar *generatingImages = new ProgressBar( 18*68, "Loading images" );

  int index10 = 0;

  for( int i = 0; i < 18; ++i )
  {
    for( int j = 0; j < 68; ++j )
    {
      shapeData.read( (char *)&len, sizeof(size_t) );

      void * data = (void *) malloc (sizeof(unsigned char)*len);

      shapeData.read( (char *)data, len*sizeof(char) );

      VipsBlob *blob = vips_blob_new( NULL, data, len );

      VImage image = VImage::pngload_buffer(blob).resize(si);

      unsigned char * c = ( unsigned char * )image.data();
      dimensions[i] = pair< int, int >( image.width(), image.height() );
      brickImages[i].push_back( vector< unsigned char >(c, c + 3*image.width()*image.height() ) );

      c = (unsigned char *)image.rot(VIPS_ANGLE_D270).data();
      dimensions[i+18] = pair< int, int >( image.height(), image.width() );
      brickImages[i+18].push_back( vector< unsigned char >(c, c + 3*image.width()*image.height() ) );   

      generatingImages->Increment();
    }
  }

  generatingImages->Finish();
}

void generateOutput2( vector< vector< int > > &mosaic1, vector< vector< int > > &mosaic2, string outputImage )
{
  vector< vector< vector< unsigned char > > > brickImages(36);
  vector< pair< int, int > > dimensions(36);

  generateImages2( brickImages, dimensions );

  size_t lastindex = outputImage.find_last_of("."); 
  string outputFile = outputImage.substr(0, lastindex).append(".ldr"); 

  ofstream output( outputFile );

  int width = mosaic1[0].size();
  int height = mosaic1.size();

  int width2 = (int)ceil(width*round(77*si));
  int height2 = (int)ceil(mosaic2.size()*round(77*si));

  vector< vector< int > > brickData;

  for( int brick = 0; brick < 16; ++brick )
  {
    int w = brickSizes2[brick][0];
    int h = brickSizes2[brick][1];

    vector< int > indices( (width-w+1)*(height-h+1) );
    iota( indices.begin(), indices.end(), 0 );

    // Shuffle the points so that patterens do not form
    shuffle( indices.begin(), indices.end(), default_random_engine(0));

    for( int k = 0; k < indices.size(); ++k )
    {
      int i = indices[k]/(width-w+1);
      int j = indices[k]%(width-w+1);

      int squareColor = mosaic1[i][j];

      if( squareColor < 0 ) continue;

      bool badColor = false;

      for( int x = j; x < j + w && !badColor; ++x )
      {
        if( brick < 11 && ( i < 1 || mosaic1[i-1][x] == -1 ) )
        {
          badColor = true; 
        }
      }

      for( int y = i; y < i + h && !badColor; ++y )
      {
        for( int x = j; x < j + w && !badColor; ++x )
        {
          if( mosaic1[y][x] != squareColor )
          {
            badColor = true;
          }
        }
      }

      if( !badColor )
      {
        for( int y = i; y < i + h; ++y )
        {
          for( int x = j; x < j + w; ++x )
          {
            mosaic1[y][x] = -2;
          }
        }

        int bh = dimensions[brick].second;
        int bw = dimensions[brick].first;

        int bi = floor(30.8*si*i);
        int bj = floor(77.0*si*j);

        int colorIndex = 0;

        for( int l = 0; l < 68; ++l )
        {
          if( color_codes[l] == squareColor )
          {
            colorIndex = l;
            break;
          }
        }

        brickData.push_back({brick,colorIndex,bh,bw,bi,bj});

        output << " 1 " << squareColor << " " << j*20+w*10 << " " << i*8 << " 0" << " 1 0 0 0 1 0 0 0 1 " << brickNames3[brick] << endl;
      }
    }
  }

  for( int i = 0; i < height; ++i )
  {
    for( int j = 0; j < width; ++j )
    {
      if( mosaic1[i][j] >= 0 )
      {
        if( i < 1 || mosaic1[i-1][j] == -1 )
        {
          int bh = dimensions[16].second;
          int bw = dimensions[16].first;

          int bi = floor(30.8*si*i);
          int bj = floor(77.0*si*j);

          int colorIndex = 0;

          for( int l = 0; l < 68; ++l )
          {
            if( color_codes[l] == mosaic1[i][j] )
            {
              colorIndex = l;
              break;
            }
          }

          brickData.push_back({16,colorIndex,bh,bw,bi,bj});

          output << " 1 " << mosaic1[i][j] << " " << j*20+10 << " " << i*8 << " 0" << " 1 0 0 0 1 0 0 0 1 30039.dat\n";
        }
        else
        {
          int bh = dimensions[17].second;
          int bw = dimensions[17].first;

          int bi = floor(30.8*si*i);
          int bj = floor(77.0*si*j);

          int colorIndex = 0;

          for( int l = 0; l < 68; ++l )
          {
            if( color_codes[l] == mosaic1[i][j] )
            {
              colorIndex = l;
              break;
            }
          }

          brickData.push_back({17,colorIndex,bh,bw,bi,bj});

          output << " 1 " << mosaic1[i][j] << " " << j*20+10 << " " << i*8 << " 0" << " 1 0 0 0 1 0 0 0 1 3024.dat\n";
        }
      }
    }
  }

  width = mosaic2[0].size();
  height = mosaic2.size();

  for( int brick = 0; brick < 16; ++brick )
  {
    int w = brickSizes2[brick][1];
    int h = brickSizes2[brick][0];

    vector< int > indices( (width-w+1)*(height-h+1) );
    iota( indices.begin(), indices.end(), 0 );

    // Shuffle the points so that patterens do not form
    shuffle( indices.begin(), indices.end(), default_random_engine(0));

    for( int k = 0; k < indices.size(); ++k )
    {
      int i = indices[k]/(width-w+1);
      int j = indices[k]%(width-w+1);

      int squareColor = mosaic2[i][j];

      if( squareColor < 0 ) continue;

      bool badColor = false;

      for( int y = i; y < i + h && !badColor; ++y )
      {
        if( brick < 11 && ( j < 1 || mosaic2[y][j-1] == -1 ) )
        {
          badColor = true; 
        }
      }

      for( int y = i; y < i + h && !badColor; ++y )
      {
        for( int x = j; x < j + w && !badColor; ++x )
        {
          if( mosaic2[y][x] != squareColor )
          {
            badColor = true;
          }
        }
      }

      if( !badColor )
      {
        for( int y = i; y < i + h; ++y )
        {
          for( int x = j; x < j + w; ++x )
          {
            mosaic2[y][x] = -2;
          }
        }


        int bh = dimensions[brick+18].second;
        int bw = dimensions[brick+18].first;

        int bi = floor(77.0*si*i);
        int bj = floor(30.8*si*j);

        int colorIndex = 0;

        for( int l = 0; l < 68; ++l )
        {
          if( color_codes[l] == squareColor )
          {
            colorIndex = l;
            break;
          }
        }

        brickData.push_back({brick+18,colorIndex,bh,bw,bi,bj});

        output << " 1 " << squareColor << " " << j*8 << " " << i*20+h*10 << " 0" << " 0 1 0 -1 0 0 0 0 1 " << brickNames3[brick] << endl;
      }
    }
  }

  for( int i = 0; i < height; ++i )
  {
    for( int j = 0; j < width; ++j )
    {
      if( mosaic2[i][j] >= 0 )
      {
        if( j < 1 || mosaic2[i][j-1] == -1 )
        {
          int bh = dimensions[34].second;
          int bw = dimensions[34].first;

          int bi = floor(77.0*si*i);
          int bj = floor(30.8*si*j);

          int colorIndex = 0;

          for( int l = 0; l < 68; ++l )
          {
            if( color_codes[l] == mosaic2[i][j] )
            {
              colorIndex = l;
              break;
            }
          }

          brickData.push_back({34,colorIndex,bh,bw,bi,bj});

          output << " 1 " << mosaic2[i][j] << " " << j*8 << " " << i*20+10 << " 0" << " 0 1 0 -1 0 0 0 0 1 30039.dat\n";
        }
        else
        {
          int bh = dimensions[35].second;
          int bw = dimensions[35].first;

          int bi = floor(77.0*si*i);
          int bj = floor(30.8*si*j);

          int colorIndex = 0;

          for( int l = 0; l < 68; ++l )
          {
            if( color_codes[l] == mosaic2[i][j] )
            {
              colorIndex = l;
              break;
            }
          }

          brickData.push_back({35,colorIndex,bh,bw,bi,bj});

          output << " 1 " << mosaic2[i][j] << " " << j*8 << " " << i*20+10 << " 0" << " 0 1 0 -1 0 0 0 0 1 3024.dat\n";
        }
      }
    }
  }

  output.close();

  buildImage( outputImage, width2, height2, brickData, brickImages );
}

void bestSquare( int i1, int j1, int width, int height, vector< vector< vector< float > > > &floatData, vector< vector< vector< float > > > &floatData2, vector< vector< int > > &mosaic1, vector< vector< int > > &mosaic2, bool dither, float colors[][3], vector< int > colorsToUse )
{
  vector< vector< float > > tempFloatData, tempFloatData1, tempFloatData2;

  stringstream ss1, ss2;

  double sumDifferences1 = 0;
  double sumDifferences2 = 0;

  vector< vector< int > > results1;
  vector< vector< int > > results2;

  for( int i = 0; i < 5; ++i )
  {
    for( int j = 0; j < 2; ++j )
    {
      int bestColor = -1;
      double bestDifference = 1000000000;
      double bestDifference2 = 1000000000;
      
      float averageColor[3] = {0,0,0};

      for( int y = i*2; y < i*2+2; ++y )
      {
        for( int x = j*5; x < j*5+5; ++x )
        {
          for( int i10 = 0; i10 < 3; ++i10 )
          {
            averageColor[i10] += floatData[i1+y][j1+x][i10];
          }
        }
      }

      for( int i10 = 0; i10 < 3; ++i10 )
      {
        averageColor[i10] /= 10.0;
      }

      for( int k1 = 0; k1 < colorsToUse.size(); ++k1 )
      {
        int k = colorsToUse[k1];

        int r = 0;
        int g = 0;
        int b = 0;

        double difference = 0;

        for( int y = i*2; y < i*2+2; ++y )
        {
          for( int x = j*5; x < j*5+5; ++x )
          {
            difference += de00( color_values[k][0], color_values[k][1], color_values[k][2], round(floatData2[i1+y][j1+x][0]), round(floatData2[i1+y][j1+x][1]), round(floatData2[i1+y][j1+x][2]) );
          }
        }

        double difference2 = vips_col_dE00( colors[k][0], colors[k][1], colors[k][2], averageColor[0], averageColor[1], averageColor[2] );

        if( difference2 < bestDifference2 )
        {
          bestDifference2 = difference2;
          bestColor = k;
        }

        if( difference < bestDifference )
        {
          bestDifference = difference;
        }
      }

      if( dither )
      {
        float error[3];

        for( int i10 = 0; i10 < 3; ++i10 )
        {
          error[i10] = averageColor[i10]-colors[bestColor][i10];
        }

        if( j1+j*5 < width - 5 )
        {
          for( int y = i*2; y < i*2+2; ++y )
          {
            for( int x = j*5; x < j*5+5; ++x )
            {
              int x1 = j1+x+5;
              int y1 = i1+y;

              tempFloatData.push_back({(float)x1,(float)y1,floatData[y1][x1][0],floatData[y1][x1][1],floatData[y1][x1][2]});

              for( int i10 = 0; i10 < 3; ++i10 )
              {
                floatData[y1][x1][i10] += error[i10] * 7.0 / 16.0;
                fit(floatData[y1][x1][i10],i10);
              }

              tempFloatData1.push_back({(float)x1,(float)y1,floatData[y1][x1][0],floatData[y1][x1][1],floatData[y1][x1][2]});
            }
          }
        }

        if( i1+i*2 < height - 2 )
        {
          if( j1+j > 0 )
          {
            for( int y = i*2; y < i*2+2; ++y )
            {
              for( int x = j*5; x < j*5+5; ++x )
              {
                int x1 = j1+x-5;
                int y1 = i1+y+2;

                tempFloatData.push_back({(float)x1,(float)y1,floatData[y1][x1][0],floatData[y1][x1][1],floatData[y1][x1][2]});

                for( int i10 = 0; i10 < 3; ++i10 )
                {
                  floatData[y1][x1][i10] += error[i10] * 3.0 / 16.0;
                  fit(floatData[y1][x1][i10],i10);
                }

                tempFloatData1.push_back({(float)x1,(float)y1,floatData[y1][x1][0],floatData[y1][x1][1],floatData[y1][x1][2]});
              }
            }
          }

          {
            for( int y = i*2; y < i*2+2; ++y )
            {
              for( int x = j*5; x < j*5+5; ++x )
              {
                int x1 = j1+x;
                int y1 = i1+y+2;

                tempFloatData.push_back({(float)x1,(float)y1,floatData[y1][x1][0],floatData[y1][x1][1],floatData[y1][x1][2]});

                for( int i10 = 0; i10 < 3; ++i10 )
                {
                  floatData[y1][x1][i10] += error[i10] * 5.0 / 16.0;
                  fit(floatData[y1][x1][i10],i10);
                }

                tempFloatData1.push_back({(float)x1,(float)y1,floatData[y1][x1][0],floatData[y1][x1][1],floatData[y1][x1][2]});
              }
            }
          }

          if( j1+j*5 < width - 5 )
          {
            for( int y = i*2; y < i*2+2; ++y )
            {
              for( int x = j*5; x < j*5+5; ++x )
              {
                int x1 = j1+x+5;
                int y1 = i1+y+2;

                tempFloatData.push_back({(float)x1,(float)y1,floatData[y1][x1][0],floatData[y1][x1][1],floatData[y1][x1][2]});

                for( int i10 = 0; i10 < 3; ++i10 )
                {
                  floatData[y1][x1][i10] += error[i10] * 1.0 / 16.0;
                  fit(floatData[y1][x1][i10],i10);
                }

                tempFloatData1.push_back({(float)x1,(float)y1,floatData[y1][x1][0],floatData[y1][x1][1],floatData[y1][x1][2]});
              }
            }
          }
        }
      }

      sumDifferences1 += bestDifference;

      results1.push_back( {j1/5+j,i1/2+i,color_codes[bestColor]} );
    }
  }

  for( int i = tempFloatData.size()-1; i >= 0; --i )
  {
    floatData[(int)tempFloatData[i][1]][(int)tempFloatData[i][0]][0] = tempFloatData[i][2];
    floatData[(int)tempFloatData[i][1]][(int)tempFloatData[i][0]][1] = tempFloatData[i][3];
    floatData[(int)tempFloatData[i][1]][(int)tempFloatData[i][0]][2] = tempFloatData[i][4];
  }

  tempFloatData.clear();

  for( int i = 0; i < 2; ++i )
  {
    for( int j = 0; j < 5; ++j )
    {
      int bestColor = -1;
      double bestDifference = 1000000000;
      double bestDifference2 = 1000000000;
      
      float averageColor[3] = {0,0,0};

      for( int y = i*5; y < i*5+5; ++y )
      {
        for( int x = j*2; x < j*2+2; ++x )
        {
          for( int i10 = 0; i10 < 3; ++i10 )
          {
            averageColor[i10] += floatData[i1+y][j1+x][i10];
          }
        }
      }

      for( int i10 = 0; i10 < 3; ++i10 )
      {
        averageColor[i10] /= 10.0;
      }

      for( int k1 = 0; k1 < colorsToUse.size(); ++k1 )
      {
        int k = colorsToUse[k1];

        int r = 0;
        int g = 0;
        int b = 0;

        double difference = 0;

        for( int y = i*5; y < i*5+5; ++y )
        {
          for( int x = j*2; x < j*2+2; ++x )
          {
            difference += de00( color_values[k][0], color_values[k][1], color_values[k][2], round(floatData2[i1+y][j1+x][0]), round(floatData2[i1+y][j1+x][1]), round(floatData2[i1+y][j1+x][2]) );
          }
        }

        double difference2 = vips_col_dE00( colors[k][0], colors[k][1], colors[k][2], averageColor[0], averageColor[1], averageColor[2] );

        if( difference2 < bestDifference2 )
        {
          bestDifference2 = difference2;
          bestColor = k;
        }

        if( difference < bestDifference )
        {
          bestDifference = difference;
        }
      }

      if( dither )
      {
        float error[3];

        for( int i10 = 0; i10 < 3; ++i10 )
        {
          error[i10] = averageColor[i10]-colors[bestColor][i10];
        }

        if( j1+j*2 < width - 2 )
        {
          for( int y = i*5; y < i*5+5; ++y )
          {
            for( int x = j*2; x < j*2+2; ++x )
            {
              int x1 = j1+x+2;
              int y1 = i1+y;

              tempFloatData.push_back({(float)x1,(float)y1,floatData[y1][x1][0],floatData[y1][x1][1],floatData[y1][x1][2]});

              for( int i10 = 0; i10 < 3; ++i10 )
              {
                floatData[y1][x1][i10] += error[i10] * 7.0 / 16.0;
                fit(floatData[y1][x1][i10],i10);
              }

              tempFloatData2.push_back({(float)x1,(float)y1,floatData[y1][x1][0],floatData[y1][x1][1],floatData[y1][x1][2]});
            }
          }
        }

        if( i1+i*5 < height - 5 )
        {
          if( j1+j > 0 )
          {
            for( int y = i*5; y < i*5+5; ++y )
            {
              for( int x = j*2; x < j*2+2; ++x )
              {
                int x1 = j1+x-2;
                int y1 = i1+y+5;

                tempFloatData.push_back({(float)x1,(float)y1,floatData[y1][x1][0],floatData[y1][x1][1],floatData[y1][x1][2]});

                for( int i10 = 0; i10 < 3; ++i10 )
                {
                  floatData[y1][x1][i10] += error[i10] * 3.0 / 16.0;
                  fit(floatData[y1][x1][i10],i10);
                }

                tempFloatData2.push_back({(float)x1,(float)y1,floatData[y1][x1][0],floatData[y1][x1][1],floatData[y1][x1][2]});
              }
            }
          }

          {
            for( int y = i*5; y < i*5+5; ++y )
            {
              for( int x = j*2; x < j*2+2; ++x )
              {
                int x1 = j1+x;
                int y1 = i1+y+5;

                tempFloatData.push_back({(float)x1,(float)y1,floatData[y1][x1][0],floatData[y1][x1][1],floatData[y1][x1][2]});

                for( int i10 = 0; i10 < 3; ++i10 )
                {
                  floatData[y1][x1][i10] += error[i10] * 5.0 / 16.0;
                  fit(floatData[y1][x1][i10],i10);
                }

                tempFloatData2.push_back({(float)x1,(float)y1,floatData[y1][x1][0],floatData[y1][x1][1],floatData[y1][x1][2]});
              }
            }
          }

          if( j1+j*2 < width - 2 )
          {
            for( int y = i*5; y < i*5+5; ++y )
            {
              for( int x = j*2; x < j*2+2; ++x )
              {
                int x1 = j1+x+2;
                int y1 = i1+y+5;

                tempFloatData.push_back({(float)x1,(float)y1,floatData[y1][x1][0],floatData[y1][x1][1],floatData[y1][x1][2]});

                for( int i10 = 0; i10 < 3; ++i10 )
                {
                  floatData[y1][x1][i10] += error[i10] * 1.0 / 16.0;
                  fit(floatData[y1][x1][i10],i10);
                }

                tempFloatData2.push_back({(float)x1,(float)y1,floatData[y1][x1][0],floatData[y1][x1][1],floatData[y1][x1][2]});
              }
            }
          }
        }
      }

      sumDifferences2 += bestDifference;

      results2.push_back( {j1/2+j,i1/5+i,color_codes[bestColor]} );
    }
  }

  for( int i = tempFloatData.size()-1; i >= 0; --i )
  {
    floatData[(int)tempFloatData[i][1]][(int)tempFloatData[i][0]][0] = tempFloatData[i][2];
    floatData[(int)tempFloatData[i][1]][(int)tempFloatData[i][0]][1] = tempFloatData[i][3];
    floatData[(int)tempFloatData[i][1]][(int)tempFloatData[i][0]][2] = tempFloatData[i][4];
  }

  if( true && sumDifferences1 <= sumDifferences2 )
  {
    for( int i = 0; i < 10; ++i )
    {
      mosaic1[results1[i][1]][results1[i][0]] = results1[i][2];
    }

    for( int i = 0; i < tempFloatData1.size(); ++i )
    {
      floatData[(int)tempFloatData1[i][1]][(int)tempFloatData1[i][0]][0] = tempFloatData1[i][2];
      floatData[(int)tempFloatData1[i][1]][(int)tempFloatData1[i][0]][1] = tempFloatData1[i][3];
      floatData[(int)tempFloatData1[i][1]][(int)tempFloatData1[i][0]][2] = tempFloatData1[i][4];
    }
  }
  else
  {
    for( int i = 0; i < 10; ++i )
    {
      mosaic2[results2[i][1]][results2[i][0]] = results2[i][2];
    }

    for( int i = 0; i < tempFloatData2.size(); ++i )
    {
      floatData[(int)tempFloatData2[i][1]][(int)tempFloatData2[i][0]][0] = tempFloatData2[i][2];
      floatData[(int)tempFloatData2[i][1]][(int)tempFloatData2[i][0]][1] = tempFloatData2[i][3];
      floatData[(int)tempFloatData2[i][1]][(int)tempFloatData2[i][0]][2] = tempFloatData2[i][4];
    }
  }
}

void generateLegoMosaic2( string inputImage, string outputFile, int numAcross, bool dither, double gamma, vector< int > &colorsToUse )
{
  // Load input image
  VImage image = VImage::thumbnail((char *)inputImage.c_str(),(numAcross-numAcross%2)*5,VImage::option()->set( "auto_rotate", true ));//.flip(VIPS_DIRECTION_HORIZONTAL);

  // Convert to a three band image
  if( image.bands() == 1 )
  {
    image = image.bandjoin(image).bandjoin(image);
  }
  if( image.bands() == 4 )
  {
    image = image.flatten();
  }

  int width = image.width();
  int height = image.height();

  width -= width%10;
  height -= height%10;

  image = image.extract_area(0,0,width,height);

  // Get image data
  unsigned char * inputData = ( unsigned char * )image.data();
  float * c2 = new float[3*width*height];

  vector< vector< vector< float > > > floatData(height,vector< vector< float > >(width,vector< float >(3,0))), floatData2;

  int threads = sysconf(_SC_NPROCESSORS_ONLN);

  ProgressBar *changingColorspace = new ProgressBar(height/threads, "Changing colorspace");

  vector<Point> points;

  double center1[3] = {0.0,0.0,0.0};

  float colors[68][3];

  int n100 = 0;

  for( int k1 = 0; k1 < colorsToUse.size(); ++k1 )
  {
    int k = colorsToUse[k1];

    ++n100;

    rgbToLab(color_values[k][0],color_values[k][1],color_values[k][2],colors[k][0],colors[k][1],colors[k][2]);

    points.push_back(Point(colors[k][0],colors[k][1],colors[k][2]));

    for( int i1 = 0; i1 < 3; ++i1 )
    {
      center1[i1] += colors[k][i1];
    }
  }

  Point center(center1[0]/(float)n100,center1[1]/(float)n100,center1[2]/(float)n100);

  Polyhedron poly;

  CGAL::convex_hull_3(points.begin(), points.end(), poly);

  Tree tree(faces(poly).first, faces(poly).second, poly);

  future< void > ret2[threads];

  for( int k = 0; k < threads; ++k )
  {
    ret2[k] = async( launch::async, &changeColorspace, ref(tree), ref(center), inputData, c2, k*height/threads, (k+1)*height/threads, width, gamma, true, changingColorspace, true, false );
  }

  // Wait for threads to finish
  for( int k = 0; k < threads; ++k )
  {
    ret2[k].get();
  }

  changingColorspace->Finish();

  for( int i = 0, index = 0; i < height; ++i )
  {
    for( int j = 0; j < width; ++j, ++index )
    {
      for( int k = 0; k < 3; ++k )
      {
        floatData[i][j][k] = c2[3*index+k];
      }
    }
  }

  floatData2 = floatData;

  vector< vector< int > > mosaic1( height/2, vector< int > ( width/5, -1 ) );
  vector< vector< int > > mosaic2( height/5, vector< int > ( width/2, -1 ) );

  ProgressBar *buildingMosaic = new ProgressBar(height*width/100, "Building Mosaic");

  for( int i = 0; i < height; i+=10 )
  {
    for( int j = 0; j < width; j+=10 )
    {
      int bestColor = -1;
      double bestDifference = 1000000000;

      bestSquare( i, j, width, height, floatData, floatData2, mosaic1, mosaic2, dither, colors, colorsToUse );
  
      buildingMosaic->Increment();
    }
  }

  buildingMosaic->Finish();

  generateOutput2( mosaic1, mosaic2, outputFile );
}

void RunLegoMosaic( string inputName, string outputName, int numHorizontal, int tileSize, bool sidesOut, bool dither, double gamma, string colorName )
{
  tileSize = max( tileSize, 5 );
  tileSize = min( tileSize, 77 );

  si = (double)tileSize/77.0;

  vector< int > colorsToUse;

  ifstream file(colorName);

  string str; 
  while (std::getline(file, str))
  {
    int n = stoi(str);

    int index = -1;

    for( int i = 0; i < 68; ++i )
    {
      if( color_codes[i] == n )
      {
        index = i;
        break;
      }
    }

    colorsToUse.push_back(index);
  }

  if( !sidesOut )
  {
    generateLegoMosaic( inputName, outputName, numHorizontal, dither, gamma, colorsToUse );
  }
  else
  {
    generateLegoMosaic2( inputName, outputName, numHorizontal, dither, gamma, colorsToUse );    
  }
}