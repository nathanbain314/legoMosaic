#include "LegoMosaic.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

double si = 1;

int numberOfCPUS()
{
  int numCPU;

#ifdef WIN32
  SYSTEM_INFO sysinfo;
  GetSystemInfo(&sysinfo);
  numCPU = sysinfo.dwNumberOfProcessors;
#else
  numCPU = sysconf(_SC_NPROCESSORS_ONLN);
#endif

  return numCPU;
}

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

void dilate( int width, int height, unsigned char * c, unsigned char * c2, int dilateDistance )
{
  for( int i = 0; i < height; ++i )
  {
    for( int j = 0; j < width; ++j )
    {
      float bestL = -1;
      for( int i1 = -dilateDistance; i1 <= dilateDistance; ++i1 )
      {
        if( i+i1 < 0 || i+i1 >= height ) continue;
        for( int j1 = -dilateDistance; j1 <= dilateDistance; ++j1 )
        {
          if( j+j1 < 0 || j+j1 >= width ) continue;

          float r = c[3*((i+i1)*width+j+j1)+0];
          float g = c[3*((i+i1)*width+j+j1)+1];
          float b = c[3*((i+i1)*width+j+j1)+2];

          float l = 0.2126*r + 0.7152*g + 0.0722*b;

          if( l <= 0.0031308 )
          {
            l = 12.92 * l;
          }
          else
          {
            l = 1.055 * pow( l, 1.0/2.4 ) - 0.055;
          }

          if( bestL < 0 || l < bestL )
          {
            bestL = l;

            c2[3*(i*width+j)+0] = r;
            c2[3*(i*width+j)+1] = g;
            c2[3*(i*width+j)+2] = b;
          }
        }
      }
    }
  }
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

void generateImages( vector< vector< vector< unsigned char > > > &brickImages, vector< pair< int, int > > &dimensions, vector< vector< bool > > &brickUsed )
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

      // If both colors are used
      if( brickUsed[index][j] && ( rotImage[i] == 0 || brickUsed[index+1][j] ) )
      {
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
      }
      // First is used, second is not
      else if( index < 24 && brickUsed[index][j] && !brickUsed[index+1][j] )
      {
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

        if( rotImage[i] > 0 )
        {
          brickImages[index+1].resize(brickImages[index+1].size()+1);
        }
      }
      // First is not used, second is
      else if( index < 24 && !brickUsed[index][j] && brickUsed[index+1][j] )
      {
        void * data = (void *) malloc (sizeof(unsigned char)*len);

        shapeData.read( (char *)data, len*sizeof(char) );

        VipsBlob *blob = vips_blob_new( NULL, data, len );

        VImage image = VImage::pngload_buffer(blob).resize(si);

        unsigned char * c;

        brickImages[index].resize(brickImages[index].size()+1);

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
      }
      // Neither are used
      else
      {
        brickImages[index].resize(brickImages[index].size()+1);
        
        if( rotImage[i] > 0 )
        {
          brickImages[index+1].resize(brickImages[index+1].size()+1);
        }

        shapeData.seekg( len, ios::cur );
      }

      generatingImages->Increment();
    }

    index += rotImage[i] > 0 ? 2 : 1;
  }

  generatingImages->Finish();
}

void buildTopLevel( string outputImage, int start, int end, int outputWidth, int outputHeight, vector< vector< int > > &brickData, vector< pair< int, int > > &dimensions, vector< vector< vector< unsigned char > > > &brickImages, int maxTileWidth, int maxTileHeight, ProgressBar *topLevel )
{
  for( int tileOffsetY = start; tileOffsetY < end; tileOffsetY += maxTileHeight )
  {
    int tileHeight = min( outputHeight - tileOffsetY, maxTileHeight );

    for( int tileOffsetX = 0; tileOffsetX < outputWidth; tileOffsetX += maxTileWidth )
    {
      if( start == 0 ) topLevel->Increment();

      int tileWidth = min( outputWidth - tileOffsetX, maxTileWidth );

      unsigned char * tileData = ( unsigned char * )calloc (tileHeight*tileWidth*3,sizeof(unsigned char));

      for( int k = 0; k < brickData.size(); ++k )
      {
        int brick = brickData[k][0];
        int colorIndex = brickData[k][1];
        int bh = dimensions[brick].second;
        int bw = dimensions[brick].first;
        int bi = brickData[k][2];
        int bj = brickData[k][3];

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

      string tileOutputName = string(outputImage).append(to_string(tileOffsetX/256)+"_"+to_string(tileOffsetY/256)+".jpeg");

      VImage::new_from_memory( tileData, tileHeight*tileWidth*3, tileWidth, tileHeight, 3, VIPS_FORMAT_UCHAR ).jpegsave((char *)tileOutputName.c_str(), VImage::option()->set( "optimize_coding", true )->set( "strip", true ) );

      free( tileData );
    }
  }
}

void buildTopLevelThread( string outputImage, int start, int end, unsigned char * tileData, int outputWidth, int outputHeight, vector< vector< int > > &brickData, vector< pair< int, int > > &dimensions, vector< vector< vector< unsigned char > > > &brickImages, int maxTileWidth, int maxTileHeight, ProgressBar *topLevel )
{
  for( int k = start; k < end; ++k )
  {
    if( start == 0 && k%10 == 0 ) topLevel->Increment();

    int brick = brickData[k][0];
    int colorIndex = brickData[k][1];
    int bh = dimensions[brick].second;
    int bw = dimensions[brick].first;
    int bi = brickData[k][2];
    int bj = brickData[k][3];

    if( bj >= outputWidth || bj + bw < 0 || bi >= outputHeight || bi + bh < 0 ) continue;

    for( int y1 = 0; y1 < bh; ++y1 )
    {
      for( int x1 = 0; x1 < bw; ++x1 )
      {
        int x = x1 + bj;
        int y = y1 + bi;

        if( !(x >= 0 && x < outputWidth && y >= 0 && y < outputHeight) ) continue;

        int index1 = 3*(y*outputWidth+x);

        int index2 = 3*(y1*bw+x1);

        tileData[index1+0] = brickImages[brick][colorIndex][index2+0];
        tileData[index1+1] = brickImages[brick][colorIndex][index2+1];
        tileData[index1+2] = brickImages[brick][colorIndex][index2+2];
      }
    }
  }
}

void buildImage( string outputImage, int outputWidth, int outputHeight, vector< vector< int > > &brickData, vector< pair< int, int > > &dimensions, vector< vector< vector< unsigned char > > > &brickImages )
{
  // Save image as static image or zoomable image
  if( vips_foreign_find_save( outputImage.c_str() ) != NULL )
  {
    int threads = numberOfCPUS();

    ProgressBar *topLevel = new ProgressBar(ceil((double)brickData.size()/(10*threads)), "Generating image");

    future< void > ret[threads];

    unsigned char * tileData = ( unsigned char * )calloc (outputHeight*outputWidth*3,sizeof(unsigned char));

    for( int k = 0; k < threads; ++k )
    {
      int start = k*brickData.size()/threads;
      int end = (k+1)*brickData.size()/threads;

      ret[k] = async( launch::async, &buildTopLevelThread, outputImage, start, end, tileData,  outputWidth, outputHeight, ref(brickData), ref(dimensions), ref(brickImages), outputWidth, outputHeight, topLevel );
    }

    // Wait for threads to finish
    for( int k = 0; k < threads; ++k )
    {
      ret[k].get();
    }

    topLevel->Finish();

    cout << "Saving " << outputImage << endl;

    VImage::new_from_memory( tileData, outputWidth*outputHeight*3, outputWidth, outputHeight, 3, VIPS_FORMAT_UCHAR ).vipssave((char *)outputImage.c_str());
  }
  else
  {
    int level = (int)ceil(log2( max(outputWidth,outputHeight) ) );

    if( outputImage.back() == '/' ) outputImage = outputImage.substr(0, outputImage.size()-1);

    g_mkdir(string(outputImage).append("_files/").c_str(), 0777);
    g_mkdir(string(outputImage).append("_files/").append(to_string(level)).c_str(), 0777);

    int threads = numberOfCPUS();

    ProgressBar *topLevel = new ProgressBar(ceil((double)outputWidth/256.0)*ceil((double)outputHeight/((double)threads*256.0)), "Building top level");

    future< void > ret[threads];

    for( int k = 0; k < threads; ++k )
    {
      int start = k*outputHeight/threads;
      int end = (k+1)*outputHeight/threads;

      start = (start / 256 ) * 256;
      if( k+1 < threads ) end = (end / 256 ) * 256;

      ret[k] = async( launch::async, &buildTopLevel, string(outputImage).append("_files/"+to_string(level)+"/"), start, end, outputWidth, outputHeight, ref(brickData), ref(dimensions), ref(brickImages), 256, 256, topLevel );
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

void generateOutput( vector< vector< int > > mosaic, string outputImage, bool randomize )
{
  vector< vector< vector< unsigned char > > > brickImages(25);
  vector< pair< int, int > > dimensions(25);

  size_t lastindex = outputImage.find_last_of("."); 
  string outputFile = outputImage.substr(0, lastindex).append(".ldr"); 

  ofstream output( outputFile );

  int width = mosaic[0].size();
  int height = mosaic.size();

  int width2 = (int)ceil(width*round(77*si));
  int height2 = (int)ceil(height*round(77*si));

  vector< vector< int > > brickData;
  vector< vector< bool > > brickUsed(25,vector< bool >(68,false));

  for( int brick = 0; brick < numBricks; ++brick )
  {
    int w = brickSizes[brick][0];
    int h = brickSizes[brick][1];

    vector< int > indices( (width-w+1)*(height-h+1) );
    iota( indices.begin(), indices.end(), 0 );

    // Shuffle the points so that patterens do not form
    if( randomize ) shuffle( indices.begin(), indices.end(), default_random_engine(0));

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

        brickData.push_back({brick,colorIndex,bi,bj});

        brickUsed[brick][colorIndex] = true;

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

        brickData.push_back({24,colorIndex,bi,bj});

        brickUsed[24][colorIndex] = true;
      }
    }
  }

  output.close();

  generateImages( brickImages, dimensions, brickUsed );

  buildImage( outputImage, width2, height2, brickData, dimensions, brickImages );
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

void legoMosaicThread( int threadIdx, int numThreads, int width, int height, float colors[68][3], vector< vector< int > > &mosaic, vector< vector< vector< float > > > &floatData, vector< int > &colorsToUse, bool dither, condition_variable *cv, ProgressBar *buildingMosaic )
{
  mutex m;
  unique_lock<mutex> lk(m);

  for( int i = threadIdx; i < height; i += numThreads )
  {
    for( int j = 0; j < width; ++j )
    {
      while( dither && i > 0 && j < width-1 && mosaic[i-1][j+1] < 0 )
      {
        cv[threadIdx].wait(lk);
      }


      if( dither )
      {
        for( int i1 = 0; i1 < 3; ++i1 )
        {
          if( j > 0 )
          {
            floatData[i][j][i1] += floatData[i][j-1][i1] * 7.0/16.0;
          }

          if( i > 0 )
          {
            if( j > 0 )
            {
              floatData[i][j][i1] += floatData[i-1][j-1][i1] * 1.0/16.0;
            }

            floatData[i][j][i1] += floatData[i-1][j][i1] * 5.0/16.0;

            if( j < width-1 )
            {
              floatData[i][j][i1] += floatData[i-1][j+1][i1] * 3.0/16.0;
            }
          }

          fit(floatData[i][j][i1],i1);
        }
      }

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
      
        for( int i1 = 0; i1 < 3; ++i1 )
        {
          floatData[i][j][i1] = error[i1];
        }
      }

      mosaic[i][j] = color_codes[bestColor];

      cv[(threadIdx+1)%numThreads].notify_one();
    }

    if( threadIdx == 0 ) buildingMosaic->Increment();
  }
}

void generateLegoMosaic( string inputImage, string outputImage, int numAcross, bool dither, bool randomize, int dilateDistance, double gamma, vector< int > &colorsToUse )
{
  // Load input image
  VImage image = VImage::vipsload( (char *)inputImage.c_str() ).autorot();

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
  unsigned char * inputData2 = ( unsigned char * )image.data();
  unsigned char * inputData3 = new unsigned char[width*height*3];

  dilate( width, height, inputData2, inputData3, dilateDistance );

  VImage image2 = VImage::new_from_memory( inputData3, width*height*3, width, height, 3, VIPS_FORMAT_UCHAR ).thumbnail_image(numAcross);//.flip(VIPS_DIRECTION_HORIZONTAL);

  width = image2.width();
  height = image2.height();

  unsigned char * inputData = ( unsigned char * )image2.data();

  float * c2 = new float[3*width*height];

  vector< vector< vector< float > > > floatData(height,vector< vector< float > >(width,vector< float >(3,0)));

  int threads = numberOfCPUS();

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

  ProgressBar *buildingMosaic = new ProgressBar(height/threads+1, "Building Mosaic");

  condition_variable *cv = new condition_variable[threads];

  for( int k = 0; k < threads; ++k )
  {
    ret2[k] = async( launch::async, &legoMosaicThread, k, threads, width, height, colors, ref(mosaic), ref(floatData), ref(colorsToUse), dither, cv, buildingMosaic );
  }

  // Wait for threads to finish
  for( int k = 0; k < threads; ++k )
  {
    ret2[k].get();
  }

  buildingMosaic->Finish();

  generateOutput( mosaic, outputImage, randomize );
}

void generateImages2( vector< vector< vector< unsigned char > > > &brickImages, vector< pair< int, int > > &dimensions, vector< vector< bool > > &brickUsed )
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

      if( brickUsed[i][j] && brickUsed[i+18][j] )
      {
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
      }
      else if( brickUsed[i][j] )
      {
        void * data = (void *) malloc (sizeof(unsigned char)*len);

        shapeData.read( (char *)data, len*sizeof(char) );

        VipsBlob *blob = vips_blob_new( NULL, data, len );

        VImage image = VImage::pngload_buffer(blob).resize(si);

        unsigned char * c = ( unsigned char * )image.data();
        dimensions[i] = pair< int, int >( image.width(), image.height() );
        brickImages[i].push_back( vector< unsigned char >(c, c + 3*image.width()*image.height() ) );

        brickImages[i+18].resize(brickImages[i+18].size()+1);
      }
      else if( brickUsed[i+18][j] )
      {
        void * data = (void *) malloc (sizeof(unsigned char)*len);

        shapeData.read( (char *)data, len*sizeof(char) );

        VipsBlob *blob = vips_blob_new( NULL, data, len );

        VImage image = VImage::pngload_buffer(blob).resize(si);

        unsigned char * c = ( unsigned char * )image.rot(VIPS_ANGLE_D270).data();
        dimensions[i+18] = pair< int, int >( image.height(), image.width() );
        brickImages[i+18].push_back( vector< unsigned char >(c, c + 3*image.width()*image.height() ) );

        brickImages[i].resize(brickImages[i].size()+1);
      }
      else
      {
        shapeData.seekg( len, ios::cur );

        brickImages[i].resize(brickImages[i].size()+1);
        brickImages[i+18].resize(brickImages[i+18].size()+1);
      }

      generatingImages->Increment();
    }
  }

  generatingImages->Finish();
}

void generateOutput2( vector< vector< int > > &mosaic1, vector< vector< int > > &mosaic2, string outputImage, bool randomize )
{
  vector< vector< vector< unsigned char > > > brickImages(36);
  vector< pair< int, int > > dimensions(36);

  size_t lastindex = outputImage.find_last_of("."); 
  string outputFile = outputImage.substr(0, lastindex).append(".ldr"); 

  ofstream output( outputFile );

  int width = mosaic1[0].size();
  int height = mosaic1.size();

  int width2 = (int)ceil(width/5*round(77*si));
  int height2 = (int)ceil(mosaic2.size()*round(77*si));

  vector< vector< int > > brickData;
  vector< vector< bool > > brickUsed(36,vector< bool >(68,false));

  for( int brick = 0; brick < 16; ++brick )
  {
    int w = brickSizes2[brick][0]*5;
    int h = brickSizes2[brick][1];

    vector< int > indices( (width-w+1)*(height-h+1) );
    iota( indices.begin(), indices.end(), 0 );

    // Shuffle the points so that patterens do not form
    if( randomize ) shuffle( indices.begin(), indices.end(), default_random_engine(0));

    for( int k = 0; k < indices.size(); ++k )
    {
      int i = indices[k]/(width-w+1);
      int j = indices[k]%(width-w+1);

      int squareColor = mosaic1[i][j];

      if( squareColor < 0 ) continue;

      bool badColor = false;

      for( int x = j; x < j + w && !badColor; x+=5 )
      {
        if( brick < 11 && ( i < 1 || mosaic1[i-1][x] == -1 ) )
        {
          badColor = true; 
        }
      }

      for( int y = i; y < i + h && !badColor; ++y )
      {
        for( int x = j; x < j + w && !badColor; x+=5 )
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
          for( int x = j; x < j + w; x+=5 )
          {
            mosaic1[y][x] = -2;
          }
        }

        int bi = floor(30.8*si*i);
        int bj = floor(15.4*si*j);

        int colorIndex = 0;

        for( int l = 0; l < 68; ++l )
        {
          if( color_codes[l] == squareColor )
          {
            colorIndex = l;
            break;
          }
        }

        brickData.push_back({brick,colorIndex,bi,bj});
        brickUsed[brick][colorIndex] = true;

        output << " 1 " << squareColor << " " << j*4+w*2 << " " << i*8 << " 0" << " 1 0 0 0 1 0 0 0 1 " << brickNames3[brick] << endl;
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
          int bi = floor(30.8*si*i);
          int bj = floor(15.4*si*j);

          int colorIndex = 0;

          for( int l = 0; l < 68; ++l )
          {
            if( color_codes[l] == mosaic1[i][j] )
            {
              colorIndex = l;
              break;
            }
          }

          brickData.push_back({16,colorIndex,bi,bj});
          brickUsed[16][colorIndex] = true;

          output << " 1 " << mosaic1[i][j] << " " << j*4+10 << " " << i*8 << " 0" << " 1 0 0 0 1 0 0 0 1 30039.dat\n";
        }
        else
        {
          int bi = floor(30.8*si*i);
          int bj = floor(15.4*si*j);

          int colorIndex = 0;

          for( int l = 0; l < 68; ++l )
          {
            if( color_codes[l] == mosaic1[i][j] )
            {
              colorIndex = l;
              break;
            }
          }

          brickData.push_back({17,colorIndex,bi,bj});
          brickUsed[17][colorIndex] = true;

          output << " 1 " << mosaic1[i][j] << " " << j*4+10 << " " << i*8 << " 0" << " 1 0 0 0 1 0 0 0 1 3024.dat\n";
        }
      }
    }
  }

  width = mosaic2[0].size();
  height = mosaic2.size();

  for( int brick = 0; brick < 16; ++brick )
  {
    int w = brickSizes2[brick][1]*2;

    int h = brickSizes2[brick][0];

//    cout << w << " " << h << endl;

    if( width < w || height < h ) continue;

    vector< int > indices( (width-w+1)*(height-h+1) );

    iota( indices.begin(), indices.end(), 0 );

    // Shuffle the points so that patterens do not form
    if( randomize ) shuffle( indices.begin(), indices.end(), default_random_engine(0));

    for( int k = 0; k < indices.size(); ++k )
    {
      int i = indices[k]/(width-w+1);
      int j = indices[k]%(width-w+1);

      int squareColor = mosaic2[i][j];

      if( squareColor < 0 ) continue;

      bool badColor = false;

      for( int y = i; y < i + h && !badColor; ++y )
      {
        if( brick < 11 && ( j < 2 || mosaic2[y][j-2] == -1 ) )
        {
          badColor = true; 
        }
      }

      if( badColor ) continue;

//      cout << "wh " << w << " " << h << endl;

      for( int y = i; y < i + h && !badColor; ++y )
      {
        for( int x = j; x < j + w && !badColor; x+=2 )
        {
//          cout << w << " " << h << " " << " | " << j << " " << i << " | " << x << " " << y << " | " << mosaic2[y][x] << " " << squareColor << endl;
          if( mosaic2[y][x] != squareColor )
          {
            badColor = true;
          }
        }
      }
//      cout << "here\n";

      if( !badColor )
      {
        for( int y = i; y < i + h; ++y )
        {
          for( int x = j; x < j + w; x+=2 )
          {
            mosaic2[y][x] = -2;
          }
        }

        int bi = floor(77.0*si*i);
        int bj = floor(15.4*si*j);

        int colorIndex = 0;

        for( int l = 0; l < 68; ++l )
        {
          if( color_codes[l] == squareColor )
          {
            colorIndex = l;
            break;
          }
        }

        brickData.push_back({brick+18,colorIndex,bi,bj});
        brickUsed[brick+18][colorIndex] = true;

        output << " 1 " << squareColor << " " << j*4 << " " << i*20+h*10 << " 0" << " 0 1 0 -1 0 0 0 0 1 " << brickNames3[brick] << endl;
      }
    }
  }

  for( int i = 0; i < height; ++i )
  {
    for( int j = 0; j < width; ++j )
    {
      if( mosaic2[i][j] >= 0 )
      {
        if( j < 1 || mosaic2[i][j-2] == -1 )
        {
          int bi = floor(77.0*si*i);
          int bj = floor(15.4*si*j);

          int colorIndex = 0;

          for( int l = 0; l < 68; ++l )
          {
            if( color_codes[l] == mosaic2[i][j] )
            {
              colorIndex = l;
              break;
            }
          }

          brickData.push_back({34,colorIndex,bi,bj});
          brickUsed[34][colorIndex] = true;

          output << " 1 " << mosaic2[i][j] << " " << j*4 << " " << i*20+10 << " 0" << " 0 1 0 -1 0 0 0 0 1 30039.dat\n";
        }
        else
        {
          int bi = floor(77.0*si*i);
          int bj = floor(15.4*si*j);

          int colorIndex = 0;

          for( int l = 0; l < 68; ++l )
          {
            if( color_codes[l] == mosaic2[i][j] )
            {
              colorIndex = l;
              break;
            }
          }

          brickData.push_back({35,colorIndex,bi,bj});
          brickUsed[35][colorIndex] = true;

          output << " 1 " << mosaic2[i][j] << " " << j*4 << " " << i*20+10 << " 0" << " 0 1 0 -1 0 0 0 0 1 3024.dat\n";
        }
      }
    }
  }

  output.close();

  generateImages2( brickImages, dimensions, brickUsed );

  buildImage( outputImage, width2, height2, brickData, dimensions, brickImages );
}

void reorder( int width, int height, vector< vector< int > > &mosaic1, vector< vector< int > > &mosaic2 )
{
  int start = 0;
  int end = start;

  while( start < width )
  {
    int newStart = -1;
    int startColor = -1;

    for( int j = start; j < width+1; ++j )
    {
      if( j == width )
      {
        end = width;
        newStart = width;
        break;
      }

      for( int i = 0; i < 2; ++i )
      {
        if( mosaic2[i][j] < 0 ) break;
        if( startColor < 0 ) startColor = mosaic2[i][j];
        if( mosaic2[i][j] != startColor )
        {
          end = j;
          newStart = j;
          if( newStart == start ) newStart += 2;
          break;
        }
      }
      if( newStart != -1 ) break;

      for( int i = 0; i < 5; ++i )
      {
        if( mosaic1[i][j] < 0 ) break;
        if( startColor < 0 ) startColor = mosaic1[i][j];
        if( mosaic1[i][j] != startColor )
        {
          end = j;
          newStart = j;
          if( newStart == start ) newStart += 5;
          break;
        }
      }
      if( newStart != -1 ) break;
    }

    int length = end - start;

    int firstEnd;

    for( int l = length, fe = 0; l >= 0; l -= 5, fe += 5 )
    {
      if( l == 0 )
      {
        firstEnd = end;
      }
      else if( l < 5 && l%2 == 0 )
      {
        firstEnd = start + fe;
      }
      else if( l < 5 && l%2 == 1 )
      {
        firstEnd = start + fe - 5;
      }
    }

    for( int j = start; j < end; ++j )
    {
      for( int i = 0; i < 2; ++i )
      {
        mosaic2[i][j] = -1;
      }

      for( int i = 0; i < 5; ++i )
      {
        mosaic1[i][j] = -1;
      }
    }

    for( int j = start; j < firstEnd; j+=5 )
    {
      for( int i = 0; i < 5; ++i )
      {
        mosaic1[i][j] = startColor;
      }
    }

    for( int j = firstEnd; j < end; j+=2 )
    {
      for( int i = 0; i < 2; ++i )
      {
        mosaic2[i][j] = startColor;
      }
    }

    start = newStart;
  }
}

void bestRow( int i1, int width, int height, vector< vector< vector< float > > > &floatData, vector< vector< int > > &mosaic1, vector< vector< int > > &mosaic2, bool dither, float colors[][3], vector< int > colorsToUse )
{  
  vector< int > allIndices( width+1, -1 );
  vector< float > allDifferences( width+1, -1 );

  vector< vector< int > > allMosaic1( width+1, vector< int >( 5, -1 ) ), allMosaic2( width+1, vector< int >( 5, -1 ) );

  vector< vector< vector< vector< float > > > > allFloatData(width+1, vector< vector< vector< float > > >( 15, vector< vector< float > >( 15 ) ) );

  allDifferences[0] = 0;

  for( int i = 0; i < 15; ++i )
  {
    if( i1+i >= height ) break;
    for( int j = 0; j < 10; ++j )
    {
      allFloatData[0][i][j+5] = floatData[i1+i][j];
    }
  }

  for( int k = 1; k < width-1; ++k )
  {
    for( int i = 0; i < 15; ++i )
    {
      if( i1+i >= height ) break;
      for( int j = 5; j < 10; ++j )
      {
        if( k+j >= width ) break;

        allFloatData[k][i][j+5] = floatData[i1+i][k+j];
      }
    }
  }

  for( int j1 = 0; j1 < width-1; ++j1 )
  {
    if( allDifferences[j1] < 0 ) continue;
  
    if( j1 == 1 || j1 == 3 ) continue;

    float tempDifference = allDifferences[j1];
    vector< vector< vector< float > > > tempFloatData = allFloatData[j1];

    vector< int > tempMosaic(5);

    if( j1 < width-4 )
    {
      for( int i = 0; i < 5; ++i )
      {
        int bestColor = -1;

        float averageColor[3] = {0,0,0};

        float bestDifference = -1;

        if( dither )
        {
          for( int y = i*2; y < i*2+2; ++y )
          {
            for( int x = 0; x < 5; ++x )
            {
              for( int i10 = 0; i10 < 3; ++i10 )
              {
                averageColor[i10] += tempFloatData[y][x+5][i10];
              }
            }
          }

          for( int i10 = 0; i10 < 3; ++i10 )
          {
            averageColor[i10] /= 10.0;
          }
        }

        for( int k1 = 0; k1 < colorsToUse.size(); ++k1 )
        {
          int k = colorsToUse[k1];

          int r = 0;
          int g = 0;
          int b = 0;

          float difference = 0;

          for( int y = i*2; y < i*2+2; ++y )
          {
            for( int x = 0; x < 5; ++x )
            {
              float tempD = 0;
              for( int i10 = 0; i10 < 3; ++i10 )
              {
                float diff = colors[k][i10] - tempFloatData[y][x+5][i10];
                tempD += diff*diff;
              }
              difference += sqrt(tempD);
            }
          }

          if( bestDifference < 0 || difference < bestDifference )
          {
            bestDifference = difference;
            bestColor = k;
          }
        }

        tempDifference += bestDifference;
        tempMosaic[i] = bestColor;

        if( dither )
        {
          float error[3];

          for( int i10 = 0; i10 < 3; ++i10 )
          {
            error[i10] = averageColor[i10]-colors[bestColor][i10];
          }

          if( j1 < width - 10 )
          {
            for( int y = i*2; y < i*2+2; ++y )
            {
              for( int x = 0; x < 5; ++x )
              {
                int x1 = x+10;
                int y1 = y;

                for( int i10 = 0; i10 < 3; ++i10 )
                {
                  tempFloatData[y1][x1][i10] += error[i10] * 7.0 / 16.0;
                  fit(tempFloatData[y1][x1][i10],i10);
                }
              }
            }
          }

          if( i1+i*2 < height - 2 )
          {
            if( j1 > 4 )
            {
              for( int y = i*2; y < i*2+2; ++y )
              {
                for( int x = 0; x < 5; ++x )
                {
                  int x1 = x;
                  int y1 = y+2;

                  for( int i10 = 0; i10 < 3; ++i10 )
                  {
                    tempFloatData[y1][x1][i10] += error[i10] * 3.0 / 16.0;
                    fit(tempFloatData[y1][x1][i10],i10);
                  }
                }
              }
            }

            {
              for( int y = i*2; y < i*2+2; ++y )
              {
                for( int x = 0; x < 5; ++x )
                {
                  int x1 = x+5;
                  int y1 = y+2;

                  for( int i10 = 0; i10 < 3; ++i10 )
                  {
                    tempFloatData[y1][x1][i10] += error[i10] * 5.0 / 16.0;
                    fit(tempFloatData[y1][x1][i10],i10);
                  }
                }
              }
            }

            if( j1 < width - 10 )
            {
              for( int y = i*2; y < i*2+2; ++y )
              {
                for( int x = 0; x < 5; ++x )
                {
                  int x1 = x+10;
                  int y1 = y+2;

                  for( int i10 = 0; i10 < 3; ++i10 )
                  {
                    tempFloatData[y1][x1][i10] += error[i10] * 1.0 / 16.0;
                    fit(tempFloatData[y1][x1][i10],i10);
                  }
                }
              }
            }
          }
        }
      }

      if( allDifferences[j1+5] < 0 || tempDifference < allDifferences[j1+5] )
      {
        allDifferences[j1+5] = tempDifference;

        for( int i = 0; i < 5; ++i )
        {
          allMosaic1[j1+5][i] = color_codes[tempMosaic[i]];
        }

        allIndices[j1+5] = j1;

        for( int i = 0; i < 15; ++i )
        {
          for( int j = 0; j < 10; ++j )
          {
            allFloatData[j1+5][i][j] = tempFloatData[i][j+5];
          }
        }

        if( dither )
        {
          tempFloatData = allFloatData[j1];
        }
      }

      tempDifference = allDifferences[j1];
    }

    for( int i = 0; i < 2; ++i )
    {
      int bestColor = -1;

      float averageColor[3] = {0,0,0};

      double bestDifference = -1;

      if( dither )
      {
        for( int y = i*5; y < i*5+5; ++y )
        {
          for( int x = 0; x < 2; ++x )
          {
            for( int i10 = 0; i10 < 3; ++i10 )
            {
              averageColor[i10] += tempFloatData[y][x+5][i10];
            }
          }
        }

        for( int i10 = 0; i10 < 3; ++i10 )
        {
          averageColor[i10] /= 10.0;
        }
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
          for( int x = 0; x < 2; ++x )
          {
            float tempD = 0;
            for( int i10 = 0; i10 < 3; ++i10 )
            {
              float diff = colors[k][i10] - tempFloatData[y][x+5][i10];
              tempD += diff*diff;
            }
            difference += sqrt(tempD);
          }
        }

        if( bestDifference < 0 || difference < bestDifference )
        {
          bestDifference = difference;
          bestColor = k;
        }
      }

      tempDifference += bestDifference;
      tempMosaic[i] = bestColor;

      if( dither )
      {
        float error[3];

        for( int i10 = 0; i10 < 3; ++i10 )
        {
          error[i10] = averageColor[i10]-colors[bestColor][i10];
        }

        if( j1 < width - 4 )
        {
          for( int y = i*5; y < i*5+5; ++y )
          {
            for( int x = 0; x < 2; ++x )
            {
              int x1 = x+7;
              int y1 = y;

              for( int i10 = 0; i10 < 3; ++i10 )
              {
                tempFloatData[y1][x1][i10] += error[i10] * 7.0 / 16.0;
                fit(tempFloatData[y1][x1][i10],i10);
              }
            }
          }
        }

        if( i1+i*5 < height - 5 )
        {
          if( j1 > 1 )
          {
            for( int y = i*5; y < i*5+5; ++y )
            {
              for( int x = 0; x < 2; ++x )
              {
                int x1 = x+3;
                int y1 = y+5;

                for( int i10 = 0; i10 < 3; ++i10 )
                {
                  tempFloatData[y1][x1][i10] += error[i10] * 3.0 / 16.0;
                  fit(tempFloatData[y1][x1][i10],i10);
                }
              }
            }
          }

          {
            for( int y = i*5; y < i*5+5; ++y )
            {
              for( int x = 0; x < 2; ++x )
              {
                int x1 = x+5;
                int y1 = y+5;

                for( int i10 = 0; i10 < 3; ++i10 )
                {
                  tempFloatData[y1][x1][i10] += error[i10] * 5.0 / 16.0;
                  fit(tempFloatData[y1][x1][i10],i10);
                }
              }
            }
          }

          if( j1 < width - 4 )
          {
            for( int y = i*5; y < i*5+5; ++y )
            {
              for( int x = 0; x < 2; ++x )
              {
                int x1 = x+7;
                int y1 = y+5;

                for( int i10 = 0; i10 < 3; ++i10 )
                {
                  tempFloatData[y1][x1][i10] += error[i10] * 1.0 / 16.0;
                  fit(tempFloatData[y1][x1][i10],i10);
                }
              }
            }
          }
        }
      }
    }

    if( allDifferences[j1+2] < 0 || tempDifference < allDifferences[j1+2] )
    {
      allDifferences[j1+2] = tempDifference;

      for( int i = 0; i < 2; ++i )
      {
        allMosaic2[j1+2][i] = color_codes[tempMosaic[i]];
      }

      allIndices[j1+2] = j1;

      for( int i = 0; i < 15; ++i )
      {
        for( int j = 0; j < 10; ++j )
        {
          allFloatData[j1+2][i][j] = tempFloatData[i][j+2];
        }
      }
    }
  }

  vector< vector< int > > tempMosaic1( 5, vector< int >( width, -1 ) ), tempMosaic2( 2, vector< int >( width, -1 ) );

  vector< int > ind;

  for( int k = width; k > 0; k = allIndices[k] )
  {
    ind.push_back(k);
  }
  ind.push_back(0);

  for( int k1 = ind.size()-1; k1 > 0; --k1 )
  {
    int ik1 = ind[k1];
    int ik2 = ind[k1-1];

    int length = ik2-ik1;
    if( length == 2 )
    {
      for( int i = 0; i < 2; ++i )
      {
        tempMosaic2[i][ik1] = allMosaic2[ik2][i];
      }
    }
    else
    {
      for( int i = 0; i < 5; ++i )
      {
        tempMosaic1[i][ik1] = allMosaic1[ik2][i];
      }
    }

    if( dither )
    {
      for( int i = 10; i < 15; ++i )
      {
        if( i1+i >= height ) break;
        for( int j = 0; j < 5; ++j )
        {
          if( ik2+j-5 < 0 ) continue;
          if( ik2+j-5 >= width ) break;
          floatData[i1+i][ik2+j-5] = allFloatData[ik2][i][j];
        }
      }
    }
  }

  reorder( width, height, tempMosaic1, tempMosaic2 );

  for( int i = 0; i < 5; ++i )
  {
    mosaic1[i1/2+i] = tempMosaic1[i];
  }

  for( int i = 0; i < 2; ++i )
  {
    mosaic2[i1/5+i] = tempMosaic2[i];
  }
}

void legoMosaicThread2( int threadIdx, int numThreads, int width, int height, float colors[68][3], vector< vector< int > > &mosaic1, vector< vector< int > > &mosaic2, vector< vector< vector< float > > > &floatData, vector< vector< vector< float > > > &floatData2, vector< int > &colorsToUse, bool dither, ProgressBar *buildingMosaic )
{
  for( int i = threadIdx*10; i < height; i += numThreads*10 )
  {
    bestRow( i, width, height, floatData, mosaic1, mosaic2, dither, colors, colorsToUse );

    if( threadIdx == 0 ) buildingMosaic->Increment();
  }
}

void generateLegoMosaic2( string inputImage, string outputFile, int numAcross, bool dither, bool randomize, float dilateDistance, double gamma, vector< int > &colorsToUse )
{
  // Load input image
  VImage image = VImage::vipsload( (char *)inputImage.c_str() ).autorot();

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
  unsigned char * inputData2 = ( unsigned char * )image.data();
  unsigned char * inputData3 = new unsigned char[width*height*3];

  dilate( width, height, inputData2, inputData3, dilateDistance );

  VImage image2 = VImage::new_from_memory( inputData3, width*height*3, width, height, 3, VIPS_FORMAT_UCHAR ).thumbnail_image((numAcross-numAcross%2)*5);//.flip(VIPS_DIRECTION_HORIZONTAL);

  width = image2.width();
  height = image2.height();

  width -= width%10;
  height -= height%10;

  image2 = image2.extract_area(0,0,width,height);

  unsigned char * inputData = ( unsigned char * )image2.data();

  float * c2 = new float[3*width*height];

  vector< vector< vector< float > > > floatData(height,vector< vector< float > >(width,vector< float >(3,0))), floatData2;

  int threads = numberOfCPUS();

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

  if( dither ) threads = 1;

  vector< vector< int > > mosaic1( height/2, vector< int > ( width, -1 ) );
  vector< vector< int > > mosaic2( height/5, vector< int > ( width, -1 ) );

  ProgressBar *buildingMosaic = new ProgressBar(ceil((double)height/(10*threads)), "Building Mosaic");
  
  for( int k = 0; k < threads; ++k )
  {
    ret2[k] = async( launch::async, &legoMosaicThread2, k, threads, width, height, colors, ref(mosaic1), ref(mosaic2), ref(floatData), ref(floatData2), ref(colorsToUse), dither, buildingMosaic );
  }

  // Wait for threads to finish
  for( int k = 0; k < threads; ++k )
  {
    ret2[k].get();
  }

  buildingMosaic->Finish();

  generateOutput2( mosaic1, mosaic2, outputFile, randomize );
}

void testGrid( int j1, int i1, int width, int height, float tempDifference, vector< vector< vector< float > > > &floatData, vector< string > &allMosaic, vector< float > &allDifferences, vector< int > &allIndices, float colors[][3], vector< int > colorsToUse, vector< vector< int > > &grid, vector< vector< int > > &offsets, bool doFlipHorizontal, bool doFlipVertical, bool doRotate, bool doReplace )
{
  int tileWidth = grid[0].size();
  int tileHeight = grid.size();

  if( j1 > width-tileWidth ) return;

  stringstream output;

  set< int > blockSet;

  for( int y = 0; y < tileHeight; ++y )
  {
    for( int x = 0; x < tileWidth; ++x )
    {
      blockSet.insert( grid[y][x] );
    }
  }

  vector< int > blocks( blockSet.begin(), blockSet.end() );

  for( int b2 = 0; b2 < blocks.size(); ++b2 )
  {
    int b1 = blocks[b2];

    int bestColor = -1;

    float bestDifference = -1;

    for( int k1 = 0; k1 < colorsToUse.size(); ++k1 )
    {
      int k = colorsToUse[k1];

      int r = 0;
      int g = 0;
      int b = 0;

      float difference = 0;

      for( int y = 0; y < tileHeight; ++y )
      {
        for( int x = 0; x < tileWidth; ++x )
        {
          if( grid[y][x] != b1 ) continue;

          float tempD = 0;
          for( int i10 = 0; i10 < 3; ++i10 )
          {
            int x2 = doFlipHorizontal^doRotate ? j1+tileWidth-x-1 : j1+x;
            int y2 = doFlipVertical^doRotate ? i1+tileHeight-y-1 : i1+y;
            float diff = colors[k][i10] - floatData[y2][x2][i10];

            tempD += diff*diff;
          }
          difference += sqrt(tempD);
        }
      }

      if( bestDifference < 0 || difference < bestDifference )
      {
        bestDifference = difference;
        bestColor = k;
      }
    }

    tempDifference += bestDifference;

    int j2 = doFlipHorizontal^doRotate ? j1*4+tileWidth*4-offsets[b1][0] : j1*4+offsets[b1][0];
    int i2 = doFlipVertical^doRotate ? i1*4+tileHeight*4-offsets[b1][1] : i1*4+offsets[b1][1];

    output << " 1 " << color_codes[bestColor] << " " << j2 << " " << i2 << " 0";

    for( int o1 = 3; o1 < 12; ++o1 )
    {
      int i = (o1-3)/3;

      int m1 = offsets[b1][o1];

      if( doFlipHorizontal && i != 1 )
      {
        m1 *= -1;
      }

      if( doFlipVertical && i > 0 )
      {
        m1 *= -1;
      }

      if( doRotate && i < 2 )
      {
        m1 *= -1;
      }

      output << " " << m1;
    }

    output << " " << offsets[b1][12] << ".dat\n";
  }

  if( allDifferences[j1+tileWidth] < 0 || tempDifference < allDifferences[j1+tileWidth]+(doReplace?0.1:-0.1) )
  {
    allDifferences[j1+tileWidth] = tempDifference;

    allMosaic[j1+tileWidth] = output.str();

    allIndices[j1+tileWidth] = j1;
  }
}

void bestSmallRow( int i1, int width, int height, int numAllGrids, int numTallGrids, int numShortGrids, vector< vector< vector< float > > > &floatData, vector< vector< vector< string > > > &allMosaic2, vector< vector< float > > &allDifferences, vector< vector< vector< int > > > &allIndices2, float colors[][3], vector< int > colorsToUse, vector< vector< vector< int > > > &grid, vector< vector< int > > &offsets, vector< bool > &doFlipHorizontal, vector< bool > &doFlipVertical, vector< bool > &doRotate )
{
  vector< vector< vector< float > > > allDifferences2( 2, vector< vector< float > >( width+1, vector< float >( width+1, -1 ) ) );

  for( int j1 = 0; j1 < width-1; ++j1 )
  {
    for( int g1 = numTallGrids; g1 < numAllGrids; ++g1 )
    {
      for( int dfh = 0; dfh <= doFlipHorizontal[g1]; ++dfh )
      {
        for( int dfv = 0; dfv <= doFlipVertical[g1]; ++dfv )
        {
          for( int dr = 0; dr <= doRotate[g1]; ++dr )
          {
            for( int level = 0; level < 2; ++level )
            {
              testGrid( j1, i1+level*5, width, height, 0, floatData, allMosaic2[level][j1], allDifferences2[level][j1], allIndices2[level][j1], colors, colorsToUse, grid[g1], offsets, dfh, dfv, dr, g1 == numTallGrids );
            }
          }
        }
      }
    }
  }

  for( int start1 = width; start1 >= 0; --start1 )
  {
    for( int end = start1; end <= width; ++end )
    {
      for( int start2 = start1; start2 <= width; ++start2 )
      {
        for( int level = 0; level < 2; ++level )
        {
          float difference = allDifferences2[level][start1][start2];

          if( difference < 0 ) continue;
          if( allDifferences2[level][start2][end] < 0 ) continue;

          if( allDifferences2[level][start1][end] < 0 || difference + allDifferences2[level][start2][end] < allDifferences2[level][start1][end] + ((end-start2)>2?0.1:-0.1) )
          {
            allDifferences2[level][start1][end] = difference + allDifferences2[level][start2][end];
            allMosaic2[level][start1][end] = allMosaic2[level][start2][end];
            allIndices2[level][start1][end] = allIndices2[level][start2][end];
          }
        }
      }
    }
  }

  for( int start = 0; start < width; ++start )
  {
    for( int end = start+1; end <= width; ++end )
    {
      if( allDifferences2[0][start][end] < 0 ) continue;

      allDifferences[start][end] = allDifferences2[0][start][end] + allDifferences2[1][start][end];
    }
  };
}

void bestRow3( int i1, int width, int height, bool useSlopes, vector< vector< vector< float > > > &floatData, vector< string > &outputMosaic, float colors[][3], vector< int > colorsToUse )
{
  vector< int > allIndices( width+1, -1 );
  vector< float > allDifferences( width+1, -1 );

  vector< string > allMosaic( width+1 );

  vector< vector< float > > allDifferencesSmall( width+1, vector< float >( width+1, -1 ) );
  vector< vector< vector< string > > > allMosaicSmall( 2, vector< vector< string > >( width+1, vector< string >( width+1 ) ) );
  vector< vector< vector< int > > > allIndicesSmall( 2, vector< vector< int > >( width+1, vector< int >( width+1, -1 ) ) );

  allDifferences[0] = 0;

  int numTallGrids = 6;
  int numShortGrids = 6;
  int numAllGrids = numTallGrids + numShortGrids;

  vector< vector< vector< int > > > grid(numAllGrids);

  for( int i = 0; i < numTallGrids; ++i )
  {
    grid[i] = vector< vector< int > >(10);
  }

  for( int i = numTallGrids; i < numAllGrids; ++i )
  {
    grid[i] = vector< vector< int > >(5);
  }

  grid[0][0] = {0,0,0,0,0};
  grid[0][1] = {0,0,0,0,0};
  grid[0][2] = {1,1,1,1,1};
  grid[0][3] = {1,1,1,1,1};
  grid[0][4] = {2,2,2,2,2};
  grid[0][5] = {2,2,2,2,2};
  grid[0][6] = {3,3,3,3,3};
  grid[0][7] = {3,3,3,3,3};
  grid[0][8] = {4,4,4,4,4};
  grid[0][9] = {4,4,4,4,4};

  grid[1][0] = {5,6,6,6,6};
  grid[1][1] = {5,6,6,6,6};
  grid[1][2] = {5,5,6,6,6};
  grid[1][3] = {5,5,6,6,6};
  grid[1][4] = {5,5,6,6,6};
  grid[1][5] = {5,5,5,6,6};
  grid[1][6] = {5,5,5,6,6};
  grid[1][7] = {5,5,5,6,6};
  grid[1][8] = {5,5,5,5,6};
  grid[1][9] = {5,5,5,5,6};

  grid[2][0] = {13,13, 7, 7, 9, 9, 9, 9, 9,15,15};
  grid[2][1] = {13,13, 7, 7, 9, 9, 9, 9, 9,15,15};
  grid[2][2] = {13,13, 7, 7,10,10,10,10,10,15,15};
  grid[2][3] = {13,13, 7, 7,10,10,10,10,10,15,15};
  grid[2][4] = {13,13, 7, 7, 8, 8, 8, 8, 8,15,15};
  grid[2][5] = {14,14, 7, 7, 7, 7, 7, 8, 8,16,16};
  grid[2][6] = {14,14,11,11,11,11,11, 8, 8,16,16};
  grid[2][7] = {14,14,11,11,11,11,11, 8, 8,16,16};
  grid[2][8] = {14,14,12,12,12,12,12, 8, 8,16,16};
  grid[2][9] = {14,14,12,12,12,12,12, 8, 8,16,16};

  grid[3][0] = {17,17,19,19,21,21,21,21,21};
  grid[3][1] = {17,17,19,19,21,21,21,21,21};
  grid[3][2] = {17,17,19,19,22,22,22,22,22};
  grid[3][3] = {17,17,19,19,22,22,22,22,22};
  grid[3][4] = {17,17,19,19,18,18,18,18,18};
  grid[3][5] = {17,17,17,17,17,20,20,18,18};
  grid[3][6] = {23,23,23,23,23,20,20,18,18};
  grid[3][7] = {23,23,23,23,23,20,20,18,18};
  grid[3][8] = {24,24,24,24,24,20,20,18,18};
  grid[3][9] = {24,24,24,24,24,20,20,18,18};

  grid[4][0] = {27,27,27,27,27,28,28,28,28,28,29,29,29,29,29};
  grid[4][1] = {27,27,27,27,27,28,28,28,28,28,29,29,29,29,29};
  grid[4][2] = {25,25,25,25,25,25,25,25,25,25,25,25,25,25,25};
  grid[4][3] = {25,25,25,25,25,25,25,25,25,25,25,25,25,25,25};
  grid[4][4] = {25,25,25,25,25,25,26,26,26,25,25,25,25,25,25};
  grid[4][5] = {25,25,25,25,25,26,26,26,26,26,25,25,25,25,25};
  grid[4][6] = {25,25,25,25,25,26,26,26,26,26,25,25,25,25,25};
  grid[4][7] = {25,25,25,25,25,26,26,26,26,26,25,25,25,25,25};
  grid[4][8] = {30,30,30,30,30,31,31,31,31,31,32,32,32,32,32};
  grid[4][9] = {30,30,30,30,30,31,31,31,31,31,32,32,32,32,32};

  grid[5][0] = {41,41,41,41,41,33,35,35,35,35,39,39,37,37,37,37,37,37};
  grid[5][1] = {41,41,41,41,41,33,33,35,35,35,39,39,37,37,37,37,37,37};
  grid[5][2] = {42,42,42,42,42,33,33,33,35,35,39,39,37,45,45,45,45,45};
  grid[5][3] = {42,42,42,42,42,33,33,33,35,35,39,39,37,45,45,45,45,45};
  grid[5][4] = {43,43,43,43,43,33,33,33,33,35,39,39,37,46,46,46,46,46};
  grid[5][5] = {43,43,43,43,43,38,40,40,34,36,36,36,36,46,46,46,46,46};
  grid[5][6] = {44,44,44,44,44,38,40,40,34,34,36,36,36,47,47,47,47,47};
  grid[5][7] = {44,44,44,44,44,38,40,40,34,34,34,36,36,47,47,47,47,47};
  grid[5][8] = {38,38,38,38,38,38,40,40,34,34,34,36,36,48,48,48,48,48};
  grid[5][9] = {38,38,38,38,38,38,40,40,34,34,34,34,36,48,48,48,48,48};

  grid[6][0] = {49,49};
  grid[6][1] = {49,49};
  grid[6][2] = {49,49};
  grid[6][3] = {49,49};
  grid[6][4] = {49,49};

  grid[7][0] = {50,50,50,50,50,50,50,50,50,50};
  grid[7][1] = {50,50,50,50,50,50,50,50,51,51};
  grid[7][2] = {50,50,50,50,50,51,51,51,51,51};
  grid[7][3] = {50,50,51,51,51,51,51,51,51,51};
  grid[7][4] = {51,51,51,51,51,51,51,51,51,51};

  grid[8][0] = {52,52,52,52,53};
  grid[8][1] = {52,52,52,53,53};
  grid[8][2] = {52,52,52,53,53};
  grid[8][3] = {52,52,53,53,53};
  grid[8][4] = {52,53,53,53,53};

  grid[9][0] = {54,54,54,54,54};
  grid[9][1] = {54,54,54,54,55};
  grid[9][2] = {54,54,54,55,55};
  grid[9][3] = {54,55,55,55,55};
  grid[9][4] = {55,55,55,55,55};

  grid[10][0] = {56,56,56,56,56};
  grid[10][1] = {56,56,56,56,56};
  grid[10][2] = {56,56,56,56,56};
  grid[10][3] = {56,56,56,56,56};
  grid[10][4] = {57,57,57,57,57};

  grid[11][0] = {58,58,58,58,59};
  grid[11][1] = {58,58,58,58,59};
  grid[11][2] = {58,58,58,58,59};
  grid[11][3] = {58,58,58,58,59};
  grid[11][4] = {58,58,58,58,59};

  vector< vector< int > > offsets(60);

  offsets[0] = {10,0,0,1,0,0,0,1,0,0,0,1,3070};
  offsets[1] = {10,8,0,1,0,0,0,1,0,0,0,1,3070};
  offsets[2] = {10,16,0,1,0,0,0,1,0,0,0,1,3070};
  offsets[3] = {10,24,0,1,0,0,0,1,0,0,0,1,3070};
  offsets[4] = {10,32,0,1,0,0,0,1,0,0,0,1,3070};

  offsets[5] = {0,20,0,0,-1,0,0,0,1,-1,0,0,61409};
  offsets[6] = {20,20,0,0,1,0,0,0,-1,-1,0,0,61409};

  offsets[7] = {8,10,0,0,1,0,0,0,-1,-1,0,0,36841};
  offsets[8] = {36,30,0,0,-1,0,0,0,1,-1,0,0,36841};
  offsets[9] = {26,0,0,1,0,0,0,1,0,0,0,1,3070};
  offsets[10] = {26,8,0,1,0,0,0,1,0,0,0,1,3070};
  offsets[11] = {18,32,0,-1,0,0,0,-1,0,0,0,1,3070};
  offsets[12] = {18,40,0,-1,0,0,0,-1,0,0,0,1,3070};
  offsets[13] = {0,10,0,0,1,0,-1,0,0,0,0,1,3070};
  offsets[14] = {0,30,0,0,1,0,-1,0,0,0,0,1,3070};
  offsets[15] = {44,10,0,0,-1,0,1,0,0,0,0,1,3070};
  offsets[16] = {44,30,0,0,-1,0,1,0,0,0,0,1,3070};

  offsets[17] = {8,10,0,0,-1,0,0,0,-1,1,0,0,36840};
  offsets[18] = {28,30,0,0,1,0,0,0,1,1,0,0,36840};
  offsets[19] = {16,10,0,0,-1,0,1,0,0,0,0,1,3070};
  offsets[20] = {20,30,0,0,1,0,-1,0,0,0,0,1,3070};
  offsets[21] = {26,0,0,1,0,0,0,1,0,0,0,1,3070};
  offsets[22] = {26,8,0,1,0,0,0,1,0,0,0,1,3070};
  offsets[23] = {10,32,0,-1,0,0,0,-1,0,0,0,1,3070};
  offsets[24] = {10,40,0,-1,0,0,0,-1,0,0,0,1,3070};

  offsets[25] = {30,8,0,1,0,0,0,1,0,0,0,1,4490}; 
  offsets[26] = {30,32,0,0,0,-1,0,1,0,1,0,0,49307}; 
  offsets[27] = {10,0,0,1,0,0,0,1,0,0,0,1,3070}; 
  offsets[28] = {30,0,0,1,0,0,0,1,0,0,0,1,3070};
  offsets[29] = {50,0,0,1,0,0,0,1,0,0,0,1,3070};
  offsets[30] = {10,32,0,1,0,0,0,1,0,0,0,1,3070}; 
  offsets[31] = {30,32,0,1,0,0,0,1,0,0,0,1,3070};
  offsets[32] = {50,32,0,1,0,0,0,1,0,0,0,1,3070};

  offsets[33] = {20,10,0,0,-1,0,0,0,1,-1,0,0,54200};
  offsets[34] = {32,30,0,0,-1,0,0,0,1,-1,0,0,54200};
  offsets[35] = {40,10,0,0,1,0,0,0,-1,-1,0,0,54200};
  offsets[36] = {52,30,0,0,1,0,0,0,-1,-1,0,0,54200};
  offsets[37] = {62,8,0,0,0,1,0,-1,0,1,0,0,36840};
  offsets[38] = {10,32,0,0,0,-1,0,1,0,1,0,0,36840};
  offsets[39] = {40,10,0,0,1,0,-1,0,0,0,0,1,3024};
  offsets[40] = {32,30,0,0,-1,0,1,0,0,0,0,1,3024};
  offsets[41] = {10,0,0,1,0,0,0,1,0,0,0,1,3070};
  offsets[42] = {10,8,0,1,0,0,0,1,0,0,0,1,3070};
  offsets[43] = {10,16,0,1,0,0,0,1,0,0,0,1,3070};
  offsets[44] = {10,24,0,1,0,0,0,1,0,0,0,1,3070};
  offsets[45] = {62,16,0,-1,0,0,0,-1,0,0,0,1,3070};
  offsets[46] = {62,24,0,-1,0,0,0,-1,0,0,0,1,3070};
  offsets[47] = {62,32,0,-1,0,0,0,-1,0,0,0,1,3070};
  offsets[48] = {62,40,0,-1,0,0,0,-1,0,0,0,1,3070};

  offsets[49] = {0,10,0,0,1,0,-1,0,0,0,0,1,3070};

  offsets[50] = {20,0,0,0,0,-1,0,-1,0,-1,0,0,61409};
  offsets[51] = {20,20,0,0,0,1,0,1,0,-1,0,0,61409};

  offsets[52] = {0,10,0,0,-1,0,0,0,-1,1,0,0,54200};
  offsets[53] = {20,10,0,0,1,0,0,0,1,1,0,0,54200};

  offsets[54] = {10,0,0,0,0,-1,0,-1,0,-1,0,0,54200};
  offsets[55] = {10,20,0,0,0,1,0,1,0,-1,0,0,54200};

  offsets[56] = {10,0,0,1,0,0,0,-1,0,0,0,-1,54200};
  offsets[57] = {10,20,0,1,0,0,0,1,0,0,0,1,54200};

  offsets[58] = {0,10,0,0,-1,0,-1,0,0,0,0,-1,54200};
  offsets[59] = {20,10,0,0,1,0,-1,0,0,0,0,1,54200};

  vector< bool > doFlipHorizontal(numAllGrids,false);
  vector< bool > doFlipVertical(numAllGrids,false);
  vector< bool > doRotate(numAllGrids,false);

  doFlipHorizontal[0] = false;
  doFlipHorizontal[1] = true;
  doFlipHorizontal[2] = true;
  doFlipHorizontal[3] = true;
  doFlipHorizontal[4] = false;
  doFlipHorizontal[5] = true;
  doFlipHorizontal[6] = false;
  doFlipHorizontal[7] = true;
  doFlipHorizontal[8] = true;
  doFlipHorizontal[9] = true;
  doFlipHorizontal[10] = false;
  doFlipHorizontal[11] = false;

  doFlipVertical[0] = false;
  doFlipVertical[1] = false;
  doFlipVertical[2] = false;
  doFlipVertical[3] = false;
  doFlipVertical[4] = true;
  doFlipVertical[5] = false;
  doFlipVertical[6] = false;
  doFlipVertical[7] = false;
  doFlipVertical[8] = false;
  doFlipVertical[9] = false;
  doFlipVertical[10] = false;
  doFlipVertical[11] = false;

  doRotate[0] = false;
  doRotate[1] = false;
  doRotate[2] = false;
  doRotate[3] = false;
  doRotate[4] = false;
  doRotate[5] = false;
  doRotate[6] = false;
  doRotate[7] = false;
  doRotate[8] = false;
  doRotate[9] = false;
  doRotate[10] = true;
  doRotate[11] = true;

  if( !useSlopes )
  {
    numShortGrids = 1;
    numAllGrids = numTallGrids + numShortGrids;
  }

  bestSmallRow( i1, width, height, numAllGrids, numTallGrids, numShortGrids, floatData, allMosaicSmall, allDifferencesSmall, allIndicesSmall, colors, colorsToUse, grid, offsets, doFlipHorizontal, doFlipVertical, doRotate );

  for( int j1 = 0; j1 < width-1; ++j1 )
  {
    if( allDifferences[j1] < 0 ) continue;

    float tempDifference = allDifferences[j1];

    for( int g1 = 0; g1 < numTallGrids; ++g1 )
    {
      for( int dfh = 0; dfh <= doFlipHorizontal[g1]; ++dfh )
      {
        for( int dfv = 0; dfv <= doFlipVertical[g1]; ++dfv )
        {
          for( int dr = 0; dr <= doRotate[g1]; ++dr )
          {
            testGrid( j1, i1, width, height, tempDifference, floatData, allMosaic, allDifferences, allIndices, colors, colorsToUse, grid[g1], offsets, dfh, dfv, dr, g1 == 0 );
          }
        }
      }
    }

    for( int end = width; end <= width; ++end )
    {
      if( allDifferencesSmall[j1][end] < 0 ) continue;

      if( allDifferences[end] < 0 || tempDifference + allDifferencesSmall[j1][end] < allDifferences[end] )
      {
        allDifferences[end] = tempDifference + allDifferencesSmall[j1][end];

        allIndices[end] = j1;

        stringstream output;

        for( int level = 0; level < 2; ++level )
        {
          for( int k = end; k > j1; )
          {
            int nextK = allIndicesSmall[level][j1][k];
            output << allMosaicSmall[level][nextK][k];

            k = nextK;
          }  
        }

        allMosaic[end] = output.str();
      }
    }
  }

  stringstream output;

  for( int k = width; k > 0; k = allIndices[k] )
  {
    output << allMosaic[k];
  }

  outputMosaic[i1/10] = output.str();
}

void legoMosaicThread3( int threadIdx, int numThreads, int width, int height, bool useSlopes, float colors[68][3], vector< string > &allMosaic, vector< vector< vector< float > > > &floatData, vector< int > &colorsToUse, ProgressBar *buildingMosaic )
{
  for( int i = threadIdx*10; i < height; i += numThreads*10 )
  {
    bestRow3( i, width, height, useSlopes, floatData, allMosaic, colors, colorsToUse );

    if( threadIdx == 0 ) buildingMosaic->Increment();
  }
}

void generateLegoMosaic3( string inputImage, string outputFile, int numAcross, bool useSlopes, bool randomize, float dilateDistance, double gamma, vector< int > &colorsToUse )
{
  // Load input image
  VImage image = VImage::vipsload( (char *)inputImage.c_str() ).autorot();

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
  unsigned char * inputData2 = ( unsigned char * )image.data();
  unsigned char * inputData3 = new unsigned char[width*height*3];

  dilate( width, height, inputData2, inputData3, dilateDistance );

  VImage image2 = VImage::new_from_memory( inputData3, width*height*3, width, height, 3, VIPS_FORMAT_UCHAR ).thumbnail_image((numAcross-numAcross%2)*5);//.flip(VIPS_DIRECTION_HORIZONTAL);

  width = image2.width();
  height = image2.height();

  width -= width%10;
  height -= height%10;

  image2 = image2.extract_area(0,0,width,height);

  unsigned char * inputData = ( unsigned char * )image2.data();

  float * c2 = new float[3*width*height];

  vector< vector< vector< float > > > floatData(height,vector< vector< float > >(width,vector< float >(3,0))), floatData2;

  int threads = numberOfCPUS();

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

  vector< string > allMosaic(height/10);

  ProgressBar *buildingMosaic = new ProgressBar(ceil((double)height/(10*threads)), "Building Mosaic");
  
  for( int k = 0; k < threads; ++k )
  {
    ret2[k] = async( launch::async, &legoMosaicThread3, k, threads, width, height, useSlopes, colors, ref(allMosaic), ref(floatData), ref(colorsToUse), buildingMosaic );
  }

  // Wait for threads to finish
  for( int k = 0; k < threads; ++k )
  {
    ret2[k].get();
  }

  buildingMosaic->Finish();

  size_t lastindex = outputFile.find_last_of("."); 
  string outputFile2 = outputFile.substr(0, lastindex).append(".ldr"); 

  ofstream output( outputFile2 );

  for( int i = 0; i < allMosaic.size(); ++i )
  {
    output << allMosaic[i];
  }

  output.close();
}

void RunLegoMosaic( string inputName, string outputName, int numHorizontal, int tileSize, int style, bool dither, bool randomize, int dilateDistance, double gamma, string colorName )
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

  switch( style )
  {
    case 0:
      generateLegoMosaic( inputName, outputName, numHorizontal, dither, randomize, dilateDistance, gamma, colorsToUse );
      break;
    case 1:
      generateLegoMosaic2( inputName, outputName, numHorizontal, dither, randomize, dilateDistance, gamma, colorsToUse );
      break;
    case 2:
      generateLegoMosaic3( inputName, outputName, numHorizontal, false, randomize, dilateDistance, gamma, colorsToUse );
      break;
    case 3:
      generateLegoMosaic3( inputName, outputName, numHorizontal, true, randomize, dilateDistance, gamma, colorsToUse );
      break;
  }
}