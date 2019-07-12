#include <iostream>
#include <vector>
#include <cmath>
#include <climits>
#include <dirent.h>
#include <unistd.h>
#include <future>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <random>
#include <set>
#include <tuple>
#include <time.h>

#include <vips/vips8>
#include "progress_bar.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

using namespace vips;
using namespace std;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polyhedron_3<K>                     Polyhedron;
typedef K::Point_3                                Point;
typedef CGAL::Surface_mesh<Point>               Surface_mesh;
typedef K::Segment_3 Segment;

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

#define numBricks 24

const int brickSizes2[16][2] = 
{
  {8,3},
  {6,3},
  {4,3},
  {3,3},
  {8,1},
  {6,1},
  {2,3},
  {4,1},
  {3,1},
  {1,3},
  {2,1},
  {8,1},
  {6,1},
  {4,1},
  {3,1},
  {2,1}
};

const string brickNames3[16] =
{
  "3008.dat",
  "3009.dat",
  "3010.dat",
  "3622.dat",
  "3460.dat",
  "3666.dat",
  "3004.dat",
  "3710.dat",
  "3623.dat",
  "3005.dat",  
  "3023.dat",
  "4162.dat",
  "6636.dat",
  "2431.dat",
  "63864.dat",
  "3069b.dat"
};

const int brickSizes[numBricks][2] = 
{
  {8,4},
  {4,8},
  {4,6},
  {6,4},
  {4,4},
  {2,8},
  {8,2},
  {6,2},
  {2,6},
  {2,4},
  {4,2},
  {1,8},
  {8,1},
  {3,2},
  {2,3},
  {1,6},
  {6,1},
  {2,2},
  {4,1},
  {1,4},
  {1,3},
  {3,1},
  {2,1},
  {1,2}
};

const string brickNames2[14] =
{
  "3035.dat",
  "3032.dat",
  "3031.dat",
  "3034.dat",
  "3795.dat",
  "3020.dat",
  "3460.dat",
  "3021.dat",
  "3666.dat",
  "3022.dat",
  "3710.dat",
  "3623.dat",
  "3023.dat",
  "3024.dat"
};

const string brickNames[numBricks] =
{
  " 1 0 0 0 1 0 0 0 1 3035.dat",
  " 0 0 1 0 1 0 -1 0 0 3035.dat",
  " 0 0 1 0 1 0 -1 0 0 3032.dat",
  " 1 0 0 0 1 0 0 0 1 3032.dat",
  " 1 0 0 0 1 0 0 0 1 3031.dat",
  " 0 0 1 0 1 0 -1 0 0 3034.dat",
  " 1 0 0 0 1 0 0 0 1 3034.dat",
  " 1 0 0 0 1 0 0 0 1 3795.dat",
  " 0 0 1 0 1 0 -1 0 0 3795.dat",
  " 0 0 1 0 1 0 -1 0 0 3020.dat",
  " 1 0 0 0 1 0 0 0 1 3020.dat",
  " 0 0 1 0 1 0 -1 0 0 3460.dat",
  " 1 0 0 0 1 0 0 0 1 3460.dat",
  " 1 0 0 0 1 0 0 0 1 3021.dat",
  " 0 0 1 0 1 0 -1 0 0 3021.dat",
  " 0 0 1 0 1 0 -1 0 0 3666.dat",
  " 1 0 0 0 1 0 0 0 1 3666.dat",
  " 1 0 0 0 1 0 0 0 1 3022.dat",
  " 1 0 0 0 1 0 0 0 1 3710.dat",
  " 0 0 1 0 1 0 -1 0 0 3710.dat",
  " 0 0 1 0 1 0 -1 0 0 3623.dat",
  " 1 0 0 0 1 0 0 0 1 3623.dat",
  " 1 0 0 0 1 0 0 0 1 3023.dat",
  " 0 0 1 0 1 0 -1 0 0 3023.dat"
};

const int color_codes[68] =
{
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    13,
    14,
    15,
    17,
    18,
    19,
    20,
    22,
    23,
    25,
    26,
    27,
    28,
    29,
    68,
    69,
    70,
    71,
    72,
    73,
    74,
    77,
    78,
    85,
    86,
    89,
    92,
    100,
    110,
    112,
    115,
    118,
    120,
    125,
    151,
    191,
    212,
    216,
    226,
    232,
    272,
    288,
    308,
    313,
    320,
    321,
    323,
    335,
    351,
    373,
    378,
    379,
    450,
    462,
    484,
    503
};

const int color_values[68][3] =
{
  {11,10,9},
  {21,48,98},
  {11,74,19},
  {6,96,93},
  {142,4,15},
  {145,25,81},
  {50,28,18},
  {102,102,95},
  {50,51,46},
  {78,141,155},
  {15,92,17},
  {28,113,123},
  {180,122,125},
  {173,127,11},
  {171,164,157},
  {105,143,108},
  {179,149,74},
  {143,118,70},
  {127,130,140},
  {105,41,91},
  {8,34,89},
  {169,74,17},
  {96,4,60},
  {93,118,7},
  {86,69,45},
  {161,94,119},
  {145,25,81},
  {115,28,79},
  {52,14,8},
  {103,98,96},
  {45,46,46},
  {50,82,118},
  {69,142,93},
  {144,24,80},
  {175,134,105},
  {56,22,72},
  {75,45,30},
  {47,62,150},
  {137,86,57},
  {174,114,96},
  {62,52,114},
  {91,87,142},
  {118,136,19},
  {98,138,125},
  {150,155,83},
  {173,102,66},
  {155,151,135},
  {162,95,12},
  {89,120,147},
  {116,23,15},
  {177,173,89},
  {72,127,133},
  {7,19,34},
  {6,29,10},
  {26,9,6},
  {64,106,132},
  {75,4,15},
  {6,82,121},
  {142,162,148},
  {83,55,51},
  {171,78,105},
  {81,54,74},
  {58,77,55},
  {58,67,79},
  {114,60,41},
  {169,79,18},
  {102,34,14},
  {159,154,147}
};

void RunLegoMosaic( string inputName, string outputName, int numHorizontal, int tileSize, bool sidesOut, bool dither, double gamma, string colorName );