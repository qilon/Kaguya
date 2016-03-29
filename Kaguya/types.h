#include <vector>

typedef double Coordinate;
typedef double Intensity;
typedef unsigned int VertexIndex;
typedef unsigned int FaceIndex;

typedef std::vector<Coordinate> Vertex;
typedef std::vector<Coordinate> Normal;
typedef std::vector<Intensity> Color;
typedef std::vector<VertexIndex> Face;

#define RGB2GRAY_R_FACTOR 0.299
#define RGB2GRAY_G_FACTOR 0.587
#define RGB2GRAY_B_FACTOR 0.114
