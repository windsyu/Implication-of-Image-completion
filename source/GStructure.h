#ifndef _GRAPH_H_
#define _GRAPH_H_
#include<opencv2/opencv.hpp>
#include<iostream>

using namespace cv;
using namespace std;
//需要调试
#define ks 50
#define ki 2
#define PatchSizeRow 13
#define PatchSizeCol 13

// BP and DP Realization
typedef enum{
    MISS,
    KNOWN,
    BORDER
}PointType;
//Record all the points in an operating patch(左端点、 右端点 和 中心锚点)
class Anchor{
private:

public:
    int head_point, end_point,  center_point;
//indicate the position of points    
    PointType type;
//the record of the points n
    vector <int > neighbors;
//Constructor
    Anchor(int head, int end, int center, PointType t);
    ~Anchor();
};




class GStructure{
private:
//Graph_Data
    Mat image, mask, regions, image_loss;
//Curve Record
    vector < vector<Anchor> > sample_anchor, unknown_anchor;
    vector < vector<Point2i> > lines;
//Beginning with the end_anchor, return the point of the last index inthe patch and the type of it
    int PutPointInPatch(int end_anchor, int anchor, PointType& type, int lineindex);
//get the nearly index for every unknown anchor
    void InputNeighbor(int index);
//Return the head point of the certain patch
    Point2i LeftTopPoint(Anchor anchor, int lineindex);
//Is Intersect?
    bool IsIntersect(int a1, int a2);
//Return the patch of the image with a given anchor
    Mat GetPatch(Anchor anchor, int lineindex);
//Mininize the Energy ,  unknown anchor , sample anchor, point index
    double EI(Anchor unknow, Anchor sample, int p_index);
    double ES(Anchor unknow, Anchor sample, int p_index);
    double E1(Anchor unknow, Anchor sample, int p_index);
    double E2(Anchor unknow1, Anchor unknow2, Anchor sample1, Anchor sample2, int p_index);

//Dynamic Program method
    vector<int > DP(vector<Anchor > unknow, vector<Anchor > sample, int p_index);

//Belief Propagation
    vector<int > BP(vector<Anchor > unknow, vector<Anchor > sample, int p_index);
//Complete image along the given curve
    void CompleteLine(vector<Anchor > sample, vector<Anchor> unknow, int p_index, bool flag);
//
    void ComputeAnchor(int p_index, vector<Anchor >& sample, vector<Anchor >& unknow);

//Given unknow 
    void CopyPic(Anchor unknow, Mat patch, Mat fullPic, int p_index);
public:
    GStructure(Mat input);
    ~GStructure(){}
//Add the unknown area
    void AddMask();
//put sample points on the curve into the sample anchor list 
//record the result in vector<vector<Anchor > > sample
    void AddAnchor();
//input and record the given curve
    void AddCurve();
//Specify the numbers of regions
    void DrawNewRegion();
};


#endif