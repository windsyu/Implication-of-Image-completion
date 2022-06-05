#ifndef _PHOTOMETRIC_H
#define _PHOTOMETRIC_H

#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace cv;
using namespace std;

class Photometric {
private:

public:
	static Mat dst;
	static Mat mask;
	Photometric(){}
	~Photometric(){}
	static void initMask(Mat image, Mat mask);
	static void modify(Mat& patch, int xoffset, int yoffset);
};

#endif
