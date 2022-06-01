#ifndef _TEXTURESYNTHESIS_H
#define _TEXTURESYNTHESIS_H


#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <ctime>
#include <math.h>
#include <opencv2/highgui.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp> 
using namespace cv;

#include <iostream>
#include <fstream>

using namespace std;

class TextureCompletion
{

private:


public:
	TextureCompletion(){}
	~TextureCompletion(){}
	void getTexture(Mat srcMat, Mat1b mask, Mat1b LineMask, Mat &result);

};

#endif