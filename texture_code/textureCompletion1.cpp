#include <iostream>
#include <opencv2/highgui.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp> 
#include <opencv2/opencv.hpp>
using namespace cv;
using namespace std;
Mat TextureCompletion1(Mat img);
int main() // sp3 , sp5 表现还可以
{
	Mat img = imread("./texture_test_img/sp3.png"); // img是画好线的图
	if (img.empty()) {
		cout << "请检查文件名称是否有误！" << endl;
		return -1;
	}
	
	// imshow("oringin img", img);
 
	//图像修复
	Mat imgInpaint;
	imgInpaint = TextureCompletion1(img);

	return 0;
}

Mat TextureCompletion1(Mat img){
	Mat gray;
	cvtColor(img, gray, COLOR_BGR2GRAY);
  
	//通过阈值处理生成Mask掩码
	Mat imgMask;
	// threshold(gray, imgMask, 245, 255, THRESH_BINARY);
	threshold(gray, imgMask, 0, 255,  THRESH_BINARY_INV);
	//对Mask掩码膨胀处理，增加Mask的面积
	Mat Kernel = getStructuringElement(MORPH_RECT, Size(3, 3));
	dilate(imgMask, imgMask, Kernel);
 
	//图像修复
	Mat imgInpaint;
	inpaint(img, imgMask, imgInpaint, 3, INPAINT_TELEA);
 	// inpaint(img, imgMask, imgInpaint, 5, INPAINT_NS);
	//显示处理结果
	// imshow("imgMask", imgMask);
	imshow("After Inpaint", imgInpaint);
	waitKey(0);
	return imgInpaint;
}