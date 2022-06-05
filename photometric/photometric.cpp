#define _CRT_SECURE_NO_WARNINGS
#include "photometric.h"
#define CNT_THRESHOLD 1000
#define DIFF_THRESHOLD 1e-15
#define OMEGA 1.8
using namespace cv;
using namespace std;

int offset[4][2] = { {1, 0}, {-1, 0}, {0, 1}, {0, -1} };
Mat Photometric::dst;
Mat Photometric::mask;

void Photometric::initMask(Mat image, Mat mask) {
	Mat temp;
	dst = Mat(mask.size().height, mask.size().width, CV_64FC3);
	image.convertTo(temp, CV_64FC3);
	temp.copyTo(dst);
	mask = Mat(mask.size().height, mask.size().width, CV_8U);
	mask.setTo(Scalar(0));
	mask.copyTo(mask);
	Mat unknown = mask = 0;
	Mat known = mask = 255;
	mask.setTo(Scalar(2), unknown);
	mask.setTo(Scalar(3), known);
	//imshow("_mask", mask);
	waitKey(0);
}

void Photometric::modify(Mat& patch, int xoffset, int yoffset) {
	int width, height;
	width = patch.size().width;
	height = patch.size().width;
	Mat temp;
	patch.convertTo(temp, CV_64FC3);
	Mat result = Mat(height, width, CV_64FC3);
	Mat source = Mat(height, width, CV_64FC3);
	temp.copyTo(result);
	temp.copyTo(source);
	Mat bitmap = Mat(height, width, CV_8U);
	Rect rect = Rect(xoffset, yoffset, width, height);
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			if (mask.at<uchar>(y + yoffset, x + xoffset) == 0) {
				result.at<double>(y, x) = dst.at<double>(y + yoffset, x + xoffset);
				source.at<double>(y, x) = dst.at<double>(y + yoffset, x + xoffset);
				bitmap.at<uchar>(y, x) = 0;
			}
			else if (mask.at<uchar>(y + yoffset, x + xoffset) == 2) {
				mask.at<uchar>(y + yoffset, x + xoffset) = 1;
				bitmap.at<uchar>(y, x) = 1;
			}
			if (x == 0 || y == 0 || x == width - 1 || y == height - 1)
				mask.at<uchar>(y + yoffset, x + xoffset) = 3;
		}
	}
	int cnt = 0;
	while (cnt < CNT_THRESHOLD) {
		double difference = 0, abstot = 0, previous_epsilon = 1e30;
		for (int x = 1; x < width - 1; x++) {
			for (int y = 1; y < height - 1; y++) {
				int adjacent = 0;
				Vec3d sum_fq = Vec3d(0, 0, 0);
				Vec3d sum_vpq = Vec3d(0, 0, 0);
				Vec3d sum_boundary = Vec3d(0, 0, 0);
				Vec3d new_val = Vec3d(0, 0, 0);
				for (int i = 0; i < 4; i++) {
					if (mask.at<uchar>(y + yoffset + offset[i][0], x + xoffset + offset[i][1]) == 0
						|| mask.at<uchar>(y + yoffset + offset[i][0], x + xoffset + offset[i][1]) == 1) {
						sum_fq += result.at<double>(y + offset[i][0], x + offset[i][1]);
						if (mask.at<uchar>(y + yoffset, x + xoffset) == mask.at<uchar>(y + yoffset + offset[i][0], x + xoffset + offset[i][1]))
							sum_vpq += source.at<double>(y, x) - source.at<double>(y + offset[i][0], x + offset[i][1]);
						adjacent++;
					}
					if (mask.at<uchar>(y + yoffset + offset[i][0], x + xoffset + offset[i][1]) == 3) {
						sum_boundary += source.at<double>(y + offset[i][0], x + offset[i][1]);
						if (bitmap.at<uchar>(y, x) == bitmap.at<uchar>(y + offset[i][0], x + offset[i][1]))
							sum_vpq += source.at<double>(y, x) - source.at<double>(y + offset[i][0], x + offset[i][1]);
						adjacent++;
					}
				}
				new_val = (sum_fq + sum_vpq + sum_boundary) / (1.0 * adjacent);
				difference += abs(new_val(0) - result.at<double>(y, x)(0))
					+ abs(new_val(1) - result.at<double>(y, x)(1))
					+ abs(new_val(2) - result.at<double>(y, x)(2));
				abstot += abs(new_val(0)) + abs(new_val(1)) + abs(new_val(2));
				result.at<double>(y, x) = (1 - OMEGA) * result.at<double>(y, x) + OMEGA * new_val;
			}
			cnt++;
			if (abs(difference / abstot) < DIFF_THRESHOLD || abs(previous_epsilon - difference / abstot) < 1e-12)	break;
			else previous_epsilon = difference / abstot;
		}
	}
	mask(rect).setTo(Scalar(0));
	Mat _result;
	result.convertTo(_result, CV_8UC3);
	_result.copyTo(patch);
	result.copyTo(dst(rect));
	return;
}