#include "textureCompletion2.h"

using namespace cv;
using namespace std;

int main() 
{
	Mat img = imread("./texture_test_img/sp6.png");
	Mat srcLineMask  = imread("./texture_test_img/mask_s6.png");
	if (img.empty() || srcLineMask.empty()) {
		cout << "请检查文件名称是否有误！" << endl;
		return -1;
	}
	
    TextureCompletion origin ;
    
    Mat gray1,gray2;
	cvtColor(img, gray1, COLOR_BGR2GRAY);
	cvtColor(srcLineMask, gray2, COLOR_BGR2GRAY);
	//通过阈值处理生成Mask掩码
	
	Mat1b LineMask,imgMask;
	threshold(gray1, imgMask, 0, 255,  THRESH_BINARY_INV);
	threshold(gray2, LineMask, 0, 255,  THRESH_BINARY_INV);
	// imshow("imgMask", imgMask);
	// imshow("LineMask", LineMask);

	Mat finalresult;
    origin.getTexture(img , imgMask, LineMask,finalresult);
	imshow("final", finalresult);
	waitKey(0);
	return 0;
}

int sqr(int x)
{
	return x * x;
}

int dist(Vec3b V1, Vec3b V2)
{
	return sqr(int(V1[0]) - int(V2[0])) + sqr(int(V1[1]) - int(V2[1])) + sqr(int(V1[2]) - int(V2[2]));
}
void TextureCompletion::getTexture(Mat srcMat, Mat1b mask, Mat1b LineMask, Mat &result)
{								// srcMat是带LineMask绘制的原图  mask是带线的遮挡部分
	cout << "Texture Synthesis begin!" << endl;
	result = srcMat.clone();
	int N = mask.rows;
	int M = mask.cols;

	vector<vector<int> >my_mask(N, vector<int>(M, 0)), sum_diff(N, vector<int>(M, 0));

	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			//mymask对应于mask（mask中的白色遮挡部分为待填充部分mymask为1，
			//mask黑色部分为无需修改部分，mymask为0
			my_mask[i][j] = (mask.at<uchar>(i, j) == 255);
			//如果LineMask中的灰色部分，则标注为2	
			if (LineMask.at<uchar>(i, j) == 0)
			{
				my_mask[i][j] = 2;
			}
		}
	/*
	my_mask的结构
	0 0 0 0 0 0 0
	0 0 1 1 1 0 0
	0 1 1 1 1 1 0
	0 2 2 2 2 2 0  ---结构线用2表示，白色等待填充部分为1，黑色不填充为0
	0 1 1 1 1 0 0  
	0 0 1 1 0 0 0
	0 0 0 0 0 0 0
	*/
	int bs = 3;
	int step = 6 * bs;
	int to_fill = 0;	
	int filled = 0;		
	auto usable(my_mask);	//自动生成了一个和mymask相同类型的变量
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			to_fill += (my_mask[i][j] == 1);
		}

	for (int i = 0; i < N; i++){
		for (int j = 0; j < M; j++)
		{
			//遍历全图，如果my_mask[i][j] == 0说明不需要填充则继续
			if (my_mask[i][j] == 0)
				continue;
			//对于mymask中需要被填充的地方
			//在一个step的矩形邻域内，需要把usable标记为2
			//usable[k][l] == 2说明需要被填充
			//（我的理解是在原来的mask周围扩大了需要补全纹理的范围，缩小了可用的纹理的范围）
			int k0 = max(0, i - bs), k1 = min(N - 1, i + bs);
			int l0 = max(0, j - bs), l1 = min(M - 1, j + bs);
			for (int k = k0; k <= k1; k++)
				for (int l = l0; l <= l1; l++)
					usable[k][l] = 2;
		}

	}

	//按照usable中2的地方生成一个黑白图，其中白色是需要填充的地方值为2
	//也就是说实际要填充的部分比正常的黑白图要大
	Mat use = mask.clone();
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
			if (usable[i][j] == 2)
				use.at<uchar>(i, j) = 255;
			else use.at<uchar>(i, j) = 0;

			int itertime = 0;
			Mat match;
			Mat output;
		    match = result.clone();
			while (true)
			{
				itertime++;
				int x, y, cnt = -1;
				for (int i = 0; i < N; i++)
					for (int j = 0; j < M; j++)//找到边界点
					{
						//略过不需要填充的地方以及轮廓线部分
						if (my_mask[i][j] != 1) continue;
						//此时my_mask[i][j] == 1 是白色待填充部分
						//首先要找到需要填充的区域的边界点
						//edge用于判断这个点是不是边界
						bool edge = false;
						int k0 = max(0, i - 1), k1 = min(N - 1, i + 1);
						int l0 = max(0, j - 1), l1 = min(M - 1, j + 1);
						//取到像素点的一个小邻域8个像素点，如果这个邻域内的点有一个是0则最后edge==true
						/*
						0 0 0
						0 1 0
						0 0 0
						*/
						for (int k = k0; k <= k1; k++)
							for (int l = l0; l <= l1; l++)
								edge |= (my_mask[k][l] == 0);	//或等于 edge = edge | (my_mask==0);
						if (!edge) continue;//如果不是边界点则继续找
						//如果edge==true说明当前像素点是边界点
						k0 = max(0, i - bs), k1 = min(N - 1, i + bs);
						l0 = max(0, j - bs), l1 = min(M - 1, j + bs);
						int tmpcnt = 0;
						//此时取到当前像素点周围的一个step大小的矩形邻域
						//tmpcnt计算了这个矩形邻域内不需要填充的像素点的个数
						for (int k = k0; k <= k1; k++)
							for (int l = l0; l <= l1; l++)
								tmpcnt += (my_mask[k][l] == 0);
						if (tmpcnt > cnt)
						{
							cnt = tmpcnt;
							x = i;
							y = j;
						}
						//结束for循环的时候xy记录了边界点
					}
				//如果cnt==-1说明所有edge都是false，也就是说所有mymask[i，j]都是0都是不需要填充，跳出while
				if (cnt == -1) break;

				bool debug = false;
				bool debug2 = false;


				//再次遍历全图；比较一个邻域内和整张图片其他邻域内是否有相似的块
				int k0 = min(x, bs), k1 = min(N - 1 - x, bs);
				int l0 = min(y, bs), l1 = min(M - 1 - y, bs);
				//这里使用p0q0使得本身就在对应点的邻域寻找
				int p0 = max(x - step, bs), p1 = max(N - 1 - x - step, bs);
				int q0 = max(y - step, bs), q1 = max(M - 1 - y - step, bs);
				int p2 = min(x + step, N);
				int q2 = min(y + step, M);
				int sx = 1000000;
				int sy = 1000000;
				int min_diff = 1000000;	//最大的int值
				for (int j = q0; j + bs < M-q1; j += bs)
					for (int i = p0; i + bs < N-p1; i += bs)
					{
				
						//printf("%d\n", tmp);
						//通过usable找到最近的不需要填充的像素点
						//如果==2说明这里没有纹理
						//match.at<Vec3b>(i, j) = Vec3b(255, 0, 0);
						if (my_mask[i][j] == 2) {
							break;
						}
						if (usable[i][j] == 2)	continue;//需要填充
						
						int tmp_diff = 0;
						//取到xy和ij周围step的矩形邻域
						for (int k = -k0; k <= k1; k++)
							for (int l = -l0; l <= l1; l++)
							{
								//printf("%d %d %d %d %d %d\n", i + k, j + l, x + k, y + l, N, M);
								//ij表示可以用来比较的不需要填充纹理的坐标点
								//xy表示当前需要被填充的点，由之前的for循环生成
								//[x + k][y + l]表示xy的step邻域内的某点
								//[i + k][j + l]表示ij的step邻域内的某点
								if (my_mask[x + k][y + l] != 1)   
									tmp_diff += dist(result.at<Vec3b>(i + k, j + l), result.at<Vec3b>(x + k, y + l));
								//tmp_diff计算了这两个对应点之间，RGB值的差异；显然需要全图搜索找到一个最小的tmpdiff，这说明这两块邻域最像

							}
						//printf("tmp_diff = %d", tmp_diff);
						sum_diff[i][j] = tmp_diff;
						if (min_diff > tmp_diff)
						{	
							sx = i;
							sy = j;
							min_diff = tmp_diff;
						}
						sum_diff[i][j] = tmp_diff;
						
						//结束循环的时候，得到的是对比xy有最小tmpdiff的点的坐标sx，sy
						
					}
				// imshow("iii", match);
				// waitKey(10);
				
				cout << "current xy：" << x << y << endl;
				if (sx == 1000000 && sy == 1000000) {
					//这种点实际上特别多！！！要保证可以获取到能用的texture！！
					//这里已经是触发异常的点，进行全局搜索,
					cout << "exception xy" << endl;
					match.at<Vec3b>(x, y) = Vec3b(0, 0, 255);
					for (int j = M - step; j - bs > step; j -= bs)
						for (int i = N - step; i - bs > step; i -= bs)
						
						{
							int tmp_diff = 0;
							/*if (my_mask[i][j] == 2) {
								cout << i << " , " << j << endl;
								break;
							}*/
							if (usable[i][j] == 2)	continue;
							
							for (int k = -k0; k <= k1; k++)
								for (int l = -l0; l <= l1; l++)
									if (my_mask[x + k][y + l] != 1)
										tmp_diff += dist(result.at<Vec3b>(i + k, j + l), result.at<Vec3b>(x + k, y + l));
							sum_diff[i][j] = tmp_diff;
							if (min_diff > tmp_diff)
							{

								sx = i;
								sy = j;
								min_diff = tmp_diff;
							}
							sum_diff[i][j] = tmp_diff;
						}
					if (usable[sx][sy] == -1) {

						printf("------have already filled -----");
					}
					
				}
				if (sx == 1000000 && sy == 1000000) {
					sx = x;
					sy = y;
					printf("still no");
				}

				cout << "dui ying de dian xy " << sx << sy << endl;

				usable[x][y] = -1;
				//用（sx，sy）周围的点的RGB值填充xy周围需要被填充的点
				for (int k = -k0; k <= k1; k++)
					for (int l = -l0; l <= l1; l++)
						if (my_mask[x + k][y + l] == 1)
						{
							result.at<Vec3b>(x + k, y + l) = result.at<Vec3b>(sx + k, sy + l);
							my_mask[x + k][y + l] = 0;
							//usable[x + k][y + l] = 1;
							filled++;
						}

				printf("filled: %d , to-fill: %d \n", filled ,to_fill);
				imwrite("final.png", result);
				// imshow("final", result);
				// waitKey(0);
			}
}