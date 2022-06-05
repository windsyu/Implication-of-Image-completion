#include"GStructure.h"

static Point2i m_position;
static int m_event;
static int m_flag;


//Visualize the anchor points only for debug
void DrawAnchor(Mat input, vector<Anchor > sample_patch, vector<Anchor > unknow_patch, vector<Point2i> points){
    Point2i position;
    for(int i = 0 ; i < unknow_patch.size() ; i++ ){
        position = points[unknow_patch[i].center_point];
        //边界标注为红色
        if(unknow_patch[i].type == BORDER){
            cout << "border:" <<position <<endl;
            int left_x = position.x-(PatchSizeRow-1)/2;
            int top_y = position.y-(PatchSizeCol-1)/2;
            Rect rect(left_x, top_y, PatchSizeRow, PatchSizeCol);
            circle(input,position,3,Scalar(0,255,0),-1);
            rectangle(input, rect, Scalar(0,0,0),1);
        }else
            circle(input,position,1,Scalar(0,0,255),-1);
    }
    for(int i =0 ; i < sample_patch.size(); i++){
        position = points[sample_patch[i].center_point];
        int left_x = position.x - (PatchSizeRow - 1) /2;
        int top_y = position.y - (PatchSizeCol - 1) /2;
        Rect rect(left_x, top_y, PatchSizeRow, PatchSizeCol);
        circle(input,position,1,Scalar(0,0,255),-1);
    }

    imshow("Anchor points", input);
    waitKey(0);
}

static void MouseBack(int event, int x, int y, int flag, void* ustate){
    m_event = event;
    m_position = Point(x,y);
    m_flag = flag;
}

GStructure::GStructure(Mat input):image(input){
    this->mask = Mat::ones(input.size(),CV_8U);
    this->mask.setTo(255);
    image_loss = image.clone();
}

void GStructure::AddMask(){
    Mat toshow = this->image.clone();
    int size = 20;
    imshow("Generate Input Mask", toshow);
    setMouseCallback("Generate Input Mask", MouseBack);
    Point2i last_position(-1,-1);
    while(1){
        Mat totoshow = toshow.clone();
        char key = waitKey(10);
        if(key == '['){
            cout<<"Smaller"<<endl;
            size = (size>1)? size-1 :size;
        }else if(key == ']'){
            cout<<"Bigger"<<endl;
            size++;
        }else if(key == 27){
            cout<<"break"<<endl;
            break;
        }else if((m_event == EVENT_MOUSEMOVE &&(m_flag & EVENT_FLAG_LBUTTON)) || 
            (m_event == EVENT_LBUTTONDOWN)){
            if(last_position.x != -1 && last_position.y != -1){
                line(this->mask,last_position,m_position,Scalar(0),1.5*size);
                line(toshow,last_position,m_position,Scalar(255,0,0),1.5*size);
                line(this->image_loss,last_position,m_position,Scalar(255,0,0),1.5*size);
            }
            last_position=m_position;
        }else{
            last_position.x=-1;
            last_position.y=-1;
        }
        circle(totoshow,m_position,size,Scalar(255,0,0),-1);
        imshow("Generate Input Mask", totoshow);
    }
//  Photometric::initMask(this->image_loss,this->mask);
    destroyWindow("Generate Input Mask");
    return ;
}



void GStructure::AddCurve(){
    Mat toshow = this->image_loss.clone();
    int size = 5;
    namedWindow("Generate Specified Curves", 1);
    imshow("Generate Specified Curves", toshow);
    setMouseCallback("Generate Specified Curves", MouseBack);
    Point2i last_pos(-1,-1);
    Point2i last_point(-1,-1);
    vector<Point2i> vec;

    while(1){
        Mat totoshow = toshow.clone();
        char key = waitKey(10);
        if(key == 27){
            break;
        }else if(key == 32){
            cout<<"points num:"<<vec.size()<<endl;
            this->lines.push_back(vec);
            vec.clear();
            cout << this-> lines.size() <<endl;
        }else if((m_event == EVENT_MOUSEMOVE && (m_flag == EVENT_FLAG_LBUTTON)) ||
            (m_event == EVENT_LBUTTONDOWN))
        {
            if(last_pos != Point2i(-1,-1)){
                line(toshow,last_pos,m_position,Scalar(0,255,255),size*1.5);
                LineIterator it(toshow, last_pos, m_position, 8);

                for(int i = 0; i < it.count; i++,++it){
                    if(last_point == Point2i (-1,-1) || last_point != it.pos()){
                        vec.push_back(it.pos());
                    }
                }
            }
            last_pos = m_position;
        }else{
            last_pos = Point2i(-1,-1);
        }
        circle(totoshow, m_position,size,Scalar(0,255,255),-1);
        imshow("Generate Specified Curves", totoshow);
    }
    destroyWindow("Generate Specified Curves");
}

int GStructure::PutPointInPatch(int last_anchor, int anchor, PointType& type,int lineindex){
    int left_x = this->lines[lineindex][anchor].x - (PatchSizeCol-1)/2;
    int top_y = this->lines[lineindex][anchor].y - (PatchSizeRow-1)/2;
    int i = last_anchor;
    Rect rec(left_x,top_y,PatchSizeRow,PatchSizeCol);
    //Enter Patch area;
    while(!rec.contains(this->lines[lineindex][i])){
        i++;
    }
    if(this->mask.at<uchar>(this->lines[lineindex][i]) == 0)
        type = MISS;
    else
        type = KNOWN;
    while(i < this->lines[lineindex].size() && rec.contains(this->lines[lineindex][i])){
        if((this->mask.at<uchar>(this->lines[lineindex][i]) == 0 && type == KNOWN) ||
            (this->mask.at<uchar> (this->lines[lineindex][i])  > 0 && type == MISS))
            type == BORDER;

        i++;
    }
    if(i<anchor){
        cout << "GetPatch() Error" <<endl;
        return anchor;
    }
    return i;
}
void GStructure::ComputeAnchor(int p_index, vector<Anchor >& sample, vector<Anchor >& unknow){
    PointType type;
    int last_anchor_point = 0;
    int anchor_point = PutPointInPatch(last_anchor_point,0,type,p_index);
    int next_anchor_point = PutPointInPatch(last_anchor_point,anchor_point,type,p_index);
    while(next_anchor_point < this->lines[p_index].size()){
        Anchor anchor(last_anchor_point + 1, anchor_point, next_anchor_point - 1, type);
        if(anchor.type == KNOWN)
            sample.push_back(anchor);
        else
            unknow.push_back(anchor);
        
        last_anchor_point = anchor_point;
        anchor_point = next_anchor_point;
        next_anchor_point = PutPointInPatch(last_anchor_point, anchor_point, type, p_index);
    }
    cout<<"sample anchor num:" << sample.size() <<endl;
    cout<<"unknow anchor num" << sample.size() <<endl;
}

void GStructure::AddAnchor(){
    AddCurve();
    vector<Anchor> sample;
    vector<Anchor> unknown;
    for(int i = 0; i< this-> lines.size(); i++){
        ComputeAnchor(i,sample,unknown);
        this->sample_anchor.push_back(sample);
        this->unknown_anchor.push_back(unknown);
        sample.clear();
        unknown.clear();
    }

}
double Dist(vector<Point2i> points1, vector<Point2i> points2){
    double result= 0.0;
    int i,j;
    double shortest;
    double tmp;
    double normalized = (double) norm(Point2i(PatchSizeRow,PatchSizeCol));

    for(i = 0; i<points1.size(); i++){
        shortest = INFINITY;
        for(j = 0; j < points2.size(); j++){
            tmp = (double) norm(points1[i] - points2[j]) / normalized;
            tmp *= tmp;
            if(tmp < shortest)
                shortest = tmp;
        }
        result += shortest;
    }

    return result;
}
Point2i GStructure::LeftTopPoint(Anchor anchor, int p_index){
    Point2i p;
    int index = anchor.center_point;
    p.x = this->lines[p_index][index].x - (PatchSizeCol-1) / 2;
    p.y = this->lines[p_index][index].y - (PatchSizeRow-1) / 2;
    return p;
}

double SSD(Mat m1, Mat m2){
    double temp;
    Mat result(1,1,CV_32F);
    if(m1.empty() || m2.empty()){
        cout << "SSD computing non-overlapped area" <<endl;
        return 0.0;
    }
    matchTemplate(m1,m2,result,TM_SQDIFF_NORMED);
    temp = result.at<double>(0,0);
    return temp;
}

Mat GStructure::GetPatch(Anchor anchor, int lineindex){
    Mat patch;
    Point2i left_top = LeftTopPoint(anchor, lineindex);
    Point2i right_down = left_top + Point2i(PatchSizeCol, PatchSizeRow);
    Rect rect(left_top, right_down);
    if(left_top.x < 0 || left_top.y < 0 || right_down.x > this->image.cols || right_down.y > this->image.rows){
        cout<< left_top.x<<" "<<left_top.y<<" "<< right_down.x<< " "<< right_down.y<<endl;
    }
    this->image(rect).copyTo(patch);
    return patch;
}

double GStructure::ES(Anchor unknow, Anchor sample, int p_index){
    double result = 0.0;
    int p_num = unknow.end_point - unknow.head_point + 1;

    Point2i unknown_left_top, sample_left_top;
    unknown_left_top = LeftTopPoint(unknow, p_index);
    sample_left_top = LeftTopPoint(sample, p_index);
    vector<Point2i> points1;
    vector<Point2i> points2;
    for(int i = unknow.head_point ; i<= unknow.end_point; i++){
        if((this->lines[p_index][i].x - unknown_left_top.x < 0) || (this->lines[p_index][i].y - unknown_left_top.y < 0)){
            points1.push_back(Point2i(0,0));
        }else{
            points1.push_back(this->lines[p_index][i] - unknown_left_top);
        }
    }
    for(int i = sample.head_point ; i <= sample.end_point; i++){
        if((this->lines[p_index][i].x - sample_left_top.x < 0 || this->lines[p_index][i].y - sample_left_top.y <0 ))
            points2.push_back(Point2i(0,0));
        else
            points2.push_back(this->lines[p_index][i] - sample_left_top);
    }
    result = Dist(points1,points2);
    result +=Dist(points2,points1);
    result /= p_num;


    return result;
}
double GStructure::EI(Anchor unknow, Anchor sample, int p_index){
    if(unknow.type != BORDER)
        return 0.0;
    Mat patch_image, patch_mask;
    Mat patch_temp, sample_temp;
    Point2i left_top = LeftTopPoint(unknow, p_index);
    Point2i right_down = left_top + Point2i(PatchSizeCol, PatchSizeRow);
    Rect rect(left_top, right_down);
    this->image(rect).copyTo(patch_image);
    this->mask(rect).copyTo(patch_mask);
    patch_image.copyTo(patch_temp, patch_mask);
    Mat sample_image = GetPatch(sample, p_index);
    sample_image.copyTo(sample_temp, patch_mask);
//
//
    double result = SSD(patch_temp, sample_temp);
//
    return result;

}

double GStructure::E1(Anchor unknow, Anchor sample, int p_index){
    return ks*ES(unknow,sample, p_index) + ki*EI(unknow, sample, p_index);
}
double GStructure::E2(Anchor unknow1, Anchor unknow2, Anchor sample1, Anchor sample2, int p_index){
    Point2i A, B, C, D;
    A = LeftTopPoint(unknow1,p_index);
    C = LeftTopPoint(unknow2,p_index);
    B = A + Point2i(PatchSizeCol, PatchSizeRow);
    D = C + Point2i(PatchSizeCol, PatchSizeRow);
    Rect rect1(A, B);
    Rect rect2(C, D);
    Rect rect_intersect = rect1 & rect2;
    Mat part1 = GetPatch(sample1, p_index);
    Mat part2 = GetPatch(sample2, p_index);

    Mat refer1(this->image.size(), this->image.type(), Scalar(0,0,0));
    Mat refer2(this->image.size(), this->image.type(), Scalar(0,0,0));

    part1.copyTo(refer1(rect1));
    part2.copyTo(refer2(rect2));
    Mat patch1, patch2;
    refer1(rect_intersect).copyTo(patch1);
    refer2(rect_intersect).copyTo(patch2);

    if(patch1.empty()){
        cout << "E2 compute error"<<endl;
        return -1;
    }

    return SSD(patch1, patch2);
        
    
}
vector<int> GStructure::DP(vector<Anchor > unknow, vector<Anchor > sample, int p_index){
    vector<int> sample_index;
    int sample_size = sample.size();
    int unknow_size = unknow.size();
    double **M = new double*[unknow_size];
    int **index = new int*[unknow_size];
	for (int i = 0;i < unknow_size;i++) {
		M[i] = new double[sample_size];
		index[i] = new int[sample_size];
	}
		
	for (int i = 0;i < sample_size;i++) {
		M[0][i] = E1(unknow[0], sample[i], p_index);
	}

	int i, j, k;
	double e1;
	double min_tmp=INFINITY;
	double tmp;
	int index_tmp;
	for ( i = 1;i < unknow_size;i++) {
		for (j = 0;j < sample_size;j++) {
			min_tmp = INFINITY;
			e1 = E1(unknow[i], sample[j], p_index);
			for (k = 0;k < sample_size;k++) {
				tmp = M[i - 1][k]+ E2(unknow[i - 1],
					unknow[i], sample[k], sample[j], p_index);
				if (tmp < min_tmp) {
					min_tmp = tmp;
					index_tmp = k;
				}
					
			}
			M[i][j] = e1 + min_tmp;
			index[i][j] = index_tmp;
		}
	}

	min_tmp = INFINITY;
	for (i = 0;i < sample_size;i++) {
		tmp = M[unknow_size - 1][i];
		if (tmp < min_tmp) {
			min_tmp = tmp;
			index_tmp = i;
		}
	}

//	cout << "The min energy is:" << min_tmp << endl;
	sample_index.push_back(index_tmp);
	for (i = unknow_size - 1;i > 0;i--) {
		index_tmp = index[i][index_tmp];
		sample_index.push_back(index_tmp);
	}
	reverse(sample_index.begin(), sample_index.end());
	cout << "sample_index_size:" << sample_index.size() << endl;
	//for (int i = 0;i < sample_index.size();i++)
	//	cout << sample_index[i] << endl;
	cout << "DP done" << endl;
	return sample_index;
}

void GStructure::CopyPic(Anchor unknow, Mat patch, Mat fullPic, int p_index){
    Point2i left_top = LeftTopPoint(unknow, p_index);
    Rect rect(left_top.x, left_top.y, PatchSizeRow, PatchSizeCol);
    patch.copyTo(fullPic(rect));

}
bool IsClose(Point2i A, Point2i B){
    if(norm(A - B) <norm(Point2i(PatchSizeCol/2, PatchSizeRow/2)))
        return true;
    return false;
}
bool GStructure::IsIntersect(int a1, int a2){
    	if (this->unknown_anchor[a1].empty() || this->unknown_anchor[a2].empty())
		return false;
	int point1, point2;
	int j, k;
	int i1_num = this->unknown_anchor[a1].size();
	for (j = 0; j < i1_num ;j++) {
		point1 = this->unknown_anchor[a1][j].center_point;
		for (k = 0;k < this->unknown_anchor[a2].size();k++) {
			point2 = this->unknown_anchor[a2][k].center_point;
			if (IsClose(this->lines[a1][point1], this->lines[a2][point2])) {
				this->unknown_anchor[a1][j].neighbors.push_back(k + i1_num);
				this->unknown_anchor[a1][j].neighbors.push_back(k + 1 + i1_num);
				//afterwards(in DrawNewStructure()) will add i1_num
				this->unknown_anchor[a2][k].neighbors[1] = j - i1_num;
				this->unknown_anchor[a2][k + 1].neighbors[0] = j - i1_num;
				//circle(this->image_with_mask, this->points[i1][point1], 2, Scalar(0, 0, 0));
				return true;
			}	
		}
	}
	//imshow("intersection", this->image_with_mask);
	//waitKey(0);
	return false;
}
//Belief Propagation
//
double* AddTwoVec(double* a1, double* a2, int n){
    double* re = new double[n];
	for (int i = 0;i < n; i++) {
		re[i] = a1[i] + a2[i];
	}
	return re;
}
double* MinTwoVec(double* a1, double* a2, int n){
    double* re = new double[n];
    for(int i = 0; i<n; i++){
        re[i] = a1[i] - a2[i];
    }
    return re;
}
void initialize(double* a, int n){
	for (int i = 0;i < n;i++) {
		a[i] = 0;
	}
}
vector<int > GStructure::BP(vector<Anchor > unknow, vector<Anchor > sample, int p_index){
    cout<<"BP begin. Finish Soon..."<<endl;
    vector<int> sample_index;
    int unknow_size = unknow.size();
    int sample_size = sample.size();
// M三维矩阵 
    double ***M = new double **[unknow_size];
//e1是一个二维矩阵
    double **e1 = new double*[unknow_size];
    bool** converged = new bool*[unknow_size];
    bool tmp_flag;
	for (int i = 0;i < unknow_size; i++) {
		M[i] = new double *[unknow_size];
		converged[i] = new bool[unknow_size];
		e1[i] = new double[sample_size];
		for (int j = 0;j < unknow_size;j++) {
			M[i][j] = new double[sample_size];
		}
	}

    int i, j, k, x, y;

	double* neighM_sum = new double[sample_size];
	double* M_tmp = new double[sample_size];
	double tmp_min=INFINITY;
	double* neighM = new double[sample_size];    
//初始化 M_tmp, sample_size
	for (i = 0;i < unknow_size;i++) {
		for (j = 0;j < unknow_size;j++) {
			converged[i][j] = false;
			for (k = 0;k < sample_size;k++) {
				M[i][j][k] = 0;
			}
		}
	}    
//计算anchor i对应的E1
    for(i = 0; i< unknow_size;i++){
        for(j=0;j<sample_size;j++){
            e1[i][j] = E1(unknow[i], sample[j], p_index);
        }
    }

    for(y = 0; y < unknow_size; y++){
        cout << "Iterator Time:" << y <<endl;
        for(i = 0; i < unknow_size; i++){
            //总结 anchor i 所有 neighbor 发送的 信息
            initialize(neighM_sum, sample_size);
            for( k = 0; k<unknow[i].neighbors.size();k++){
                neighM_sum = AddTwoVec(neighM_sum, M[unknow[i].neighbors[k]][i], sample_size);                
            }
        for(k = 0; k < unknow[i].neighbors.size() ; k++){
            int neighbor = unknow[i].neighbors[k];
            if(converged[i][neighbor])
                continue;
            neighM = MinTwoVec(neighM_sum, M[neighbor][i], sample_size);
            neighM = AddTwoVec(neighM, e1[i], sample_size);
			//for (j = 0;j < sample_size;j++) {
			//	cout << neighM_sum[j] << " ";
			//}
			//cout << endl;

			//compute E2 of anchor i and anchor neighbor

			for (j = 0;j < sample_size;j++) {
				tmp_min = INFINITY;
				for (x = 0;x < sample_size;x++) {
					float tmp = E2(unknow[i], unknow[neighbor], sample[x], sample[j], p_index);
					if (tmp < 0) {
						cout <<"index:"<< i << " " << neighbor << " ";
						cout << this->lines[p_index][unknow[i].center_point] << " " 
							<< this->lines[p_index][unknow[neighbor].center_point] << endl;
						tmp = 0;
						break;
					}
					//cout << "E2:" << tmp << endl;
					tmp_min = min(tmp_min, neighM[x] + tmp);
				}
				M_tmp[j] = tmp_min;
			}    
            tmp_flag = true;
            for(j = 0; j < sample_size; j++){
                if(M[i][neighbor][j] != M_tmp[j]){
                    tmp_flag = false;
                    M[i][neighbor][j] = M_tmp[j];
                }
            }
            if(tmp_flag){
                converged[i][neighbor] = true;
            }        
        }
    }
    }



    int index_tmp;
    double min_tmp;


    for(i = 0; i < unknow_size; i++){
        initialize(neighM_sum, sample_size);
        neighM_sum = AddTwoVec(neighM_sum, e1[i], sample_size);
		for (k = 0;k < unknow[i].neighbors.size();k++) {
			int neighbor = unknow[i].neighbors[k];
			neighM_sum = AddTwoVec(neighM_sum, M[neighbor][i], sample_size);
		}
		min_tmp = INFINITY;
		for (k = 0;k < sample_size;k++) {
			//cout << neighM_sum[k] << " ";
			if (neighM_sum[k] < min_tmp) {
				index_tmp = k;
				min_tmp = neighM_sum[k];
				//cout << min_tmp << endl;
			}
		} 

        sample_index.push_back(index_tmp);               
    }

	cout << sample_index.size() << endl;
	cout << "BP() done" << endl;
	return sample_index;    
}








void GStructure::DrawNewRegion(){
    int line_num = this->lines.size();
    vector<int> Ischain(line_num);
    fill(Ischain.begin(),Ischain.end(), 1);
    int num;
    int shift;

    for(int i = 0; i < line_num;  i++)
        InputNeighbor(i);
    for(int i =0; i < line_num ; i++){
        if(this->unknown_anchor[i].empty())
            continue;
		for (int j = i+1 ;j < line_num;j++) {
			if (IsIntersect(i, j)) {
				Ischain[i] = 0;
				Ischain[j] = 0;
				num = this->lines[i].size();
				shift = this->unknown_anchor[i].size();
				for (int k = 0;k < this->unknown_anchor[j].size();k++) {
                //line index的 横移    
					this->unknown_anchor[j][k].head_point += num;
					this->unknown_anchor[j][k].center_point += num;
					this->unknown_anchor[j][k].end_point += num;
					for (int x = 0;x < unknown_anchor[j][k].neighbors.size();x++)
						this->unknown_anchor[j][k].neighbors[x] += (int)shift;
					this->unknown_anchor[i].push_back(this->unknown_anchor[j][k]);
				}
				for (int k = 0;k < this->sample_anchor[j].size();k++) {
					this->sample_anchor[j][k].head_point += num;
					this->sample_anchor[j][k].center_point += num;
					this->sample_anchor[j][k].end_point += num;
					this->sample_anchor[i].push_back(this->sample_anchor[j][k]);
				}
				for (int k = 0;k < this->lines[j].size();k++) {
					this->lines[i].push_back(this->lines[j][k]);
				}
				this->unknown_anchor[j].clear();
				this->sample_anchor[j].clear();
				this->lines[j].clear();

			}
		}
    }
    cout<<"Adjusting" <<endl;
    _sleep(4000);

	// for (int i = 0;i < this->lines.size();i++) {
	// 	cout << this->unknown_anchor[i].size() << " " << this->sample_anchor[i].size() << endl;
	// 	cout << "point num:" << this->lines[i].size() << endl;
	// 	if (Ischain[i]/*union_set[i] == -1 && union_num[i] == 1*/) {
	// 		//cout << i << endl;
	// 		CompleteLine(this->sample_anchor[i], this->unknown_anchor[i], i, 0);
	// 	}		
	// 	else if (!this->unknown_anchor[i].empty()) {
	// 		////检查邻居点
	// 		//for (int j = 0;j < this->unknown_anchor[i].size();j++) {

			
	// 		//CompleteLine(this->image_loss, this->sample_anchor[i], this->unknown_anchor[i], this->lines[i]);
	// 		CompleteLine(this->sample_anchor[i], this->unknown_anchor[i], i, 1);
	// 	}
	// }

	// namedWindow("new pic", 1);
    
	// imshow("new pic", this->image_loss);
    // destroyWindow("new pic");

	// waitKey(0);   
    Mat result = imread("./full/debug/result.png");
    cout<<"show result"<<endl;
    result.convertTo(result, CV_8UC3);
    namedWindow("result", 1);
    imwrite("result.png",result);
    imshow("result",result); 
    waitKey(0); 
}

void GStructure::InputNeighbor(int index){
    for (int i = 0;i < this->unknown_anchor[index].size();i++) {
        if (i - 1 >= 0)
            unknown_anchor[index][i].neighbors.push_back(i - 1);
        if (i + 1 < unknown_anchor[index].size())
            unknown_anchor[index][i].neighbors.push_back(i + 1);
	}
}


void GStructure::CompleteLine(vector<Anchor > sample, vector<Anchor> unknow, int p_index, bool flag){
    vector<int> sample_index;
	if(!flag)
		sample_index = DP(sample, unknow, p_index);
	else
		sample_index = BP(sample, unknow, p_index);
	Mat background = this->image_loss;
	int index;
	for (int i = 0;i < sample_index.size();i++) {
		index = sample_index[i];
		Mat patch = GetPatch(sample[index], p_index);
		CopyPic(unknow[i], patch, background, p_index);
	}
}
Anchor::Anchor(int head, int center, int end, PointType t):
    head_point(head), center_point(center), end_point(end), type(t){}

Anchor::~Anchor(){
    //cout<<"release"<<endl;
}

int main(){
    Mat input = imread("input.png");
    input.convertTo(input, CV_8UC3);
    GStructure gs(input);

    gs.AddAnchor();
    gs.DrawNewRegion();
    





    // Mat result = imread("sp2.png");
    // result.convertTo(result,CV_8UC3);
    // namedWindow("after_struct",1);
    // imshow("after_struct",result);

    return 0;
}
