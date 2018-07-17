#include <string>
#include <iostream>
#include <opencv2\opencv.hpp>
#include "H5Cpp.h"
#include "datastruct.h"
//#include "light_field_proc.h"
using namespace std;
using namespace H5;

//#define FILE            "lf.h5"
#define DATASET         "LF"
#define GROUNDTRUTH     "GT_DEPTH"

//read HCI *.h5 light field image
int* light_field_read_HCI(const char * path, int* dims)
{
	
	H5File file(path, H5F_ACC_RDONLY);
	DataSet dataset = file.openDataSet(DATASET);
	H5T_class_t type_class = dataset.getTypeClass();
	/*              
	* Get filespace for rank and dimension
	*/
	DataSpace filespace = dataset.getSpace();
	/*
	* Get number of dimensions in the file dataspace
	*/
	int rank = filespace.getSimpleExtentNdims();
	/*
	* Get and print the dimension sizes of the file dataspace
	*/
	hsize_t _dims[5];    // dataset dimensions
	rank = filespace.getSimpleExtentDims(_dims);
	dims[0] = _dims[0];
	dims[1] = _dims[1];
	dims[2] = _dims[2];
	dims[3] = _dims[3];
	dims[4] = _dims[4];
	cout << "dataset rank = " << rank << ", dimensions "
		<< (unsigned long)(dims[0]) << " x "
		<< (unsigned long)(dims[1]) << " x "
		<< (unsigned long)(dims[2]) << " x "
		<< (unsigned long)(dims[3]) << " x "
		<< (unsigned long)(dims[4]) << endl;
	DataSpace mspace1(rank, _dims);
	int *data = new int[dims[0] * dims[1] * dims[2] * dims[3] * dims[4]];  
	dataset.read(data, PredType::NATIVE_INT, mspace1, filespace);
	
	int length_per_line = dims[1] * dims[2] * dims[3] * dims[4];
	int *data_out = new int[length_per_line];
	memcpy(data_out, data + length_per_line*4, sizeof(int)*length_per_line);
	free(data);
	return data_out;
}

//read HCI *.h5 light field image
double* light_field_read_GT(const char * path, int* dims)
{

	H5File file(path, H5F_ACC_RDONLY);
	DataSet dataset = file.openDataSet(GROUNDTRUTH);
	H5T_class_t type_class = dataset.getTypeClass();
	/*
	* Get filespace for rank and dimension
	*/
	DataSpace filespace = dataset.getSpace();
	/*
	* Get number of dimensions in the file dataspace
	*/
	int rank = filespace.getSimpleExtentNdims();
	/*
	* Get and print the dimension sizes of the file dataspace
	*/
	hsize_t _dims[4];    // dataset dimensions
	rank = filespace.getSimpleExtentDims(_dims);
	dims[0] = _dims[0];
	dims[1] = _dims[1];
	dims[2] = _dims[2];
	dims[3] = _dims[3];

	cout << "GT rank = " << rank << ", dimensions "
		<< (unsigned long)(dims[0]) << " x "
		<< (unsigned long)(dims[1]) << " x "
		<< (unsigned long)(dims[2]) << " x "
		<< (unsigned long)(dims[3]) << endl;
	DataSpace mspace1(rank, _dims);
	double *data = new double[dims[0] * dims[1] * dims[2] * dims[3]];
	dataset.read(data, PredType::NATIVE_DOUBLE, mspace1, filespace);

	int length_per_line = dims[1] * dims[2] * dims[3];
	int length_per_view = dims[2] * dims[3];
	double *data_out = new double[length_per_view];
	memcpy(data_out, data + length_per_line * 4 + length_per_view * 4, sizeof(double)*length_per_view);
	free(data);
	return data_out;
}

///////////////////////////////////////////////////////////////////////////////////////
// for one row LF data
int* light_field_downsample(int X, int Y, int S, int T, int* image)
{
	int view_offset = X * Y * 3;
	int newX = (X + 1) / 2;
	int newY = (Y + 1) / 2;
	long newlength = S * newX *newY * 3;
	int *image_downsampled = new int[newlength];
	long off_out = 0, off_in = 0;
	
	for (int s = 0; s < S; s++)
	{
		off_in = s * view_offset;
		for (int y = 0; y < Y; y+=2)
		{
			for (int x = 0; x < X; x+=2)
			{
				
				image_downsampled[off_out] = image[off_in];
				off_out ++;
				off_in += 2;
			}
			//cout << off_in << ' ' << off_out << endl;
		}
		
	}
	return image_downsampled;
}

int* light_field_extract_yt_slice(int X, int Y, int S, int T, int* image, int y, int t)
{
	int view_offset = X * Y * 3;
	int *EPI = new int[Y *S * 3];
	long off_in = y*Y * 3;
	long off_out = (S - 1)*Y * 3;
	for (int s = 0; s < S; s++)
	{
		memcpy(EPI + off_out, image + off_in, sizeof(int)*Y * 3);
		off_in += view_offset;
		off_out -= Y * 3;
		//cout << off_in << ' ' << off_out << endl;
	}

	return EPI; 
}



vector<Peaks> light_field_Lisad2_peaks(vector<vector<vector <float> > > Lisad2, int sigma_num, int phi_num, int X)
{
	vector<Peaks> peaks;
	
	for (int sigma = 1; sigma < sigma_num - 1; sigma++)
	{
		for (int phi = 1; phi < phi_num - 1; phi++)
		{
			for (int x = 1; x < X - 1; x++)
			{
				if (//Lisad2[sigma][phi][x] <= Lisad2[sigma - 1][phi - 1][x - 1] &&
					//Lisad2[sigma][phi][x] <= Lisad2[sigma - 1][phi - 1][x] &&
					//Lisad2[sigma][phi][x] <= Lisad2[sigma - 1][phi - 1][x + 1] &&//
					//Lisad2[sigma][phi][x] <= Lisad2[sigma][phi - 1][x - 1] &&
					Lisad2[sigma][phi][x] <= Lisad2[sigma][phi - 1][x] &&
					//Lisad2[sigma][phi][x] <= Lisad2[sigma][phi - 1][x + 1] &&//
					//Lisad2[sigma][phi][x] <= Lisad2[sigma + 1][phi - 1][x - 1] &&
					//Lisad2[sigma][phi][x] <= Lisad2[sigma + 1][phi - 1][x] &&
					//Lisad2[sigma][phi][x] <= Lisad2[sigma + 1][phi - 1][x + 1] &&//
					//Lisad2[sigma][phi][x] <= Lisad2[sigma - 1][phi][x - 1] &&
					Lisad2[sigma][phi][x] <= Lisad2[sigma - 1][phi][x] &&
					//Lisad2[sigma][phi][x] <= Lisad2[sigma - 1][phi][x + 1] &&//
					Lisad2[sigma][phi][x] <= Lisad2[sigma][phi][x - 1] &&
					Lisad2[sigma][phi][x] <= Lisad2[sigma][phi][x + 1] &&//
					//Lisad2[sigma][phi][x] <= Lisad2[sigma + 1][phi][x - 1] &&
					Lisad2[sigma][phi][x] <= Lisad2[sigma + 1][phi][x] &&
					//Lisad2[sigma][phi][x] <= Lisad2[sigma + 1][phi][x + 1] &&//
					//Lisad2[sigma][phi][x] <= Lisad2[sigma - 1][phi + 1][x - 1] &&
					//Lisad2[sigma][phi][x] <= Lisad2[sigma - 1][phi + 1][x] &&
					//Lisad2[sigma][phi][x] <= Lisad2[sigma - 1][phi + 1][x + 1] &&//
					//Lisad2[sigma][phi][x] <= Lisad2[sigma][phi + 1][x - 1] &&
					Lisad2[sigma][phi][x] <= Lisad2[sigma][phi + 1][x] //&&
					//Lisad2[sigma][phi][x] <= Lisad2[sigma][phi + 1][x + 1] &&//
					//Lisad2[sigma][phi][x] <= Lisad2[sigma + 1][phi + 1][x - 1] &&
					//Lisad2[sigma][phi][x] <= Lisad2[sigma + 1][phi + 1][x] &&
					//Lisad2[sigma][phi][x] <= Lisad2[sigma + 1][phi + 1][x + 1] //
					){
					Peaks peak = { phi, sigma, x, Lisad2[sigma][phi][x] };
					peaks.push_back(peak);
				}

				else if (
					//Lisad2[sigma][phi][x] >= Lisad2[sigma - 1][phi - 1][x - 1] &&
					//Lisad2[sigma][phi][x] >= Lisad2[sigma - 1][phi - 1][x] &&
					//Lisad2[sigma][phi][x] >= Lisad2[sigma - 1][phi - 1][x + 1] &&//
					//Lisad2[sigma][phi][x] >= Lisad2[sigma][phi - 1][x - 1] &&
					Lisad2[sigma][phi][x] >= Lisad2[sigma][phi - 1][x] &&
					//Lisad2[sigma][phi][x] >= Lisad2[sigma][phi - 1][x + 1] &&//
					//Lisad2[sigma][phi][x] >= Lisad2[sigma + 1][phi - 1][x - 1] &&
					//Lisad2[sigma][phi][x] >= Lisad2[sigma + 1][phi - 1][x] &&
					//Lisad2[sigma][phi][x] >= Lisad2[sigma + 1][phi - 1][x + 1] &&//
					//Lisad2[sigma][phi][x] >= Lisad2[sigma - 1][phi][x - 1] &&
					Lisad2[sigma][phi][x] >= Lisad2[sigma - 1][phi][x] &&
					//Lisad2[sigma][phi][x] >= Lisad2[sigma - 1][phi][x + 1] &&//
					Lisad2[sigma][phi][x] >= Lisad2[sigma][phi][x - 1] &&
					Lisad2[sigma][phi][x] >= Lisad2[sigma][phi][x + 1] &&//
					//Lisad2[sigma][phi][x] >= Lisad2[sigma + 1][phi][x - 1] &&
					Lisad2[sigma][phi][x] >= Lisad2[sigma + 1][phi][x] &&
					//Lisad2[sigma][phi][x] >= Lisad2[sigma + 1][phi][x + 1] &&//
					//Lisad2[sigma][phi][x] >= Lisad2[sigma - 1][phi + 1][x - 1] &&
					//Lisad2[sigma][phi][x] >= Lisad2[sigma - 1][phi + 1][x] &&
					//Lisad2[sigma][phi][x] >= Lisad2[sigma - 1][phi + 1][x + 1] &&//
					//Lisad2[sigma][phi][x] >= Lisad2[sigma][phi + 1][x - 1] &&
					Lisad2[sigma][phi][x] >= Lisad2[sigma][phi + 1][x]// &&
					//Lisad2[sigma][phi][x] >= Lisad2[sigma][phi + 1][x + 1] &&//
					//Lisad2[sigma][phi][x] >= Lisad2[sigma + 1][phi + 1][x - 1] &&
					//Lisad2[sigma][phi][x] >= Lisad2[sigma + 1][phi + 1][x] &&
					//Lisad2[sigma][phi][x] >= Lisad2[sigma + 1][phi + 1][x + 1] //
					)
				{
					Peaks peak = { phi, sigma, x, Lisad2[sigma][phi][x] };
					peaks.push_back(peak);
				}
			}
		}
		
	}
	return peaks;
}

void light_field_assign_phi(vector<Peaks> peaks, vector<vector<int> > &phi, vector<vector<float> > &value, vector<int> sigma_set, int X)
{
	int sigma, x;
	int left, right;
	phi = vector<vector<int> >(X);
	value = vector<vector<float> >(X);
	for (int i = 0; i < peaks.size(); i++)
	{
		sigma = sigma_set[peaks[i].sigma];
		x = peaks[i].x;
		left = x - sigma > 0 ? x - sigma : 0;
		right = x + sigma < X ? x + sigma : X;
		for (int j = left; j < right; j++)
		{
			phi[j].push_back(peaks[i].phi);
			value[j].push_back(peaks[i].peaks);
		}
	}
}

void light_field_assign_phi_1st(vector<Peaks> peaks, vector<vector<int> > &phi, vector<vector<float> > &value, vector<int> sigma_set, int X)
{
	int sigma, x;
	int left, right;
	phi = vector<vector<int> >(X);
	value = vector<vector<float> >(X);
	for (int i = 0; i < peaks.size(); i++)
	{
		sigma = sigma_set[peaks[i].sigma];
		x = peaks[i].x;
		phi[x].push_back(peaks[i].phi);
		value[x].push_back(peaks[i].peaks);
	}
}


void light_field_phi2depth(vector<vector<int> > phi, vector<vector<float> > value, vector<double> depth_set, int X, int y, cv::Mat &depth_map, cv::Mat &mask)
{
	vector<float>::iterator max_value;
	int phi_index;
	for (int x = 0; x < X; x++)
	{
		if (!phi[x].empty())
		{
			max_value = std::max_element(std::begin(value[x]), std::end(value[x]));
			phi_index = phi[x][max_value - std::begin(value[x])];
			depth_map.at<float>(y, x) = depth_set[phi_index];
			mask.at<double>(y, x) = 1;
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////


// //for complete LF data
//int* light_field_extract_yt_slice(int X, int Y, int S, int T, int* image, int y, int t)
//{
//	int view_offset = X * Y * 3;
//	int *EPI = new int[ X *S *3];
//	long off_in = y*X*3 + t*S*view_offset;
//	long off_out = 0;
//	for (int s = 0; s < S; s++)
//	{
//
//		memcpy(EPI + off_out, image + off_in, sizeof(int)*X * 3);
//		off_in += view_offset;
//		off_out += X * 3; 
//		//cout << off_in << ' ' << off_out << endl;
//	}
//
//	return EPI;
//}


cv::Mat light_field_to_mat(int *image, int W, int H, bool gray = false)
{
	cv::Mat img;
	if (gray) img = cv::Mat(H, W, CV_8UC1);
	else img = cv::Mat(H, W, CV_8UC3);
	int offset = 0;
	for (int row = 0; row < img.rows; row++)
	{
		for (int col = 0; col < img.cols; col++)
		{
			if (gray)img.at<uchar>(row, col) = image[offset];//opencv:BGR; H5file:RGB
			else img.at<cv::Vec3b>(row, col) = cv::Vec3b(image[offset * 3 + 2], image[offset * 3 + 1], image[offset * 3]);//opencv:BGR; H5file:RGB
			offset++;
		}
	}
	return img;
	//cv::resize(img, img, cv::Size( W,H*10));
	//cv::imshow("LF", img2);
	//cv::waitKey(0);
}

cv::Mat GT_to_mat(double *image, int W, int H, bool gray = false)
{
	cv::Mat img;
	if (gray) img = cv::Mat(H, W, CV_64FC1);
	else img = cv::Mat(H, W, CV_64FC3);
	int offset = 0;
	for (int row = 0; row < img.rows; row++)
	{
		for (int col = 0; col < img.cols; col++)
		{
			if (gray)img.at<double>(row, col) = image[offset];//opencv:BGR; H5file:RGB
			else img.at<cv::Vec3d>(row, col) = cv::Vec3d(image[offset * 3 + 2], image[offset * 3 + 1], image[offset * 3]);//opencv:BGR; H5file:RGB
			offset++;
		}
	}
	return img;
	//cv::resize(img, img, cv::Size( W,H*10));
	//cv::imshow("LF", img2);
	//cv::waitKey(0);
}

void normalize_and_show(const char* title, cv::Mat img)
{
	cv::Mat norm = cv::Mat::zeros(img.size(), CV_32FC1);
	cv::normalize(img, norm, 1.0, 0.0, cv::NORM_MINMAX);//归一到0~1之间
	cv::imshow(title, norm);
	cv::waitKey();
	norm.release();
}




///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
//test

void mark_peak_on_EPI(cv::Mat EPI, float* res, int res_size)
{
	for (int i = 1; i < res_size - 1; i++)
	{
		//cout << res[i - 1] << ' ' << res[i] << ' ' << res[i + 1] << endl;
		if ((res[i] <= res[i - 1] && res[i] <= res[i + 1]) || (res[i] >= res[i - 1] && res[i] >= res[i + 1]))
		{
			cv::circle(EPI, cv::Point(i, 4), 4, cv::Scalar(0,0,255), 1, 8, 0);
		}
	}
}