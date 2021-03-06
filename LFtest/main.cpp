﻿#include <iostream>
#include <opencv2\opencv.hpp>
#include <limits>
#include "light_field_proc.h"
#include "gaussian_filtering.h"

#define SIGMA_MIN 1
#define SIGMA_MAX 101
#define SIGMA_STEP 5

//stilllife
//#define PHI_MIN (-0.3976*pi)
//#define PHI_MAX (0.3976*pi)
#define PHI_MIN (-0.4*PI)
#define PHI_MAX (0.4*PI)
#define PHI_STEP (0.0125*PI)

#define KERNEL_WIDTH 401
#define OCTAVE 1

////stilllife
//#define Baseline 153.6000
//#define focallength 9.3750
//#define deltaX 72

////mona
//#define Baseline 38.4000
//#define focallength 9.3750
//#define deltaX 18

////buddha
//#define Baseline 48
//#define focallength 9.3750
//#define deltaX 22.4999


//medieval
#define Baseline 38.4000
#define focallength 10.3198
#define deltaX 19.8140
using namespace std;

int main()
{
	int dims[5];
	int* LF = light_field_read_HCI("D:\\Users\\cjz\\Documents\\MATLAB\\LFToolbox0.4\\HCI\\medieval\\lf.h5", dims);
	//"D:\\Users\\cjz\\Documents\\projects\\LF\\test\\LFtest\\LFtest\\lf.h5"
	//"D:\\Users\\cjz\\Documents\\MATLAB\\LFToolbox0.4\\HCI\\monasRoom\\lf.h5"
	//"D:\\Users\\cjz\\Documents\\MATLAB\\LFToolbox0.4\\HCI\\buddha\\buddha.h5"
	//"D:\\Users\\cjz\\Documents\\MATLAB\\LFToolbox0.4\\HCI\\medieval\\lf.h5"
	int X = dims[2], Y = dims[3], S = dims[0], T = dims[1];
	cout << X << Y << S << T;
	int dims_GT[4];
	double* GT = light_field_read_GT("D:\\Users\\cjz\\Documents\\MATLAB\\LFToolbox0.4\\HCI\\medieval\\lf.h5", dims_GT);
	int* EPI;
	cv::Mat lf = light_field_to_mat(LF + 3 * X * Y * 4, Y, X);//row,col
	cv::imshow("LF", lf);
	//cv::waitKey();
	cv::Mat gt = GT_to_mat(GT, Y, X, 1);

	cv::imshow("GT", gt);
	//cv::waitKey();
	vector<int> sigma_set;
	vector<double> phi_set;
	vector<double> depth_set;


	vector< vector<cv::Mat> > Kernels;

	cv::Mat img, img_c, conv_res;
	vector<Peaks> peaks[3];
	//int a[2] = { 0, 1 };
	//cout << a[3];
	//if (a[1] != 1 || a[3]==2)
	//{
	//	cout << 'O';
	//}

	for (int i = SIGMA_MIN; i <= SIGMA_MAX; i += SIGMA_STEP)
	{
		sigma_set.push_back(i);
	}
	for (double i = PHI_MIN; i <= PHI_MAX; i += PHI_STEP)
	{
		phi_set.push_back(i);
		depth_set.push_back(tan(i));
	}

	int sigma_num = sigma_set.size();
	int phi_num = phi_set.size();

	Kernels = vector< vector<cv::Mat> >(sigma_set.size(), vector<cv::Mat>(phi_set.size()));
	for (int sigma = 0; sigma < sigma_set.size(); sigma++)
	{
		for (int phi = 0; phi < phi_set.size(); phi++)
		{
			Kernels[sigma][phi] = cv::Mat::zeros(S, KERNEL_WIDTH, CV_32FC1);
		}
	}
	createFilter(KERNEL_WIDTH, S, sigma_set, phi_set, Kernels);

	//for (int sigma = 0; sigma < sigma_num; sigma++)
	//{
	//	for (int phi = 0; phi < phi_num; phi++)
	//	{
	//		cv::Mat norm_ker = cv::Mat::zeros(S, KERNEL_WIDTH, CV_32FC1);
	//		cout << sigma << ' ' << phi << endl;
	//		cv::normalize(Kernels[sigma][phi], norm_ker, 1.0, 0.0, cv::NORM_MINMAX);//¹éÒ»µ½0~1Ö®¼ä
	//		cv::imshow("KERNEL", norm_ker);
	//		cv::waitKey(0);
	//	}
	//}



	//cv::filter2D(img, conv_res, CV_32F, Kernels[0][0], cv::Point(-1, -1), 0, cv::BORDER_ISOLATED);
	//normalize_and_show("res", conv_res);
	//cout << conv_res.at<float>(0, 0);
	vector<vector<vector <float> > > Lisad2(sigma_num, vector<vector <float> >(phi_num, vector <float>(X, 0.)));
	vector<vector<int> > phi_assigned;
	vector<vector<float> > value_assigned;
	cv::Mat depth_map[3];
	depth_map[0] = cv::Mat::zeros(X,Y, CV_32FC1);// +depth_set[0] - 3;
	depth_map[1] = cv::Mat::zeros(X,Y, CV_32FC1);// +depth_set[0] - 3;
	depth_map[2] = cv::Mat::zeros(X,Y, CV_32FC1);// +depth_set[0] - 3;
	cv::Mat res_223(sigma_num, phi_num, CV_32FC1);
	cv::Mat r, g, b;
	cv::Mat out[] = { b, g, r };
	cv::Mat mask = cv::Mat::zeros(X, Y, CV_64FC1);
	for (int i = 0; i < OCTAVE; i++)
	{
		for (int y = 0; y < X; y++) //350 //445
		{
			//int offset = 0;
			cout << y << endl;
			EPI = light_field_extract_yt_slice(X, Y, S, T, LF, y, 4);
			img_c = light_field_to_mat(EPI, Y, S);
			//img_c = cv::imread("D:\\Users\\cjz\\Documents\\projects\\LF\\test\\LFtest\\LFtest\\epi.jpg");
			//cv::cvtColor(img_c, img, CV_BGR2GRAY);
			//cv::imshow("epi" ,img_c);
			cv::cvtColor(img_c, img_c, CV_BGR2HSV);
			//cv::waitKey();
			cv::split(img_c, out);
			//cv::imshow("hsv", out[0]);
			//cv::waitKey();

			for (int ch = 2; ch < 3; ch++)
			{
				for (int sigma = 0; sigma < sigma_num; sigma++)
				{
					for (int phi = 0; phi < phi_num; phi++)
					{
						//cout << y << ' ' << sigma << ' ' << phi << endl;
						//cv::filter2D(img, conv_res, CV_32F, Kernels[sigma][phi], cv::Point(-1, -1), 0, cv::BORDER_ISOLATED);
						cv::filter2D(out[ch], conv_res, CV_32F, Kernels[sigma][phi], cv::Point(-1, -1), 0, cv::BORDER_ISOLATED);

						float *data = conv_res.ptr<float>(S / 2);
						//res_223.at<float>(sigma, phi) = data[171];
						Lisad2[sigma][phi].swap(vector<float>(data, data + conv_res.cols));
					}
				}
				vector<Peaks> peak = light_field_Lisad2_peaks(Lisad2, sigma_num, phi_num, Y);

				light_field_assign_phi(peak, phi_assigned, value_assigned, sigma_set, Y);

				light_field_phi2depth(phi_assigned, value_assigned, depth_set, Y, y, depth_map[ch], mask);

			}
		}
	}
	//int offset =   0;
	//for (int row = 0; row < img.rows; row++)  
	//{

	//}
	cv::Mat error[3];
	cv::Mat errorpercent[3];
	cv::Mat red[3];
	cv::Mat anchor[3];
	for (int ch = 2; ch < 3; ch++)
	{
		depth_map[ch] = (Baseline * focallength) / (depth_map[ch] + deltaX);
		//error[ch] = cv::Mat::zeros(X, Y, CV_64FC1);
		//errorpercent[ch] = cv::Mat::zeros(Y, X, CV_64FC1);
		depth_map[ch].convertTo(depth_map[ch], CV_64FC1);
		error[ch] = depth_map[ch] - gt;
		errorpercent[ch] = error[ch] / gt;
	}
	for (int ch = 2; ch < 3; ch++)
	{
		red[ch] = lf.clone();
		anchor[ch] = cv::Mat::zeros(lf.size(), CV_8UC1);
		for (int row = 0; row < lf.rows; row++)
		{
			for (int col = 0; col < lf.cols; col++)
			{
				if (errorpercent[ch].at<double>(row, col)>0.005)
				{
					red[ch].at<cv::Vec3b>(row, col) = cv::Vec3b(0, 0, 255);
				}
				else
				{
					anchor[ch].at<uchar>(row, col) = 1;
				}
			}
		}
	}

	cv::FileStorage fsFeature("./depth_map.xml", cv::FileStorage::WRITE);
	//fsFeature << "depth_map0" << depth_map[0];
	//fsFeature << "depth_map1" << depth_map[1];
	fsFeature << "depth_map2" << depth_map[2];

	//cv::FileStorage fsAchor("./anchor.xml", cv::FileStorage::WRITE);
	//fsFeature << "anchor0" << anchor[0];
	//fsFeature << "anchor1" << anchor[1];
	fsFeature << "anchor2" << anchor[2];
	fsFeature.release();
	//for (int ch = 0; ch < 3; ch++)
	//{
	//	red[ch] = lf.clone();

	//}

	//cv::imshow("epi", img);
	//cout << Kernels[0][0].at<float>(0, 0);
	cv::Mat depth_map_norm;
	cv::normalize(depth_map[0], depth_map_norm, 1.0, 0.0, cv::NORM_MINMAX);//¹éÒ»µ½0~1Ö®¼ä
	cv::imshow("depth_map", depth_map_norm);
	cv::waitKey(0);


	system("pause");
}