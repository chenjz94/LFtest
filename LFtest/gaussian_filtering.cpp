#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <opencv2\opencv.hpp>
#include "gaussian_filtering.h"

using namespace std;

void createFilter(int X, int U, vector<int> sigma_set, vector<double> phi_set, vector< vector<cv::Mat> > &Kernels)
{
	int xcenter = X / 2;
	int ucenter = U / 2;
	double a, b, c;
	int i_sigma = 0, i_phi = 0;
	for (auto sigma : sigma_set)
	{
		i_phi = 0;
		
		for (auto phi : phi_set)
		{
			for (int u = 0; u < U; u++)
			{
				for (int x = 0; x < X; x++)
				{
					a = exp(-pow(((x - xcenter) + (u - ucenter)*tan(phi)), 2) / (2 * pow(sigma, 2)));
					b = ((x - xcenter) + (u - ucenter)*tan(phi) - sigma) * ((x - xcenter) + (u - ucenter)*tan(phi) + sigma);
					c = pow(sigma, 5)*sqrt(2 * PI);
					Kernels[i_sigma][i_phi].at<float>(u,x) = a * b / c;
				}
			}
			Kernels[i_sigma][i_phi] *= sigma*sigma;
			i_phi++;
		}
		i_sigma++;
	}
 //   for (int sigma = 0; sigma < sigma_set.size(); sigma++)
	//{
	//	for (int phi = 0; phi < phi_set.size(); phi++)
	//	{
	//		cv::Mat norm_ker = cv::Mat::zeros(U, X, CV_32FC1);
	//		cv::normalize(Kernels[sigma][phi], norm_ker, 1.0, 0.0, cv::NORM_MINMAX);//归一到0~1之间
	//		cv::imshow("KERNEL", norm_ker);
	//		cv::waitKey(0);
	//	}
	//}

}

void createFilter_vertical(int X, int U, vector<int> sigma_set, vector<double> phi_set, vector< vector<cv::Mat> > &Kernels)
{
	int xcenter = X / 2;
	int ucenter = U / 2;
	double a, b, c;
	int i_sigma = 0, i_phi = 0;
	for (auto sigma : sigma_set)
	{
		i_phi = 0;

		for (auto phi : phi_set)
		{
			for (int u = 0; u < U; u++)
			{
				for (int x = 0; x < X; x++)
				{
					a = exp(-pow(((x - xcenter) + (u - ucenter)*tan(phi)), 2) / (2 * pow(sigma, 2)));
					b = ((x - xcenter) + (u - ucenter)*tan(phi) - sigma) * ((x - xcenter) + (u - ucenter)*tan(phi) + sigma);
					c = pow(sigma, 5)*sqrt(2 * PI);
					Kernels[i_sigma][i_phi].at<float>(x, u) = a * b / c;
				}
			}
			Kernels[i_sigma][i_phi] *= sigma*sigma;
			i_phi++;
		}
		i_sigma++;
	}
	//   for (int sigma = 0; sigma < sigma_set.size(); sigma++)
	//{
	//	for (int phi = 0; phi < phi_set.size(); phi++)
	//	{
	//		cv::Mat norm_ker = cv::Mat::zeros(U, X, CV_32FC1);
	//		cv::normalize(Kernels[sigma][phi], norm_ker, 1.0, 0.0, cv::NORM_MINMAX);//归一到0~1之间
	//		cv::imshow("KERNEL", norm_ker);
	//		cv::waitKey(0);
	//	}
	//}

}

void create1stFilter(int X, int U, vector<int> sigma_set, vector<double> phi_set, vector< vector<cv::Mat> > &Kernels)
{
	int xcenter = X / 2;
	int ucenter = U / 2;
	double a, b, c;
	int i_sigma = 0, i_phi = 0;
	for (auto sigma : sigma_set)
	{
		i_phi = 0;

		for (auto phi : phi_set)
		{
			for (int u = 0; u < U; u++)
			{
				for (int x = 0; x < X; x++)
				{
					a = exp(-pow(((x - xcenter) + (u - ucenter)*tan(phi)), 2) / (2 * pow(sigma, 2)));
					b = (x - xcenter) + (u - ucenter)*tan(phi) ;
					c = pow(sigma, 3)*sqrt(2 * PI);
					Kernels[i_sigma][i_phi].at<float>(u, x) = a * b / c;
				}
			}
			Kernels[i_sigma][i_phi] *= sigma;
			i_phi++;
		}
		i_sigma++;
	}
	//   for (int sigma = 0; sigma < sigma_set.size(); sigma++)
	//{
	//	for (int phi = 0; phi < phi_set.size(); phi++)
	//	{
	//		cv::Mat norm_ker = cv::Mat::zeros(U, X, CV_32FC1);
	//		cv::normalize(Kernels[sigma][phi], norm_ker, 1.0, 0.0, cv::NORM_MINMAX);//归一到0~1之间
	//		cv::imshow("KERNEL", norm_ker);
	//		cv::waitKey(0);
	//	}
	//}

}

//int main()
//{
//	double gKernel[5][5];
//	createFilter(gKernel);
//	for (int i = 0; i < 5; ++i)
//	{
//		for (int j = 0; j < 5; ++j)
//			cout << gKernel[i][j] << "\t";
//		cout << endl;
//	}
//}