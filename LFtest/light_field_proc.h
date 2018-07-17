#pragma once 
#include <string>
#include <opencv2\opencv.hpp>
#include "datastruct.h"
int* light_field_read_HCI(const char* path, int* dims);
double* light_field_read_GT(const char * path, int* dims);

int* light_field_extract_yt_slice(int X, int Y, int S, int T, int* image, int y, int t);

std::vector<Peaks> light_field_Lisad2_peaks(std::vector< std::vector< std::vector <float> > > Lisad2, int sigma_num, int phi_num, int X);
void light_field_assign_phi(std::vector<Peaks> peaks, std::vector<std::vector<int> > &phi, std::vector<std::vector<float> > &value,  std::vector<int> sigma_set, int X);
void light_field_phi2depth(std::vector<std::vector<int> > phi, std::vector<std::vector<float> > value, std::vector<double> depth_set, int X, int y, cv::Mat &depth_map, cv::Mat &mask);
void light_field_assign_phi_1st(std::vector<Peaks> peaks, std::vector<std::vector<int> > &phi, std::vector<std::vector<float> > &value, std::vector<int> sigma_set, int X);

cv::Mat light_field_to_mat(int *image, int W, int H, bool gray = false);
cv::Mat GT_to_mat(double *image, int W, int H, bool gray = false);
void normalize_and_show(const char* title, cv::Mat img);
void mark_peak_on_EPI(cv::Mat EPI, float* res, int res_size);