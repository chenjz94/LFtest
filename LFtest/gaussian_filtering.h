#define PI 3.14159265358979323846



void createFilter(int X, int U, std::vector<int> sigma_set, std::vector<double> phi_set, std::vector< std::vector<cv::Mat> > &Kernels);
void create1stFilter(int X, int U, std::vector<int> sigma_set, std::vector<double> phi_set, std::vector< std::vector<cv::Mat> > &Kernels);