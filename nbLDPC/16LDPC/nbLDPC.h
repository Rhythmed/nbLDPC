#ifndef _16_NON_BINARY_LOW_DENSITY_PARITY_CHECK_CODE_
#define _16_NON_BINARY_LOW_DENSITY_PARITY_CHECK_CODE_
#include <math.h>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
using namespace std;


typedef struct{
	int eb_;
	int ef_;
	int es_;
	int mj_;
	int frame_;
	double FER_;
	double UER_;
	double SER_;
	double BER_;
	double SNR_;
}SCTData;
typedef struct{
	vector<vector<int> > i_;
	vector<vector<int> > j_;
	vector<vector<int> > index_;
	vector<vector<int> > elemt_;
}SCTH;
class nbLDPC{
	private:
		int M_;
		int N_;
		int K_;
		int q_;
		int wR_;
		int wC_;
		int maxIter_;
		string path_;
	private:
		bool check_;
		vector<int> c_;
		vector<int> u_;
		vector<int> r_;
		vector<double> y_;
	private:
		SCTH sctH;
		vector<int> Exch_;
		vector<vector<int> > H_;
		vector<vector<int> > P_;
		vector<vector<double> > pro_;
	public:
		bool init();
		vector<int> test(double SNR);
		SCTData loop(int maxErr, double SNR);
	private:
		bool loadFile();
		bool calculateG();
	private:
		vector<int> encode();
		vector<int> decode();
		vector<int> message();
		vector<int> rearrange(vector<int> r);
		vector<int> calculateSx(vector<int> r);
		vector<int> compare(vector<int> c, vector<int> r);
		vector<double> tAWGN(double SNR);
		vector<vector<double> > calProbability(double SNR);
	public:
		 nbLDPC();
		 nbLDPC(int maxIter, string path);
		~nbLDPC();
};


static short addTable[16][16] =
{ { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
{ 1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14 },
{ 2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13 },
{ 3, 2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8, 15, 14, 13, 12 },
{ 4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11 },
{ 5, 4, 7, 6, 1, 0, 3, 2, 13, 12, 15, 14, 9, 8, 11, 10 },
{ 6, 7, 4, 5, 2, 3, 0, 1, 14, 15, 12, 13, 10, 11, 8, 9 },
{ 7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8 },
{ 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7 },
{ 9, 8, 11, 10, 13, 12, 15, 14, 1, 0, 3, 2, 5, 4, 7, 6 },
{ 10, 11, 8, 9, 14, 15, 12, 13, 2, 3, 0, 1, 6, 7, 4, 5 },
{ 11, 10, 9, 8, 15, 14, 13, 12, 3, 2, 1, 0, 7, 6, 5, 4 },
{ 12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3 },
{ 13, 12, 15, 14, 9, 8, 11, 10, 5, 4, 7, 6, 1, 0, 3, 2 },
{ 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1 },
{ 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 }
};
static short mulTable[16][16] =
{ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
{ 0, 2, 4, 6, 8, 10, 12, 14, 3, 1, 7, 5, 11, 9, 15, 13 },
{ 0, 3, 6, 5, 12, 15, 10, 9, 11, 8, 13, 14, 7, 4, 1, 2 },
{ 0, 4, 8, 12, 3, 7, 11, 15, 6, 2, 14, 10, 5, 1, 13, 9 },
{ 0, 5, 10, 15, 7, 2, 13, 8, 14, 11, 4, 1, 9, 12, 3, 6 },
{ 0, 6, 12, 10, 11, 13, 7, 1, 5, 3, 9, 15, 14, 8, 2, 4 },
{ 0, 7, 14, 9, 15, 8, 1, 6, 13, 10, 3, 4, 2, 5, 12, 11 },
{ 0, 8, 3, 11, 6, 14, 5, 13, 12, 4, 15, 7, 10, 2, 9, 1 },
{ 0, 9, 1, 8, 2, 11, 3, 10, 4, 13, 5, 12, 6, 15, 7, 14 },
{ 0, 10, 7, 13, 14, 4, 9, 3, 15, 5, 8, 2, 1, 11, 6, 12 },
{ 0, 11, 5, 14, 10, 1, 15, 4, 7, 12, 2, 9, 13, 6, 8, 3 },
{ 0, 12, 11, 7, 5, 9, 14, 2, 10, 6, 1, 13, 15, 3, 4, 8 },
{ 0, 13, 9, 4, 1, 12, 8, 5, 2, 15, 11, 6, 3, 14, 10, 7 },
{ 0, 14, 15, 1, 13, 3, 2, 12, 9, 7, 6, 8, 4, 10, 11, 5 },
{ 0, 15, 13, 2, 9, 6, 4, 11, 1, 14, 12, 3, 8, 7, 5, 10 }
};
static short divTable[16][16] =
{ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
{ 0, 1, 9, 14, 13, 11, 7, 6, 15, 2, 12, 5, 10, 4, 3, 8 },
{ 0, 2, 1, 15, 9, 5, 14, 12, 13, 4, 11, 10, 7, 8, 6, 3 },
{ 0, 3, 8, 1, 4, 14, 9, 10, 2, 6, 7, 15, 13, 12, 5, 11 },
{ 0, 4, 2, 13, 1, 10, 15, 11, 9, 8, 5, 7, 14, 3, 12, 6 },
{ 0, 5, 11, 3, 12, 1, 8, 13, 6, 10, 9, 2, 4, 7, 15, 14 },
{ 0, 6, 3, 2, 8, 15, 1, 7, 4, 12, 14, 13, 9, 11, 10, 5 },
{ 0, 7, 10, 12, 5, 4, 6, 1, 11, 14, 2, 8, 3, 15, 9, 13 },
{ 0, 8, 4, 9, 2, 7, 13, 5, 1, 3, 10, 14, 15, 6, 11, 12 },
{ 0, 9, 13, 7, 15, 12, 10, 3, 14, 1, 6, 11, 5, 2, 8, 4 },
{ 0, 10, 5, 6, 11, 2, 3, 9, 12, 7, 1, 4, 8, 14, 13, 15 },
{ 0, 11, 12, 8, 6, 9, 4, 15, 3, 5, 13, 1, 2, 10, 14, 7 },
{ 0, 12, 6, 4, 3, 13, 2, 14, 8, 11, 15, 9, 1, 5, 7, 10 },
{ 0, 13, 15, 10, 14, 6, 5, 8, 7, 9, 3, 12, 11, 1, 4, 2 },
{ 0, 14, 7, 11, 10, 8, 12, 2, 5, 15, 4, 3, 6, 13, 1, 9 },
{ 0, 15, 14, 5, 7, 3, 11, 4, 10, 13, 8, 6, 12, 9, 2, 1 },
};


#endif