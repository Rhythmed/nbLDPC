#define _NON_BINARY_LOW_DENSITY_PARITY_CHECK_CODE_SEPTERMBER_
#include "nbLDPC.h"
#include <ctime>
#include <string>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <algorithm>
using namespace std;


bool printData(string path, SCTData data);

int main(int argc, char** argv){
	srand(time(NULL));
	int maxErr = 10;
	int maxIter = 10;
	float SNR = 3.0;
	string path = "./data/H.txt";

	sscanf(argv[1], "%d", &maxErr);
	sscanf(argv[2], "%d", &maxIter);
	sscanf(argv[3], "%f", &SNR);
	nbLDPC ldpc(maxIter, path);
	if (!ldpc.init())
		return -1;
	SCTData data = ldpc.loop(maxErr, SNR);
	printData("./result/out.txt", data);

	return 0;
}


bool printData(string path, SCTData data){
	ofstream out(path.c_str(), ios::app);
	if (!out.is_open())
		return false;
	out << " SNR :" << "  " << data.SNR_ << endl;
	out << "   FER:" << "  " << data.FER_ << endl;
	out << "   UER:" << "  " << data.UER_ << endl;
	out << "   SER:" << "  " << data.SER_ << endl;
	out << "   BER:" << "  " << data.BER_ << endl;
	out << "   error bits:\t" << data.eb_ << endl;
	out << "   error symbol: " << data.es_ << endl;
	out << "   error frame:\t" << data.ef_ << endl;
	out << "   error misjud: " << data.mj_ << endl;
	out << "   total frame:\t" << data.frame_ << endl << endl;

	out.close();
	return true;
}
