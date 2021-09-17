#include "rand.h"
#include "nbLDPC.h"


vector<int> nbLDPC::encode(){
	vector<int> s(M_);
	for (int i = 0; i < M_; i++){
		int tmp = 0;
		for (int j = 0; j < K_; j++){
			if (P_[i][j] && c_[j])
				tmp = addTable[tmp][mulTable[P_[i][j]][c_[j]]];
		}
		s[i] = tmp;
	}
	this->u_.clear();
	for (int i = 0; i < N_; i++){
		if (i < M_)
			this->u_.push_back(s[i]);
		else
			this->u_.push_back(c_[i - M_]);
	}
	for (int i = (M_ - 1); i >= 0; i--){
		if (Exch_[i] > -1){
			int tmp = this->u_[i];
			this->u_[i] = this->u_[Exch_[i]];
			this->u_[Exch_[i]] = tmp;
		}
	}
	return u_;
}

vector<int> nbLDPC::decode(){
	this->r_.clear();
	this->check_ = false;
	vector<int> ck(N_);
	vector<vector<vector<double> > > proTmp(wR_);
	vector<vector<vector<double> > > proV2C(wR_);
	vector<vector<vector<double> > > proC2V(wR_);
	for (int i = 0; i < wR_; i++){
		proTmp[i].resize(q_);
		proV2C[i].resize(q_);
		proC2V[i].resize(q_);
		for (int j = 0; j < q_; j++){
			proTmp[i][j].resize(M_);
			proV2C[i][j].resize(M_);
			proC2V[i][j].resize(M_);
		}
	}
	for (int i = 0; i < N_; i++){
		for (int j = 0; j < wC_; j++){
			for (int k = 0; k < q_; k++){
				proV2C[sctH.j_[i][j]][k][sctH.i_[i][j]] = pro_[k][i];
			}
		}
	}

	for (int iter = 0; iter < maxIter_; iter++){
		if (!iter){																// first iter
			for (int i = 0; i < N_; i++){
				for (int j = 0; j < wC_; j++){
					double tmp = 0;
					int i_ = sctH.i_[i][j];
					int j_ = sctH.j_[i][j];
					for (int k = 0; k < q_; k++){
						tmp += pro_[k][i];
						proV2C[j_][k][i_] = pro_[k][i];
					}
					for (int k = 0; k < q_; k++){
						proV2C[j_][k][i_] = proV2C[j_][k][i_] / tmp;
						proTmp[j_][k][i_] = proV2C[j_][k][i_];
					}
					for (int k = 1; k < q_; k++){
						proV2C[j_][mulTable[k][H_[i_][i]]][i_] = proTmp[j_][k][i_];
					}
				}
			}
		}
		else{
			for (int i = 0; i < N_; i++){
				for (int j = 0; j < wC_; j++){
					int i_ = sctH.i_[i][j];
					int j_ = sctH.j_[i][j];
					for (int k = 0; k < q_; k++){
						proV2C[j_][k][i_] = pro_[k][i];
					}
					for (int k = 0; k < wC_; k++){
						if (k != j){
							int i1_ = sctH.i_[i][k];
							int j1_ = sctH.j_[i][k];
							for (int l = 0; l < q_; l++){
								proV2C[j_][l][i_] = proV2C[j_][l][i_] * proC2V[j1_][l][i1_];
							}
						}
					}
					double tmp = 0.0;
					for (int k = 0; k < q_; k++){
						proTmp[j_][k][i_] = proV2C[j_][k][i_];
						tmp += proTmp[j_][k][i_];
					}
					for (int k = 1; k < q_; k++)
						proV2C[j_][mulTable[k][H_[i_][i]]][i_] = proTmp[j_][k][i_];
					for (int k = 0; k < q_; k++){
						proV2C[j_][k][i_] = proV2C[j_][k][i_] / tmp;
					}
				}
			}
			for (int i = 0; i < N_; i++){
				int i_ = sctH.i_[i][1];
				int j_ = sctH.j_[i][1];
				vector<double> buffer;
				for (int j = 0; j < q_; j++){
					double tmp = proTmp[j_][j][i_] * proC2V[j_][j][i_];
					buffer.push_back(tmp);
				}
				int tmp = max_element(buffer.begin(), buffer.end()) - buffer.begin();
				ck[i] = tmp;
			}
			bool error = false;
			vector<int> Sx = calculateSx(ck);
			for (int i = 0; i < Sx.size(); i++){
				if (Sx[i]){
					error = true;
					break;
				}
			}
			if (!error){
				this->check_ = true;
				break;
			}
		}

		for (int i = 0; i < M_; i++){
			for (int j = 0; j < wR_; j++){
				int c2v = 0;
				int i_ = sctH.index_[i][j];
				for (int k = 0; k < wR_; k++){
					if (k != j){
						if (!c2v){
							for (int l = 0; l < q_; l++)
								proC2V[j][l][i] = proV2C[k][l][i];
						}
						else{
							double tmp = 0;
							for (int l = 0; l < q_; l++){
								tmp += proV2C[k][l][i];
							}
							if (tmp != 0){
								vector<vector<double> > buffer(q_);
								for (int l = 0; l < q_; l++){
									for (int m = 0; m < q_; m++)
										buffer[l].push_back(proC2V[j][l][i] * proV2C[k][m][i]);
								}
								for (int l = 0; l < q_; l++)
									proC2V[j][l][i] = buffer[l][0];
								for (int l = 1; l < q_; l++){
									for (int m = 0; m < q_; m++){
										if (buffer[m][l] > proC2V[j][addTable[m][l]][i]){
											proC2V[j][addTable[m][l]][i] = buffer[m][l];
										}
									}
								}
							}
						}
						c2v++;
					}
				}
				double tmp = 0.0;
				for (int k = 0; k < q_; k++)
					tmp += proC2V[j][k][i];
				for (int k = 0; k < q_; k++){
					proC2V[j][k][i] = proC2V[j][k][i] / tmp;
					proTmp[j][k][i] = proC2V[j][k][i];
				}
				if (i_ != -1){
					for (int k = 1; k < q_; k++){
						int div = divTable[k][H_[i][i_]];
						proC2V[j][div][i] = proTmp[j][k][i];
					}
				}
			}
		}
	}
	this->r_ = ck;
	return ck;
}

vector<int> nbLDPC::message(){
	int tmp = 0;
	this->c_.clear();
	for (int i = 0; i < K_; i++){
		tmp = rand() % (q_);
		this->c_.push_back(tmp);
	}
	return c_;
}

vector<int> nbLDPC::test(double SNR){
	message();
	encode();
	vector<int> Sx = calculateSx(u_);
	tAWGN(SNR);
	calProbability(SNR);
	cout << "  raw data is :\n" << "  ";
	for (int i = 0; i < c_.size(); i++)
		cout << setw(3) << c_[i];
	cout << endl;

	cout << "  encode :\n" << "  ";
	for (int i = 0; i < u_.size(); i++)
		cout << setw(3) << u_[i];
	cout << endl << endl;


	cout << "  Syndrome :\n" << "  ";
	for (int i = 0; i < Sx.size(); i++)
		cout << setw(3) << Sx[i];
	cout << endl << endl;


	cout << "  sign :\n" << "  ";
	for (int i = 0; i < y_.size(); i++)
		cout << y_[i] << "  ";
	cout << endl << endl;


	vector<int> r = decode();
	cout << "  decode :\n" << "  ";
	for (int i = 0; i < r.size(); i++)
		cout << setw(3) << r[i];
	cout << endl << endl;




	vector<int> data = rearrange(r);
	cout << "  actual data :\n" << "  ";
	for (int i = 0; i < data.size(); i++)
		cout << setw(3) << data[i];
	cout << endl;
	vector<int> tmp = compare(c_, data);
	cout << "  error bits :  " << tmp[1] << endl;
	cout << "  error symbols :  " << tmp[0];


	return c_;
}

SCTData nbLDPC::loop(int maxErr, double SNR){
	double FER;
	double UER;
	double SER;
	double BER;

	int frame = 0;
	int errBits = 0;
	int errFrame = 0;
	int misjudge = 0;
	int errSymbol = 0;
	while (errFrame < maxErr){
		message();
		encode();
		tAWGN(SNR);
		calProbability(SNR);
		decode();
		if (check_){
			vector<int> r = rearrange(r_);
			vector<int> tmp = compare(c_, r);
			if (tmp[0]){
				misjudge++;
				errBits += tmp[1];
				errSymbol += tmp[0];
			}
		}
		else{
			errFrame++;
			vector<int> r = rearrange(r_);
			vector<int> tmp = compare(c_, r);
			errBits += tmp[1];
			errSymbol += tmp[0];
		}
		frame++;
	}
	int nb = log(double(q_)) / log(2.0);
	FER = double(errFrame) / double(frame);
	UER = double(misjudge) / double(frame);
	SER = double(errSymbol) / double((frame*M_));
	BER = double(errBits) / double((frame*M_*nb));

	SCTData out;
	out.eb_ = errBits;
	out.mj_ = misjudge;
	out.ef_ = errFrame;
	out.es_ = errSymbol;
	out.SNR_ = SNR;
	out.BER_ = BER;
	out.FER_ = FER;
	out.SER_ = SER;
	out.UER_ = UER;
	out.frame_ = frame;
	return out;
}

vector<int> nbLDPC::rearrange(vector<int> r){
	for (int i = 0; i < Exch_.size(); i++){
		if (Exch_[i] >= 0){
			int tmp = r[i];
			r[i] = r[Exch_[i]];
			r[Exch_[i]] = tmp;
		}
	}
	vector<int> c;
	for (int i = M_; i < N_; i++)
		c.push_back(r[i]);
	return c;
}

vector<int> nbLDPC::calculateSx(vector<int> r){
	vector<int> Sx;
	for (int i = 0; i < M_; i++) {
		int tmp = 0;
		for (int j = 0; j < N_; j++){
			if (r[j] && H_[i][j])
				tmp = addTable[mulTable[r[j]][H_[i][j]]][tmp];
		}
		Sx.push_back(tmp);
	}
	return Sx;
}

vector<int> nbLDPC::compare(vector<int> c, vector<int> r){
	int bit = 0;
	int symbol = 0;
	vector<int> dif;
	vector<int> out(2);
	for (int i = 0; i < c.size(); i++){
		if (addTable[c[i]][r[i]]){
			symbol++;
			dif.push_back(addTable[c[i]][r[i]]);
		}
	}
	out[0] = symbol;
	int nb = log(q_) / log(2);
	for (int i = 0; i < dif.size(); i++){
		int tmp = dif[i];
		for (int j = 0; j < nb; j++){
			if (tmp & 1)
				bit++;
			tmp = tmp >> 1;
		}
	}
	out[1] = bit;
	return out;
}

vector<double> nbLDPC::tAWGN(double SNR){
	int nbit = log(double(q_)) / log(2.0);
	vector<int> bitSm(nbit*N_);
	for (int i = 0; i < N_; i++){										// 转为比特流
		int tmp = this->u_[i];
		for (int j = (nbit - 1); j >= 0; j--){
			bitSm[nbit*i + j] = (tmp & 1);
			tmp = tmp >> 1;
		}
	}
	for (int i = 0; i < N_*nbit; i++){									// BPSK调制
		bitSm[i] = 1 - (2 * bitSm[i]);									// 0 -->  1
	}																	// 1 --> -1
	double rc = double(K_) / double(N_);
	double tEbNO = pow(10.0, (SNR / 10.0));
	double vNoise = 1.0 / (2.0*rc*tEbNO);
	vNoise = sqrt(vNoise);

	this->y_.clear();
	for (int i = 0; i < (nbit*N_); i++){
		double tmp = vNoise * GaussRand();
		y_.push_back(tmp + bitSm[i]);
	}
	return y_;
}

vector<vector<double> > nbLDPC::calProbability(double SNR){
	int nbit = log(double(q_)) / log(2.0);
	vector<double> rec;
	vector<vector<double> > pro(q_);
	for (int i = 0; i < q_; i++)
		pro[i].resize(N_);
	double rc = double(K_) / double(N_);
	double sigma = 1.0 / (sqrt(2 * rc) * pow(10.0, SNR / 20));
	double tSigma = exp(-2.0 / pow(sigma, 2.0));
	for (int i = 0; i < y_.size(); i++){
		double tmp = 1.0 / (1 + pow(tSigma, y_[i]));
		rec.push_back(tmp);
	}
	for (int i = 0; i < N_; i++){
		for (int j = 0; j < q_; j++){
			int sym = j;
			double tmp = 1.0;
			for (int k = (nbit - 1); k >= 0; k--){
				int bit = (sym & 1);
				sym = sym >> 1;
				if (bit)
					tmp *= (1 - rec[i*nbit + k]);
				else
					tmp *= rec[i*nbit + k];
			}
			pro[j][i] = tmp;
		}
	}
	this->pro_ = pro;
	return pro;
}

bool nbLDPC::init(){
	if (!this->loadFile())
		return false;
	if (!this->calculateG())
		return false;
	return true;
}

bool nbLDPC::loadFile(){
	ifstream chk(path_.c_str());
	ifstream che("./data/elm.txt");
	if (!chk.is_open()){
		cout << "Sorry, could not find file named \"" << path_ << "\"" << endl;
		return false;
	}
	chk >> this->N_;
	chk >> this->M_;
	chk >> this->wR_;
	this->wC_ = 3;
	this->K_ = N_ - M_;
	vector<vector<int> > H(M_);
	vector<vector<int> > I_(N_);
	vector<vector<int> > J_(N_);
	vector<vector<int> > elm(M_);
	vector<vector<int> > idx(M_);
	if (!che.is_open()){								// withou element file
		for (int i = 0; i < M_; i++){
			int tmp = 0;
			H[i].resize(N_);
			elm[i].resize(wR_);
			idx[i].resize(wR_);
			for (int j = 0; j < wR_; j++){
				chk >> tmp;
				idx[i][j] = tmp - 1;
				if (tmp){
					elm[i][j] = rand() % (q_ - 1) + 1;
					H[i][(tmp - 1)] = elm[i][j];
				}
			}
		}

	}
	else{
		for (int i = 0; i < M_; i++){
			int tmp = 0;
			int emp = 0;
			H[i].resize(N_);
			elm[i].resize(wR_);
			idx[i].resize(wR_);
			for (int j = 0; j < wR_; j++){
				chk >> tmp;
				che >> emp;
				idx[i][j] = tmp - 1;
				if (tmp){
					elm[i][j] = emp;
					H[i][(tmp - 1)] = elm[i][j];
				}
			}
		}
	}
	for (int i = 0; i < N_; i++){
		for (int j = 0; j < M_; j++){
			if (H[j][i]){
				I_[i].push_back(j);
				for (int k = 0; k < wR_; k++){
					if (idx[j][k] == i){
						J_[i].push_back(k);
						break;
					}
				}
			}
		}
	}
	chk.close();
	che.close();
	this->H_ = H;
	this->sctH.i_ = I_;
	this->sctH.j_ = J_;
	this->sctH.elemt_ = elm;
	this->sctH.index_ = idx;
	return true;
}

bool nbLDPC::calculateG(){
	vector<int> Exch(M_, -1);
	vector<vector<int> > H;
	vector<vector<int> > P(M_);

	H = this->H_;
	for (int i = 0; i < this->M_; i++){
		if (!H[i][i]){
			int tmp = *max_element(H[i].begin(), H[i].end());
			if (!tmp)
				return false;
			for (int j = 0; j < N_; j++)
			if (H[i][j]){
				tmp = j;
				break;
			}
			vector<int> buffer;
			buffer.clear();
			for (int j = 0; j < M_; j++)
				buffer.push_back(H[j][i]);
			for (int j = 0; j < M_; j++){
				H[j][i] = H[j][tmp];
				H[j][tmp] = buffer[j];
			}
			Exch[i] = tmp;
		}
		int tmp = H[i][i];
		for (int j = 0; j < N_; j++)
			H[i][j] = divTable[H[i][j]][tmp];
		for (int j = 0; j < M_; j++){
			if ((j != i) && (H[j][i])){
				int tmp = H[j][i];
				for (int k = i; k < N_; k++)
					H[j][k] = addTable[mulTable[tmp][H[i][k]]][H[j][k]];
			}
		}
	}
	for (int i = 0; i < M_; i++){
		P[i].resize(K_);
		for (int j = 0; j < K_; j++)
			P[i][j] = H[i][j + M_];
	}
	this->P_ = P;
	this->Exch_ = Exch;
	return true;
}

nbLDPC::nbLDPC(){
	this->q_ = 16;
	this->path_ = "./data/H.txt";
	this->maxIter_ = 50;
	cout << "  Default path is \"" << path_ << "\" ,make sure that file exists. \n";
}

nbLDPC::nbLDPC(int maxIter, string path){
	this->q_ = 16;
	this->path_ = path;
	this->maxIter_ = maxIter;
}

nbLDPC::~nbLDPC()
{
}
