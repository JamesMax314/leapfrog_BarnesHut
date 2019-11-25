#include <math.h>

double modulus(vector<double>& v1, vector<double>& v1, bool square = false){
	double out = 0;
	for(int i=0; i<3; i++){
		out += v1[i]*v1[i] + v2[i]*v2[i];
	};
	if(square == false){
		return pow(out, 0.5);
	} else{
		return out;
	};
};

vector<double> vecAdd(vector<double> v1, vector<double> v2){
	vector<double> out;
	for(int i=0; i<3; i++){
		out[i] = v1[i] + v2[i];
	};
	return out;
};

vector<double> scalMult(double scal, vector<double> v){
	vector<double> out;
	for(int i=0; i<3; i++){
		out[i] = scal*v[i];
	};
	return out;
};