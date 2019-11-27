#include <cmath>
#include <vector>

using namespace std;

double dot(vector<double>& v1, vector<double>& v2, bool square){
	double out = 0;
	for(int i=0; i<3; i++){
		out += v1[i]*v2[i];
	}
	if(!square){
		return pow(out, 0.5);
	} else{
		return out;
	}
}

double modulus(vector<double> vec, bool square){
    double out = dot(vec, vec, square);
    return out;
}

vector<double> vecAdd(vector<double> v1, vector<double> v2){
	vector<double> out(3);
	for(int i=0; i<3; i++){
		out[i] = v1[i] + v2[i];
	}
	return out;
}

vector<double> scalMult(double scal, vector<double> &v){
	vector<double> out(30);
	for(int i=0; i<3; i++){
		out[i] = scal*v[i];
	}
	return out;
}