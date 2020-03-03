#include <cmath>
#include <vector>
#include <iostream>

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

double m_modulus(vector<double> vec, bool square){
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

vector<int> vecAdd(vector<int> v1, vector<int> v2){
	vector<int> out(3);
	for(int i=0; i<3; i++){
		out[i] = v1[i] + v2[i];
	}
	return out;
}

vector<double> scalMult(double scal, vector<double> v){
	vector<double> out(3);
	for(int i=0; i<3; i++){
		out[i] = scal*v[i];
	}
	return out;
}

vector<double> compMult(vector<double> v1, vector<double> v2) {
    double real = v1[0]*v2[0] - v1[1]*v2[1];
    double imag = v1[1]*v2[0] + v1[0]*v2[1];
    return {real, imag};
}


void printVec(vector<double> vec){
    cout << "{";
    for(int i=0; i<3; i++){
        cout << vec[i];
        if(i != 2) {cout << ", ";}
    }
    cout << "}";
}