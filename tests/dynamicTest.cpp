#include <iostream>
#include <vector>

using namespace std;

struct body{
    vector<double> pos;
    explicit body(int);
};

body::body(int i){
    pos = {double (i), double (i)};
}

void dymTest(){
    cout << "ok" << endl;
    double* arr = new double[5];
    for(int i=0; i<5; i++){
        arr[i] = i;
    }
    for(int i=0; i<5; i++){
        cout << arr[i] << endl;
    }
    delete arr;
}

void vecTest(){
    vector<body> bodies;
    for(int i=0; i<10; i++){
        bodies.push_back(body(i));
    }
    for(int i=0; i<10; i++){
        cout << bodies[i].pos[0] << endl;
    }
}

int main(){
    //vecTest();
    double* a = new double;
    *a = 6;
    cout << *a << endl;
    delete a;
    if (*a) cout << "ok" << endl;
    cout << *a << endl;
}