#include <fstream>
#include <iostream>
#include "save.h"

using namespace std;

std::ostream &operator<<(std::ostream &out, const body &obj) {
    out << obj.pos[0] << "\n" << obj.pos[1] << "\n" << obj.pos[2] << endl;
    return out;
}

std::istream &operator>>(std::istream &in, const body &obj) {
    return in;
}