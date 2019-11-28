#ifndef SAVE
#define SAVE

#include "bodies.h"

class bodyReader {
    friend std::ostream &operator<<(std::ostream &out, const body& obj);
    friend std::istream &operator>>(std::istream &out, const body& obj);
};

#endif //SAVE