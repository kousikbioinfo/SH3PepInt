#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <vector>
#include <stdexcept>
#include <map>
#include <set>
#include <queue>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cassert>
#include <ctime>
#include <iomanip>
#include <limits>

using namespace std;

template <class OutType, class InType>
OutType stream_cast(const InType & t)
{
 stringstream ss;
 ss << t; // first insert value to stream
 OutType result; // value will be converted to OutType
 ss >> result; // write value to result
 return result;
}

#endif
