#pragma once

struct param
{
   param(int argc, char* argv[]);
   int N;
   float gama0;
   int runs;
   int print_each;
   float delta;
   float nu;
   bool  ContinueOld;
};

unsigned int find_in_sorted(double n, std::vector<double> cdf, int size);
std::string round(float number, int digits);
bool fileExists(const char* fileName);

