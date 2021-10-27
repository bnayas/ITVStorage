#include "pch.h"
#include "Utility.h"

#define _CRT_SECURE_NO_WARNINGS

unsigned int find_in_sorted(double n, std::vector<double> cdf, int size)
{
	auto itr = cdf.begin();
	for (int i = 0; i < size; i++ )
	{
		if (*itr > n)
			return (i);
		++itr;
	}
	return(size - 1);
}


std::string round(float number, int digits)
{
    double val_curr = number + 0.00000000000001;
    std::string str = "";
    int digit = 0;
    double value = 0;
    if (number == 0)
    {
        str += "0.";
        for (int i = 0; i < digits; i++)
            str += "0";
    }
    else
    {
        str = std::to_string(int(number)) + ".";
        value = int(number);
        double eps = 0.01 / pow(10, digits);
        double val_curr = (number - value) * 10;
        for (int i = 0; i < digits; i++)
        {
            digit = floor(val_curr + eps);
                if (value == 0)
                    value = digit;
                str += std::to_string(int(digit));
                val_curr = (val_curr - digit + eps) * 10;
            
        }
        if (value == 0)
            val_curr = number * pow(10, digits + 1) ;
        while (value == 0)
        {
            digit = floor(val_curr+0.00001);
            value = digit;
            str += std::to_string(int(digit));
            val_curr = (val_curr + 0.00001 - digit) * 10;
        }
    }
	return(str);
}


bool fileExists(const char* fileName)
{
	std::ifstream test(fileName);
	return (test) ? true : false;
}

param::param(int argc, char* argv[])
{
    for (int i = 1; i < argc; i++)
    {
        char* str_i = argv[i];
        char* pch = std::strtok(str_i, "=");
        if (std::strcmp(str_i, "-continue") == 0)
        {
            ContinueOld = true;
        }
        else {
            if (std::strcmp(pch, "N") == 0)
            {
                pch = std::strtok(NULL, "=");
                int value = std::strtod(pch, NULL);
                N = value;
            }
            else if (std::strcmp(str_i, "gamma") == 0 || 
                std::strcmp(str_i, "gama") == 0)
            {
                pch = std::strtok(NULL, "=");
                float value = std::strtod(pch, NULL);
                gama0 = value;
            }
            else if (std::strcmp(str_i, "runs") == 0)
            {
                pch = std::strtok(NULL, "=");
                float value = std::strtod(pch, NULL);
                runs = value;
            }
            else if (std::strcmp(str_i, "print_each") == 0)
            {
                pch = std::strtok(NULL, "=");
                float value = std::strtod(pch, NULL);
                print_each = value;
            }
            if (std::strcmp(pch, "delta") == 0)
            {
                pch = std::strtok(NULL, "=");
                float value = std::strtod(pch, NULL);
                delta = value;
            }
            if (std::strcmp(pch, "nu") == 0)
            {
                pch = strtok(NULL, "=");
                float value = std::strtod(pch, NULL);
                nu = value;
            }

        }
    }

}



/*int x = size >> i;
auto itr = cdf.begin() + x;
while (x > 0 && x < size)
{
	i++;
	if (n < *itr)
	{
		if (n > * (itr - 1))
			return x;
		else
		{
			x -= (size >> i) + 1;
			itr -= (size >> i) + 1;
		}
	}
	else
	{
		if (n < *(itr + 1))
			return (x + 1);
		else
		{
			x += (size >> i) + 1;
			itr += (size >> i) + 1;
		}
	}
}
if (x < size)
{
	if (n < *cdf.begin())
		return 0;
	else
		return 1;
}
else
	return(size - 1);*/
