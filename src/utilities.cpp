#include "utilities.hpp"
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include <iostream>

//Returns an index using roulette selection (disproportionately selecting those with highest fitness).
unsigned int rouletteSelect(const std::vector<double>& _fitness, double _sum)
{
    if (_sum<0)
        return rouletteSelect(_fitness, calcSum(_fitness));
    else
    {
        double randPoint = randDouble()*_sum;
        double currentPoint=0.0;
        for (unsigned int i=0; i<_fitness.size(); i++)
        {
            currentPoint+=_fitness[i];
            if (currentPoint >= randPoint)
                return i;
        }

        return _fitness.size(); //Otherwise return an index we know will cause an exception because something has gone horribly wrong!
    }
}

//Returns an index using inverse roulette selection (disproportionately selecting those with lowest fitness).
unsigned int inverseRouletteSelect(const std::vector<double>& _fitness)
{
    //Find (max+1-d) sum and max value at the same time.
    double sumInversed=0.0;
    double maxVal = 0.0;
    for (double d:_fitness)
    {
        sumInversed+=d;
        if (d > maxVal) //Find maxVal...
            maxVal = d;
    }
    sumInversed = (((maxVal+1)*_fitness.size()) - sumInversed); //Now add on all the maxVal+1's.

    double randPoint = randDouble()*sumInversed;
    double currentPoint=0.0;
    for (unsigned int i=0; i<_fitness.size(); i++)
    {
        currentPoint+=(maxVal+1-_fitness[i]);
        if (currentPoint >= randPoint)
            return i;
    }

    return _fitness.size(); //Otherwise return an index we know will cause an exception because something has gone horribly wrong!
}

//Returns the sum of a vector.
double calcSum(const std::vector<double>& _v)
{
    double sum=0.0;
    for (double d:_v)
        sum+=d;
    return sum;
}

//Returns a random double between [0 1]
double randDouble()
{
    return ((double) rand() / (RAND_MAX));
}

//Returns a random index from 0 to _max (excluding _max), excluding any numbers in the _exclusionList.
unsigned int randomNumberWithExclusions(unsigned int _max, std::vector<unsigned int> _exclusionList, unsigned int _minInclusive)
{
	//Remove any duplicates.
	_exclusionList.erase(std::unique(_exclusionList.begin(), _exclusionList.end()), _exclusionList.end());

	//Sort the list.
    std::sort(_exclusionList.begin(), _exclusionList.end());

    //Check _max-_exclusionList.size() isn't 0.
	if (_max - _exclusionList.size() == 0) //If so then there is only one possible value to return, so return it.
	{
		for (unsigned int val=0; val<_max; ++val)
		{
			if (_exclusionList[val] != val) //Exclusion list was sorted earlier, so this works.
				return val;
		}
	}

	//Remove any exclusions which are smaller than _minInclusive.
	_exclusionList.erase(std::remove_if(_exclusionList.begin(), _exclusionList.end(), [_minInclusive](unsigned int i){return i<_minInclusive;}), _exclusionList.end());

    unsigned int adjustedMaximum = _max-_exclusionList.size();
	unsigned int randNum = (rand() % (adjustedMaximum-_minInclusive)) + _minInclusive;

    for (unsigned int i=0; i<_exclusionList.size(); i++)
    {
        if (randNum >= _exclusionList[i] && _exclusionList[i] >= _minInclusive)
            randNum++;
        else
            break;
    }

    if (randNum == _max)
        throw std::runtime_error("_exclusionList includes all possible values in utilities.h - randomNumberWithExclusions.");

    return randNum;
}


///Vector/matrix read and write functions.
//Returns a vector of vectors containing double values of a csv file.
std::vector<std::vector<double>> csvToMatrix(const std::string &_data, unsigned int _cols)
{
    //Somewhere to store values to.
    std::vector<std::vector<double>> output;
    for (unsigned int col=0; col<_cols; col++)
        output.push_back({});

    //Create file stream and open it.
    std::ifstream file;
    file.open(_data, std::ifstream::in);

    //Read in lines
    std::string line;
    while (std::getline(file, line))//(std::getline(file, line))
    {
        //Separate each line into doubles and add to each column.
        std::stringstream ss(line);
        std::string item;
        for (unsigned int i=0; i<_cols; i++)
        {
            std::getline(ss, item, ' ');
            double value = std::stof(item); //Convert to a double.
            output[i].push_back(value);
        }
    }

    //Close stream.
    file.close();
    return output;
}
