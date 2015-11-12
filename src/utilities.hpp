#pragma once
#include <vector>
#include <string>

#include <fstream>
#include <sstream>

//Returns an index using inverse roulette selection (proportionately selecting those with lowest fitness).
unsigned int inverseRouletteSelect(const std::vector<double>& _fitness);

//Returns the sum of a vector.
double calcSum(const std::vector<double>& _v);

//Returns a random double between [0 1]
double randDouble();

//Returns a random index from 0 to _max (excluding _max), excluding any numbers in the _exclusionList.
unsigned int randomNumberWithExclusions(unsigned int _max, std::vector<unsigned int> _exclusionList, unsigned int _minInclusive=0);

///Vector/matrix read and write functions.
std::vector<std::vector<double>> csvToMatrix(const std::string &_data, unsigned int _cols);

//Template function.
template <typename T>
void vectorToFile(std::vector<T> _vector, std::string _filename)
{
    std::ofstream file;
    file.open(_filename, std::ofstream::out | std::ofstream::trunc);

    for (unsigned int i=0; i<_vector.size(); i++)
        file << _vector.at(i) << "\n";

    file.flush();
    file.close();
}

template <typename T>
void matrixToFile(const std::vector<std::vector<T>>& _m, std::string _filename, std::string _delim)
{
    std::ofstream file;
    file.open(_filename, std::ofstream::out | std::ofstream::trunc);

    //Assumes all subarrays are the same size
    const unsigned int cols = _m.size();
    const unsigned int rows = _m.empty() ? 0 : _m[0].size();

    for (unsigned int row=0; row<rows; row++)
    {
        for (unsigned int col=0; col<cols; col++)
        {
            if (col != 0)
                file << _delim;
            if (_m[col].size() <= row) //Substitute 0
                file << 0;
            else
                file << _m[col][row];
        }
        file << "\n";
    }

    file.flush();
    file.close();
}

//Writes a single variable to file, e.g. a string or a double or anything with overloaded << operator.
template <typename T>
void valToFile(const T val, std::string filename)
{
    std::ofstream file;
    file.open(filename, std::ofstream::out | std::ofstream::trunc);
    file << val;
    file.flush();
    file.close();
}

//Returns an index using roulette selection (proportionately selecting those with highest fitness).
template <class T>
unsigned int rouletteSelect(const std::vector<T>& _fitness, double _sum)
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
