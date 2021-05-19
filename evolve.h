//
// Created by Haki Dehari on 5/16/21.
//
#include <string>
#include <map>
#include <vector>
#ifndef CPPEVOLUTIONARYMODELS_EVOLVE_H
#define CPPEVOLUTIONARYMODELS_EVOLVE_H

#endif //CPPEVOLUTIONARYMODELS_EVOLVE_H

using namespace std;

class Kimura80 {
public:
    Kimura80(double);
    string evolve(string);
    double getAlpha();
    double getBeta();
private:
    void calcPrbMatrix(double, double, int);
    map<char, vector<double>> prbMatrix;
    int t;
    double alpha;
    double beta;
    vector<char> seqList;
};

class JukesCantor69 {
public:
    JukesCantor69(double);
    string evolve(string);
    double getAlpha();
private:
    void calcPrbMatrix(double, int);
    map<char, vector<double>> prbMatrix;
    int t;
    double alpha;
    double beta;
    vector<char> seqList;
};