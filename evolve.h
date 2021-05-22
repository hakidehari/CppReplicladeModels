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
    map<char, int> nucPosMap = {{'A', 0}, {'T', 1}, {'C', 2}, {'G', 3}};
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
    vector<char> seqList;
    map<char, int> nucPosMap = {{'A', 0}, {'T', 1}, {'C', 2}, {'G', 3}};
};

class HKY85 {
public:
    HKY85(string);
    string evolve(string);
private:
    string sequence;
    void calcPrbMatrix();
    map<char, vector<double>> prbMatrix;
    int t;
    vector<char> seqList;
    map<char, int> nucPosMap = {{'A', 0}, {'T', 1}, {'C', 2}, {'G', 3}};
    map<char, int> frequencies;
    void calculateBaseFrequencies(string);
};

class F81 {
public:
    F81();
    string evolve(string);
private:
    void calcPrbMatrix(double, int);
    map<char, vector<double>> prbMatrix;
    int t;
    vector<char> seqList;
    map<char, int> nucPosMap = {{'A', 0}, {'T', 1}, {'C', 2}, {'G', 3}};
};