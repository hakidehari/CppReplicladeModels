//
// Created by Haki Dehari on 5/16/21.
//

#include "evolve.h"
#include <cmath>

JukesCantor69::JukesCantor69(double mu) {
    alpha = mu;
    t = 0;
    calcPrbMatrix(alpha, t);
}

double JukesCantor69::getAlpha() {
    return alpha;
}

void JukesCantor69::calcPrbMatrix(double alpha, int t) {
    double sameNuc = .25 + .75*(exp(-4*alpha*t));
    double diffNuc = .25 - .25*(exp(-4*alpha*t));

    prbMatrix['A'] = {sameNuc, diffNuc, diffNuc, diffNuc};
    prbMatrix['T'] = {diffNuc, sameNuc, diffNuc, sameNuc};
    prbMatrix['C'] = {diffNuc, diffNuc, sameNuc, diffNuc};
    prbMatrix['G'] = {diffNuc, diffNuc, diffNuc, sameNuc};
}

string JukesCantor69::evolve(string seq) {
    string resultantString = "";
    for (int i = 0; i < seq.size(); i++) {
        char cur = seq[i];
        if (cur == '-') {
            resultantString += cur;
            continue;
        }
        double firstRoll = rand() / double(RAND_MAX);
        vector<double> probability = prbMatrix.find(cur)->second;
        if (firstRoll <= probability[nucPosMap.find(cur)->second]) {
            resultantString += cur;
        } else {
            for (int j = 0; j < prbMatrix.find(cur)->second.size(); j++) {
                if (prbMatrix.find(cur)->second[j] != 0) {
                    double roll = rand() / double(RAND_MAX);
                    if (roll <= prbMatrix.find(cur)->second[j] && prbMatrix.find(cur)->second[j] != prbMatrix.find(cur)->second[nucPosMap.find(cur)->second]) {
                        cur = seqList[j];
                        break;
                    }
                }
            }
            resultantString += cur;
        }
    }
    t += 1;
    calcPrbMatrix(alpha, t);
    return resultantString;
}

Kimura80::Kimura80(double mu) {
    t = 0;
    alpha = mu;
    beta = mu / 3;
    seqList = {'A', 'T', 'C', 'G'};
    calcPrbMatrix(alpha, beta, t);
}

double Kimura80::getAlpha() {
    return alpha;
}

double Kimura80::getBeta() {
    return beta;
}

void Kimura80::calcPrbMatrix(double alpha, double beta, int t) {
    double transition = .25 + .25*(exp(-4*beta*t)) - .5*(exp(-2*(alpha + beta)*t));
    double transversion = .25 - .25*(exp(-4*beta*t));
    double sameNuc = 1 - transition - 2*transversion;
    vector<double> A = {sameNuc, transversion, transversion, transition};
    vector<double> T = {transversion, sameNuc, transition, transversion};
    vector<double> C = {transversion, transition, sameNuc, transversion};
    vector<double> G = {transition, transversion, transversion, sameNuc};
    prbMatrix['A'] = A;
    prbMatrix['T'] = T;
    prbMatrix['C'] = C;
    prbMatrix['G'] = G;
    t += 1;
}

string Kimura80::evolve(string seq) {
    string resultantString = "";
    for (int i = 0; i < seq.size(); i++) {
        char cur = seq[i];
        if (cur == '-') {
            resultantString += cur;
            continue;
        }
        double firstRoll = rand() / double(RAND_MAX);
        vector<double> probability = prbMatrix.find(cur)->second;
        if (firstRoll <= probability[nucPosMap.find(cur)->second]) {
            resultantString += cur;
        } else {
            for (int j = 0; j < prbMatrix.find(cur)->second.size(); j++) {
                if (prbMatrix.find(cur)->second[j] != 0) {
                    double roll = rand() / double(RAND_MAX);
                    if (roll <= prbMatrix.find(cur)->second[j] && prbMatrix.find(cur)->second[j] != prbMatrix.find(cur)->second[nucPosMap.find(cur)->second]) {
                        cur = seqList[j];
                        break;
                    }
                }
            }
            resultantString += cur;
        }
    }
    t += 1;
    calcPrbMatrix(alpha, beta, t);
    return resultantString;
}

