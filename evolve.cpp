//
// Created by Haki Dehari on 5/16/21.
//

#include "evolve.h"
#include <cmath>


Kimura80::Kimura80(double mu) {
    t = 0;
    alpha = mu;
    beta = mu / 3;
    calcPrbMatrix();
}

void Kimura80::calcPrbMatrix() {
    double transition = .25 + .25*(exp(-4*beta*t)) - .5*(exp(-2*(alpha + beta)*t));
    double transversion = .25 - .25*(exp(-4*beta*t));
    double sameNuc = 1 - transition - 2*transversion;
    vector<double> A = {sameNuc, transversion, transversion, transition};
    vector<double> T = {transversion, sameNuc, transition, transversion};
    vector<double> C = {transversion, transition, sameNuc, transversion};
    vector<double> G = {transition, transversion, transversion, sameNuc};
    prbMatrix["A"] = A;
    prbMatrix["T"] = T;
    prbMatrix["C"] = C;
    prbMatrix["G"] = G;
    t += 1;
}

