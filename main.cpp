#include <iostream>
#include "evolve.h"

using namespace std;

int main() {

    string name = "name";
    cout << name << endl;

    Kimura80 ki(.55);
    cout << ki.getAlpha() << endl;
    return 0;
}
