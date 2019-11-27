#ifdef BUILD_PYTHON_MODULE
#include "Omega.h"
#include <vector>
#include <iostream>

using namespace std;


int main() {
    vector<double> data = {0,2,4,6,8,10,9,8,7,6,5,4,3,2,1,0};
    vector<double> states = {0,1,2,3,4,5,6,7,8,9,10};

    double penalty = 0;

    Omega omega = Omega(states, penalty, data.size());
    omega.algo(data);
    vector<int> cp = omega.GetChangepoints();

    for (auto x : cp) {
        cout << x << ", ";
    }
    cout << endl;
}
#endif //BUILD_PYTHON_MODULE
