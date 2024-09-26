#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include "include/nlohmann/json.hpp"  // Include nlohmann/json library

using json = nlohmann::json;
using namespace std;

// Function to convert a string number from a given base to decimal
long long convertToDecimal(const string& numStr, int base) {
    long long result = 0;
    for (char ch : numStr) {
        int digit;
        if (ch >= '0' && ch <= '9') digit = ch - '0';
        else if (ch >= 'A' && ch <= 'Z') digit = ch - 'A' + 10;
        else if (ch >= 'a' && ch <= 'z') digit = ch - 'a' + 10;
        else continue; // Skip invalid characters
        result = result * base + digit;
    }
    return result;
}

// Function to perform Gaussian elimination
vector<double> gaussianElimination(vector<vector<double>>& A, vector<double>& b) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        // Pivot
        int maxRow = i;
        for (int k = i+1; k < n; k++) {
            if (abs(A[k][i]) > abs(A[maxRow][i])) maxRow = k;
        }
        swap(A[i], A[maxRow]);
        swap(b[i], b[maxRow]);

        // Eliminate
        for (int k = i+1; k < n; k++) {
            double c = A[k][i] / A[i][i];
            for (int j = i; j < n; j++) {
                A[k][j] -= c * A[i][j];
            }
            b[k] -= c * b[i];
        }
    }

    // Back substitution
    vector<double> x(n);
    for (int i = n-1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i+1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
    return x;
}

int main() {
    // Read JSON file
    ifstream inputFile("input.json");
    json jsonData;
    inputFile >> jsonData;

    // Extract n and k
    int n = jsonData["keys"]["n"];
    int k = jsonData["keys"]["k"];

    // Extract and decode (x, y) pairs
    vector<pair<int, long long>> dataPairs;
    for (auto& element : jsonData.items()) {
        if (element.key() == "keys") continue;
        int x = stoi(element.key());
        int base = stoi(element.value()["base"].get<string>());
        string valueStr = element.value()["value"].get<string>();
        long long y = convertToDecimal(valueStr, base);
        dataPairs.push_back({x, y});
    }

    // Check if we have enough data points
    if (dataPairs.size() < k) {
        cout << "Not enough data points." << endl;
        return 1;
    }

    // Select first k data points (or any k points)
    vector<pair<int, long long>> selectedData(dataPairs.begin(), dataPairs.begin() + k);

    // Set up the system of equations
    int degree = k - 1;
    vector<vector<double>> A(k, vector<double>(k));
    vector<double> b(k);

    for (int i = 0; i < k; i++) {
        int x = selectedData[i].first;
        long long y = selectedData[i].second;
        b[i] = static_cast<double>(y);
        for (int j = 0; j < k; j++) {
            A[i][j] = pow(x, degree - j);
        }
    }

    // Solve the system
    vector<double> coeffs = gaussianElimination(A, b);

    // The constant term 'c' is the last coefficient
    double c = coeffs.back();

    // Since the values are integers, we can round 'c' to the nearest integer
    long long constantTerm = static_cast<long long>(round(c));

    // Output the constant term 'c'
    cout << constantTerm << endl;

    return 0;
}
