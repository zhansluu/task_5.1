#include <iostream>
#include "functions.h"
#include <boost/math/quadrature/gauss.hpp>
#include <cmath>

using namespace std;

int main()
{
    setlocale(0, "Russian");
    cout << "������� 5.1. ����������� ���������� ���������� ��� ������ ������������ ������ ��������� �������������� ������� �������� (�� ����)." << endl << endl;

    double a, b, x[2], mu[2];
    size_t n;
    cout << "������� ������ ������ ��������������: � = ";
    cin >> a;
    cout << "������� ������� ������ ��������������: b = ";
    cin >> b;
    cout << "������������� ����� �����: N = 2" << endl;
    n = 2;
    /*cin >> n;
    while (n <= 0)
    {
        cout << endl << "������� ������������ �������� N! ��������� �������.";
        cin >> n;
    }*/

    using boost::math::quadrature::gauss;
    auto exact = [](double x) { return sin(x)*pow(x, -0.5); };
    cout << "������ �������� ���������: " << gauss<double, 7>::integrate(exact, a, b) << endl;

    cout << "������� ������� ��������� ���� ������������ �������: " << endl;
    for(size_t i = 0; i < n; i++)
    {
        cout << "x[" << i << "] = ";
        cin >> x[i];
        for (size_t j = 0; j < i; j++)
        {
            while (x[i] == x[j])
            {
                cout << "����� ���� ��� ����! ������� ������ �������� ����:";
                cin >> x[i];
            }
        }
    }
    //������� ������� ������� �������
    cout << endl << "������� ������� �������:" << endl;
    mu[0] = 2*sqrt(b)-2*sqrt(a);
    mu[1] = (2*pow(b, 1.5)-2*pow(a, 1.5))/3;
    mu[2] = (2*pow(b, 2.5)-2*pow(a, 2.5))/5;
    mu[3] = (2*pow(b, 3.5)-2*pow(a, 3.5))/7;

    for (size_t i = 0; i < n; i++)
        cout << "mu[" << i << "] = " << mu[i] << "; ";

    double coeff_1, coeff_2, delta, delta_1, delta_2;
    delta = x[1]-x[0];
    delta_1 = mu[0]*x[1]-mu[1];
    delta_2 = mu[1]-x[0]*mu[0];
    coeff_1 = delta_1/delta;
    coeff_2 = delta_2/delta;
    cout << endl << endl << "������������ ����������� ���: " << endl;
    cout << "A1 = " << coeff_1 << "; A2 = " << coeff_2 << endl;

    cout << "���������� ��� ������ ��� �������� ���������: " << setprecision(12) << coeff_1*f(x[0])+coeff_2*f(x[1]) << endl;
    cout << "���������� ����������� �����������: " << fabs(coeff_1*f(x[0])+coeff_2*f(x[1]) - gauss<double, 7>::integrate(exact, a, b)) << endl;
    cout << "--------------------------------------------------------------------";

    cout << endl << "������� ������� �������:" << endl;
    for (size_t i = 0; i < 2*n; i++)
        cout << "mu[" << i << "] = " << mu[i] << "; ";

    double a1, a2, D;
    a1 = (mu[0]*mu[3]-mu[2]*mu[1])/(mu[1]*mu[1]-mu[2]*mu[0]);
    a2 = (mu[2]*mu[2]-mu[3]*mu[1])/(mu[1]*mu[1]-mu[2]*mu[0]);
    D = a1*a1-4*a2;
    if (D == 0) {cout << "���� ���������!"; return 0;}
    if (D < 0) {cout << "���� �� �����������!"; return 0;}
    x[0] = (-a1-sqrt(D))/2;
    x[1] = (-a1+sqrt(D))/2;
    if(x[0] < a || x[1] > b) {cout << "���� �� ����������� (a,b)!"; return 0;}
    cout << endl << "x[1] = " << x[0] << "; x[2] = " << x[1] << endl;

    coeff_1 = (mu[1]-mu[0]*x[1])/(x[0]-x[1]);
    coeff_2 = (mu[1]-mu[0]*x[0])/(x[1]-x[0]);
    if (coeff_1 <= 0 || coeff_2 <= 0) {cout << "������������ �� ������������!"; return 0;}
    cout << endl << endl << "������������ ����������� ���: " << endl;
    cout << "A1 = " << coeff_1 << "; A2 = " << coeff_2 << endl;
    //cout << fabs(mu[3]-(coeff_1*pow(x[0],1)+coeff_2*pow(x[1],3))) << endl;

    cout << "���������� ��� ������ �� ���� �������� ���������: " << setprecision(12) << coeff_1*f(x[0])+coeff_2*f(x[1]) << endl;
    cout << "���������� ����������� �����������: " << fabs(coeff_1*f(x[0])+coeff_2*f(x[1]) - gauss<double, 7>::integrate(exact, a, b)) << endl;


    return 0;
}
