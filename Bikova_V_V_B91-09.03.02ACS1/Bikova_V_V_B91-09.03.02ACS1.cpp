// Bikova_V_V_B91-09.03.02ACS1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <omp.h>
#include <ctime>
#include <vector>
#include <chrono>
using namespace std;

const int numElements = 1000000;
float Asgl[numElements];
double Adbl[numElements];


float parallelSumFloat(const vector<float>& arr, int numThreads)
{
    float sum = 0.0f;
#pragma omp parallel for reduction(+:sum) num_threads(numThreads)
    for (int i = 0; i < arr.size(); ++i)
    {
        sum += arr[i];
    }
    return sum;
}
double parallelSumDouble(const std::vector<double>& arr, int numThreads)
{
    double sum = 0.0;
#pragma omp parallel for reduction(+:sum) num_threads(numThreads)
    for (int i = 0; i < arr.size(); ++i)
    {
        sum += arr[i];
    }
    return sum;
}


int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "Количество процессоров: " << omp_get_num_procs() << "\n";
    for (int i = 0; i < 1000000; i++) {
        Asgl[i] = rand();
        Adbl[i] = rand();
    }

    float sum1 = 0;
    unsigned int start_time1 = clock();
    for (int i = 0; i < 1000000; i++) {
        sum1 += Asgl[i];
    }
    cout << "Сумма 1: " << sum1 << "\n";
    unsigned int end_time1 = clock();
    unsigned int search_time1 = end_time1 - start_time1;
    cout << "Время 1: " << search_time1 << "\n";

    double sum2 = 0;
    unsigned int start_time2 = clock();
    for (int i = 0; i < 1000000; i++) {
        sum2 += Adbl[i];
    }
    cout << "Сумма 2: " << sum2 << "\n";
    unsigned int end_time2 = clock();
    unsigned int search_time2 = end_time2 - start_time2;
    cout << "Время 2: " << search_time2 << "\n";


    std::cout << "Num threads\tFloat time\tDouble time" << std::endl;

    for (int numThreads = 2; numThreads <= 32; numThreads *= 2)
    {
        // Замер времени для массива одинарной точности
        auto startFloat = chrono::high_resolution_clock::now();
        float sumFloat = parallelSumFloat(vector<float>(Asgl, Asgl + numElements), numThreads);
        auto endFloat = chrono::high_resolution_clock::now();
        auto durationFloat = chrono::duration_cast<chrono::milliseconds>(endFloat - startFloat);

        // Замер времени для массива двойной точности
        auto startDouble = chrono::high_resolution_clock::now();
        double sumDouble = parallelSumDouble(vector<double>(Adbl, Adbl + numElements), numThreads);
        auto endDouble = chrono::high_resolution_clock::now();
        auto durationDouble = chrono::duration_cast<chrono::milliseconds>(endDouble - startDouble);

        std::cout << numThreads << "\t\t" << durationFloat.count() << " ms\t\t" << durationDouble.count() << " ms" << std::endl;
    }
// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
