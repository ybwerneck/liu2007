
// SRI.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <string>
#include <thread>         // std::thread
#include "utils.h"// std::thread
#include "matrix.h"
#include <condition_variable>
#include "display.h"
#include <omp.h>
#include <math.h>       /* sin */
#include <chrono>
using namespace std::chrono_literals;
#include <chrono>
#include <cstdlib>
#include <ctime>    
#include <locale>
Matrix* S, * I;

//--------Parametros simulacao
int L = 5;
double dt = 0.01, dx = 0.05;
//possivel problema com long double
int tam = L / dx;
//--------- Parametros modelo
Matrix* D11, * D22;
long double A = 3.0;
long double d = 0.3; //taxa de morte natural
long double lambda = 0.65; // taxa de contagio
long double gamma = 0.8;
long double r = 0.5;
//---------
FILE* Sarq, * Iarq, * Integralsarq, * D1arq, * D2arq, * Rarq;
std::string resupasta = "resultado";
std::string subpasta;

std::string filename0 = "data/Int.txt";
std::string filename = "data/S.txt";
std::string filename2 = "data/I.txt";
std::string filename3 = "data/D1.txt";
std::string filename5 = "data/D2.txt";
std::string filename4 = "data/R.txt";

int n = 40000;
int random(int min, int max) {
    return rand() % (max + 1 - min) + min;

}
void initValues() {

    
    long double R0 = (lambda * A) / (d * (d + gamma)), H = (lambda * r) / (d * (d + gamma));
    long double I0 = d / (2.0 * lambda) * (R0 - 1 - H + sqrt(pow(R0 - 1.0 -H, 2.0) - 4.0 * H));
    long double S0 = A / (d + lambda * I0);
   
    S = new Matrix(tam, tam, S0);
    I = new Matrix(tam, tam, I0);

        
        for (int i = 0; i < tam; i++) {
        for (int j = 0; j < tam; j++)
        {
            long double r1 = random(100, 100) * 0.01;
            long double r2 = random (40, 100) * 0.01;
            I ->operator()(i, j, r2 * I0);

        }
    }

    double centrox = tam / 2, centroy = tam / 2;
    double maxdistancia = pow(pow(centrox - 0, 2) + pow(centroy - 0, 2), 0.5);

    D11 = new Matrix(tam, tam, 0.02);
    D22 = new Matrix(tam, tam, 0.005);
    for (int i = 0; i < tam; i++) {

        for (int j = 0; j < tam; j++)
        {
            double distancia = pow(pow(centrox - i, 2) + pow(centroy - j, 2), 0.5);
            double frac = (distancia / maxdistancia);
            //D11 ->operator()(i, j, pow(frac, 0.5) * 0.05);
            // D22 ->operator()(i, j, 1.0 - pow(frac, 0.5) * 0.9);

        }
    }





}
long double f(long double S, long double I) {


 
   return A - (d * S )- (lambda*S*I);
 }
long double g(long double S, long double I) {

   return (lambda * S * I) -((d+gamma)*I) - r;  
  }

long double harmonica(double n1, double n2) {
    return (2.0 * n1 * n2) / (n1 + n2);
}
long double difussionHet(Matrix* ua, int x, int y, Matrix* coef) {
    Matrix* S = ua;
    long double difx, dify;
    long double Ye = harmonica(coef->operator()(x - 1, y), coef->operator()(x, y)) / (dx * dx);
    long double Yd = harmonica(coef->operator()(x + 1, y), coef->operator()(x, y)) / (dx * dx);
    long double Yc = harmonica(coef->operator()(x, y - 1), coef->operator()(x, y)) / (dx * dx);
    long double Yb = harmonica(coef->operator()(x, y + 1), coef->operator()(x, y)) / (dx * dx);
    difx = (-S->operator()(x, y) + S->operator()(x - 1, y)) * Yd + (-S->operator()(x, y) + S->operator()(x + 1, y)) * Ye;
    dify = (-S->operator()(x, y) + S->operator()(x, y + 1)) * Yc + (-S->operator()(x, y) + S->operator()(x, y - 1)) * Yb;

    return difx + dify;

}
long double difussion(Matrix* ua, int x, int y, Matrix* coef) {

    Matrix* A = ua;
    long double difx, dify;
    long double esquerda = A->operator()((x == 0) ? tam - 1 : x - 1, y);
    long double direita = A->operator()((x == tam-1) ? 0 : x + 1, y);
    long double cima = A->operator()(x, (y == 0) ? tam - 1 : y - 1);
    long double baixo = A->operator()(x,(y == tam-1) ? 0 : y + 1);
    long double ponto = A->operator()(x, y);

    long double Y = dt*coef->operator()(x,y) / (dx * dx);

    difx = (-2.0 *ponto + esquerda + direita) * Y;
    dify = (-2.0 * ponto + cima + baixo) * Y;
    return difx + dify;

}
void step() {
    int i, j;
    for (i = 0; i < tam ; i++) {
        for (j = 0; j < tam; j++) {         
            long double r1 = f(S->operator()(i, j), I->operator()(i, j));
            long double r2 = g(S->operator()(i, j), I->operator()(i, j));
            long double deltaS =  dt*(r1) + difussion(S, i, j, D11);
            long double deltaI =  dt*(r2) + difussion(I, i, j, D22);
            S->operator()(i, j, deltaS + S->operator()(i, j));
            I->operator()(i, j, deltaI + I->operator()(i, j));
    
        }
    }


}
void printIntegrals(Matrix* S, Matrix* I, FILE* Iarq, double t) {

    fprintf(Iarq, "%f %f %f", S->sum(), I->sum(), t);
    fprintf(Iarq, "\n");

};
void printCoeficients(FILE* Darq, FILE* Rarq, double t) {


             
};
void printMatrixtoFile(Matrix* S, Matrix* I, FILE* Uarq, FILE* Varq) {

    for (int v = 0; v < tam; v++)
    {
        int dx = 1;
        for (int j = 0; j < tam; j++)
        {

            fprintf(Uarq, "%f ", S->operator()(v, j)); //att
            fprintf(Varq, "%f ", I->operator()(v, j)); //att


        }
        fprintf(Uarq, "\n");
        fprintf(Varq, "\n");

    };




}
void makeResultFiles(std::string sub) {

    std::time_t t = std::time(nullptr);
    char mbstr[100],a[100],b[100];
    
    std::strftime(mbstr, sizeof(mbstr), "%H.%M-%d%b%Y", std::localtime(&t));


    std::string data = mbstr;

    subpasta = resupasta + "/" + sub + "/" + data;
    filename0 = subpasta + "/" + filename0;
    filename = subpasta + "/" + filename;
    filename2 = subpasta + "/" + filename2;
    filename3 = subpasta + "/" + filename3;
    filename4 = subpasta + "/" + filename4;
    filename5 = subpasta + "/" + filename5;

    char param[200];
    sprintf(param, "mkdir %s", resupasta);
    system(param);
    sprintf(param, "mkdir %s\\%s", resupasta, sub);
    system(param);
    sprintf(param, "mkdir %s\\%s\\%s", resupasta, sub, data);
    system(param);
    sprintf(param, "mkdir %s\\%s\\%s\\data", resupasta, sub, data);
    system(param);
    sprintf(param, "mkdir %s\\%s\\%s\\result", resupasta, sub, data);
    system(param);


    D1arq = fopen(filename3.c_str(), "w");
    D2arq = fopen(filename5.c_str(), "w");

    Sarq = fopen(filename.c_str(), "w");
    Iarq = fopen(filename2.c_str(), "w");
    Integralsarq = fopen(filename0.c_str(), "w");
    Rarq = fopen(filename4.c_str(), "w");


}

std::atomic<bool> report = false;
void reporta(int* k, int parada) {
    while (report) {
        double o = parada == -1 ? n : parada;
        system("cls");
        printf("\nDT: %.10f DX:%.10f TAM= %d \n Quadro %d \ %d   completed %.2f%% \n", dt, dx, tam, n, *k, (*k / o) * 100.0);

        std::this_thread::sleep_for(300ms);
    }
}
int kt = 0;
int main()
{
    initValues();

    makeResultFiles("Basico");
    report = true;
    int frames = 40; // quantos quadros terão na animação resultante, afeta muito o desempenho do gnuplot 

    double tick = (1.0 / frames) * n;
    int parada = -1;
    std::thread printer = std::thread(reporta, &kt, parada);

    for (int i = 0; i < n; i++) {

        //passo no tempo
        step();
        kt = i;

        //rotinas de impressão no arquivo 
        if (i % 40 == 0)
            printIntegrals(S, I, Integralsarq, i);
        printCoeficients(D1arq, Rarq, i);

        if (parada == -1 && (i % (n / frames) == 0))
        {
            printMatrixtoFile(D11, D22, D1arq, D2arq);
            printMatrixtoFile(S, I, Sarq, Iarq);

            if (i != n - 1)
            {
                fprintf(Sarq, "\n\n");
                fprintf(Iarq, "\n\n");
                fprintf(D1arq, "\n\n");
                fprintf(D2arq, "\n\n");
            }
        }
        if (parada != -1 && i == parada)
        {

            printMatrixtoFile(S, I, Sarq, Iarq);


        }

    }

    report = false;
    fclose(D1arq);
    fclose(D2arq);
    fclose(Sarq);
    fclose(Iarq);
    fclose(Integralsarq);


    printer.join();
    if (parada == -1) {
        saveGif((char*)filename.c_str(), (char*)(subpasta + "/result/S.gif").c_str(), dx, dt, tick);
        saveGif((char*)filename2.c_str(), (char*)(subpasta + "/result/I.gif").c_str(), dx, dt, tick);
        saveGif((char*)filename3.c_str(), (char*)(subpasta + "/result/D1.gif").c_str(), dx, dt, tick);
        saveGif((char*)filename5.c_str(), (char*)(subpasta + "/result/D2.gif").c_str(), dx, dt, tick);

    }
    else
    {
        saveFoto((char*)filename.c_str(), (char*)(subpasta + "/result/S.png").c_str(), dx, dt, parada);
        saveFoto((char*)filename2.c_str(), (char*)(subpasta + "/result/I.png").c_str(), dx, dt, parada);

    }
    // saveDl((char*)filename3.c_str(), (char*)(subpasta + "/result/Dif.png").c_str(), "d1", "d2");

    saveDl((char*)filename0.c_str(), (char*)(subpasta + "/result/Int.png").c_str(), "S", "I");
    saveTl((char*)filename4.c_str(), (char*)(subpasta + "/result/R.png").c_str(), "R0", "Rd", "V");


}