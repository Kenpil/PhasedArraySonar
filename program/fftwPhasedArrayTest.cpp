#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <fstream>
//#include <string>
//#include <sstream>
#include <vector>
#include <fftw3.h>
#define N 8

using namespace std;

int main(void)
{
    FILE *ampArrayFile = fopen("resultList.txt", "r");
    vector<float> ampArray(N, 0.0);
    FILE *phaseTimeFile = fopen("phaseList.txt", "r");
    vector<float> phaseTimeRealList(1000, 0.0);
    vector<float> phaseTimeImagList(1000, 0.0);
    FILE *fftResultFile = fopen("fftResult.txt", "w");

    bool judgeEOF = 0;
    int phaseTimeListLength = 0;
    while (judgeEOF = fscanf(phaseTimeFile, "%f %f", &phaseTimeRealList[phaseTimeListLength], &phaseTimeImagList[phaseTimeListLength]) != EOF)
    {
        phaseTimeListLength++;
    }
    printf("phaseTimeListLength: %d\n", phaseTimeListLength);

    fftw_complex in[N], out[N], in2[N]; /* double [2] */
    fftw_plan p, q;

    int timeSliceNum = 0;
    vector<float> powerArray(N, 0.0);
    while (judgeEOF = fscanf(ampArrayFile, "%f,%f,%f,%f,%f,%f,%f,%f", &ampArray[0], &ampArray[1], &ampArray[2], &ampArray[3], &ampArray[4], &ampArray[5], &ampArray[6], &ampArray[7]) != EOF)
    {
        for (int i = 0; i < N; i++)
        {
            in[i][0] = ampArray[i] * phaseTimeRealList[timeSliceNum];
            in[i][1] = ampArray[i] * phaseTimeImagList[timeSliceNum];
        }
        p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(p);
        printf("time:%3d   ", timeSliceNum);
        for (int i = 0; i < N; i++)
        {
            powerArray[i] = pow(out[i][0], 2) + pow(out[i][1], 2);
            //printf("freq: %3d %+9.5f %+9.5f I\n", i, out[i][0], out[i][1]);
            printf("%10f ", powerArray[i]);
        }
        printf("\n");
        fftw_destroy_plan(p);

        timeSliceNum++;
    }
    return 0;
}