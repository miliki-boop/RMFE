#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <NTL/ZZXFactoring.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/vec_ZZ_pE.h>
#include <NTL/vec_vec_ZZ_pE.h>
#include <iostream>
#include <exception>
#include <NTL/mat_ZZ_pE.h>
#include <NTL/ZZ_pE.h>
#include <NTL/pair_ZZ_pEX_long.h>
#include <vector>
using namespace std;
using namespace NTL;

#include "util.h"
int main(){

    long p = 2;
    long k = 16;
    long n1 = 2;
    long r1 = 9;
    long s = 1;
    long D = 2;
    long n2 = 10;
    long r2 = 2*10*9;


    vector<long> Interpolation1,Interpolation2;
    vector<long> Irred1,Irred2;
    RMFE_GR_com_Init_cache(ZZ(p),k,n1,n2,r1,r2,s,D,Interpolation1,Interpolation2,Irred1,Irred2);

    vector<long> result;
    vector<long> Input = {0,3,1,8,1,5,9,1,15,1,19,26,3,45,5,6,7,8,9,60000,0,1,1,54,1,1,47,1,1,41,1,2,3,4,5,6,7,45,96,10};
    RMFE_GR_com_PHI_cache(ZZ(p),k,n1,n2,r1,r2,s,D,Interpolation1,Interpolation2,Irred1,Irred2,Input,result);

    vector<long> Output;
    RMFE_GR_com_PSI(ZZ(p),k,n1,n2,r1,r2,s,D,Interpolation1,Interpolation2,Irred1,Irred2,result,Output);

    print(Output);

    return 0;
}