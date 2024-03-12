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

    ZZ p(ZZ(2));
    long pp = 2;
    long k = 2;
    long s = 3;
    long n = 4;
    long r = 4*4*3;
    long D = 4;

    ZZ q3 = power(p,k);
    ZZ q2 = power(p,s);
    long q1;
    conv(q1,q2-ZZ(1));

    vector<long> Interpolation;
    vector<long> Irred;
    RMFE_GR_Init(p,k,n,r,s,D,Interpolation,Irred);

    vector<long> result;
    vector<long> Input = {1,2,0,3,0,0,1,2,3,3,2,1,
                          1,2,0,0,2,3,1,1,2,2,3,0,
                          0,1,0,2,1,0,1,1,2,3,1,2,
                          1,0,0,0,1,0,1,0,1,2,0,3};
    //vector<long> Input = {1,2,0,1,2,3,1,0};
    RMFE_GR_PHI(ZZ(p),k,n,r,s,D,Interpolation,Irred,Input,result);

    vector<long> Output;
    RMFE_GR_PSI(ZZ(p),k,n,r,s,D,Interpolation,Irred,result,Output);

    return 0;
}