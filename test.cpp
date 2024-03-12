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
    long k = 16;
    long s = 9;
    long n = 4;
    long r = 4*9*4;
    long D = 4;


    vector<long> Interpolation;
    vector<long> Irred;    

    RMFE_GR_Init_cache(p,k,n,r,s,D,Interpolation,Irred);

    return 0;
}