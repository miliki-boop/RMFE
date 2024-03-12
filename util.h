#ifndef UTIL_H
#define UTIL_H
using namespace std;
using namespace NTL;

vector<long> FindFactor(long n);
vec_ZZ_pE GenerateExceptSet(ZZ p, long s);
vec_ZZ_pE GetUnitGroup(ZZ p, long k, long s);
ZZ_pE Inv(const ZZ_pE& element, ZZ p, long k, long s);
bool IsInv(const ZZ_pE& element, ZZ p, long k, long s);
void interpolate2(ZZ_pEX& f, const vec_ZZ_pE& a, const vec_ZZ_pE& b, ZZ p, long k, long s);
void HenselLift(ZZ_pX& g_, const ZZ_pX& f, const ZZ_pX& g, const ZZ p, long n);
void FindPrimitivePoly(ZZ_pX& g, ZZ p, long n);
void ZZ_pX2ZZ_pEX(ZZ_pEX& F, ZZ_pX& f);
void FindRoot2(ZZ_pE& a, const ZZ_pEX& F, ZZ p, long k, long s);
void print(vector<long>& vec);
void RMFE_GR_Init(ZZ p, long k, long n, long r, long s, long D, vector<long>& Interpolation, vector<long>& Irred);
void RMFE_GR_PHI(ZZ p, long k, long n, long r, long s, long D, const vector<long>& Interpolation, 
                 const vector<long>& Irred,const vector<long>& Input, vector<long>& result);
void RMFE_GR_PSI(ZZ p, long k, long n, long r, long s, long D, const vector<long>& Interpolation, 
                 const vector<long>& Irred,const vector<long>& result, vector<long>& Output);     
void RMFE_GR_com_Init(ZZ p, long k, long n1,long n2, long r1, long r2, long s, long D, vector<long>& Interpolation1, 
                vector<long>& Interpolation2, vector<long>& Irred1, vector<long>& Irred2);
void RMFE_GR_com_PHI(ZZ p, long k, long n1,long n2, long r1, long r2, long s, long D,const vector<long>& Interpolation1, 
                const vector<long>& Interpolation2, const vector<long>& Irred1, const vector<long>& Irred2,
                const vector<long>& Input, vector<long>& result);
void RMFE_GR_com_PSI(ZZ p, long k, long n1,long n2, long r1, long r2, long s, long D,const vector<long>& Interpolation1, 
                const vector<long>& Interpolation2, const vector<long>& Irred1, const vector<long>& Irred2,
                const vector<long>& result, vector<long>& Output);
void RMFE_GR_Init_cache(ZZ p, long k, long n, long r, long s, long D, vector<long>& Interpolation, vector<long>& Irred);
void RMFE_GR_PHI_cache(ZZ p, long k, long n, long r, long s, long D, const vector<long>& Interpolation, 
                 const vector<long>& Irred,const vector<long>& Input, vector<long>& result);
void RMFE_GR_com_Init_cache(ZZ p, long k, long n1,long n2, long r1, long r2, long s, long D, vector<long>& Interpolation1, 
                vector<long>& Interpolation2, vector<long>& Irred1, vector<long>& Irred2);

void RMFE_GR_com_PHI_cache(ZZ p, long k, long n1,long n2, long r1, long r2, long s, long D,const vector<long>& Interpolation1, 
                const vector<long>& Interpolation2, const vector<long>& Irred1, const vector<long>& Irred2,
                const vector<long>& Input, vector<long>& result);

#endif