#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <NTL/ZZXFactoring.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/vec_ZZ_pE.h>
#include <NTL/vec_vec_ZZ_pE.h>
#include <iostream>

using namespace std;
using namespace NTL;

/*try to find target in vec, return index if finding, else -1*/
long find(const vec_ZZ_pE& vec, const ZZ_pE& target) {
    for (long i = 0; i < vec.length(); i++) {
        if (vec[i] == target) {
            return i;
        }
    }
    return -1; // 表示未找到
}

/*sample random vector, each element is non-zero and pairwise distinct*/
vec_ZZ_pE random_nonzero_vec_ZZ_pE(long n) {
    vec_ZZ_pE v1;
    v1.SetLength(n);

    for (long i = 0; i < n; i++) {
    ZZ_pE elem;
    do {
        random(elem);
    } while (IsZero(elem) || find(v1, elem) >= 0);
    v1[i] = elem;
}

    return v1;
}


/*Random version: phi of RMFE for GF
    Degree D - F^n_{p^s} -> F_{p^r}
    a: input vector
    size of a: D*n*s; D vectors, n dim, s coeffs
*/
Vec<ZZ_pEX> RMFE_GF_PHI(ZZ p, long n, long r, long s, long D, long a[])
{   
    if(D*(n-1)>r/s){
        cout<< "Note:D(n-1)<=r/s\n";
        Vec<ZZ_pEX> a;
        return a;
    }
    ZZ q1 = power(p,s);
    ZZ q2 = power(p,r);

    ZZ_p::init(p);
    ZZ_pX F;
    BuildIrred(F,s);

    ZZ_pE::init(F);

    Vec<ZZ_pEX> V;

    vec_ZZ_pE v1 = random_nonzero_vec_ZZ_pE(n);
    cout<< v1 << '\n';

    for(int i = 0; i < D; i++){
        vec_ZZ_pE v2;
        for (int j = 0; j < n; j++){
            ZZ_pX f;
            ZZ_pE g;
            for (int k = 0; k < s; k++){
                SetCoeff(f,k,a[i*n*s+j*s+k]);
            }
            conv(g,f);
            v2.append(g);
        }
        cout<< v2 << '\n';
        ZZ_pEX h = interpolate(v1,v2);
        V.append(h);
    }

    return V;
}

/*  phi of RMFE for GF
    Degree D - F^n_{p^s} -> F_{p^r}
    a: input vector
    size of a: D*n*s; D vectors, n dim, s coeffs
    b: interpolation points
    size of b: n*s;  n dim, s coeffs
*/
Vec<ZZ_pEX> RMFE_GF_PHI(ZZ p, long n, long r, long s, long D,long b[], long a[])
{   
    if(D*(n-1)>r/s){
        cout<< "Note:D(n-1)<=r/s\n";
        Vec<ZZ_pEX> a;
        return a;
    }
    ZZ q1 = power(p,s);
    ZZ q2 = power(p,r);

    ZZ_p::init(p);
    ZZ_pX F;
    BuildIrred(F,s);

    ZZ_pE::init(F);

    Vec<ZZ_pEX> V;

    vec_ZZ_pE v1;
    for(int i = 0;i < n; i++){
        ZZ_pX f;
        ZZ_pE g;
        for (int j = 0; j < s; j++){
            SetCoeff(f,j,b[i*s+j]);
        }
        conv(g,f);
        v1.append(g);
    }
    cout<<"--interpolation points--"<< v1 << '\n';

    for(int i = 0; i < D; i++){
        vec_ZZ_pE v2;
        for (int j = 0; j < n; j++){
            ZZ_pX f;
            ZZ_pE g;
            for (int k = 0; k < s; k++){
                SetCoeff(f,k,a[i*n*s+j*s+k]);
            }
            conv(g,f);
            v2.append(g);
        }
        ZZ_pEX h = interpolate(v1,v2);
        V.append(h);
    }

    return V;
}

/*  psi of RMFE for GF
    Degree D - F^n_{p^s} -> F_{p^r}
    c: vector for poly
    size of c: D
    b: interpolation points
    size of b: n*s;  n dim, s coeffs
*/
vec_ZZ_pE RMFE_GF_PSI(const Vec<ZZ_pEX>& c, long n, long r, long s, long D, long b[]){
    
    if(D*(n-1)>r/s){
        cout<< "Note:D(n-1)<=r/s\n";
        vec_ZZ_pE a;
        return a;
    }
    vec_ZZ_pE v1;
    for(int i = 0;i < n; i++){
        ZZ_pX f;
        ZZ_pE g;
        for (int j = 0; j < s; j++){
            SetCoeff(f,j,b[i*s+j]);
        }
        conv(g,f);
        v1.append(g);
    }

    ZZ_pEX m;
    set(m);
    for(int i = 0; i < D; i++){
        m *= c[i];
    }
    cout << "--Poly--" << m << '\n';

    vec_ZZ_pE res;
    for(int i = 0; i < n; i++){
        res.append(eval(m,v1[i]));
    }
    return res;
}


int main(){
    //  Degree D - F^n_{p^s} -> F_{p^r}


    long p = 5;
    long n = 3;
    long r = 5;
    long s = 1;
    long D = 2;

    long a[]={2,3,4,4,3,2}; // size: D * n * s; D vectors, n dim, s coeffs
    long b[]={1,3,2}; // size: n * s; n dim, s coeffs
    Vec<ZZ_pEX> c = RMFE_GF_PHI(ZZ(p),n,r,s,D,b,a);
    cout<< c<< '\n';
    cout << RMFE_GF_PSI(c,n,s,D,b);

//     long p = 2;
//     long n = 3;
//     long r = 7;
//     long s = 2;
//     long D = 3;

//     long a[]={0,0,0,1,1,0,0,1,0,0,1,1,0,1,1,0,0,1}; // size: D * n * s; D vectors, n dim, s coeffs
//     long b[]={0,1,1,1,1,0}; // size: n * s; n dim, s coeffs
//     Vec<ZZ_pEX> c = RMFE_GF_PHI(ZZ(p),n,r,s,D,b,a);
//     cout<< c<< '\n';
//     cout << RMFE_GF_PSI(c,n,r,s,D,b)<<'\n';
//     return 0;
// }