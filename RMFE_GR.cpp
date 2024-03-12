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
#include <NTL/tools.h>

using namespace std;
using namespace NTL;

#include <NTL/ZZ_pE.h>

bool Inv(const ZZ_pE& element) {
    ZZ_pE result;

    try {
        // 尝试计算逆元
        inv(result, element);

        // 如果没有引发错误，说明元素可逆
        return true;
    } catch (...) {
        // 如果引发了错误，说明元素不可逆
        return false;
    }
}

bool VectorInvertible(const vec_ZZ_pE& vector) {
    long n = vector.length();

    // 遍历所有元素对
    for (long i = 0; i < n; i++) {
        for (long j = i + 1; j < n; j++) {
            ZZ_pE difference = vector[i] - vector[j];
            // 使用 Inv 函数检查差是否可逆
            if (!Inv(difference)) {
                // 如果差不可逆，返回 false
                return false;
            }
        }
    }

    // 所有差都可逆，返回 true
    return true;
}

vec_ZZ_pE GenerateInvertibleDifferenceVector(long length) {
    vec_ZZ_pE result;
    result.SetLength(length);

    int MAX_ATTEMPTS = 10;

    for (long i = 0; i < length; i++) {
        bool validElement = false;

        // 限制重新生成的次数，防止无限循环
        for (int attempt = 0; attempt < MAX_ATTEMPTS; attempt++) {
            random(result[i]);

            // 验证当前元素与之前的元素的差是否可逆
            validElement = true;
            for (long j = 0; j < i; j++) {
                if (!Inv(result[i] - result[j])) {
                    validElement = false;
                    break;
                }
            }

            if (validElement) {
                break;  // 当前元素满足条件，跳出重新生成的循环
            }
        }

        if (!validElement) {
            // 重新生成的次数超过限制，从头开始生成整个向量
            i = -1;  // 因为下一轮循环会执行i++
        }
    }

    return result;
}



/*  phi of RMFE for GR
    Degree D - GR(p^k,s)^n -> GR(p^k,r)  r/s>=D(n-1)
    a: input vector
    size of a: D*n*s; D vectors, n dim, s coeffs
    b: interpolation points
    size of b: n*s;  n dim, s coeffs
*/
Vec<ZZ_pEX> RMFE_GR_PHI(ZZ p, long k, long n, long r, long s, long D,long b[], long a[])
{   
    if(D*(n-1)>r/s){
        cout<< "Note:D(n-1)<=r/s\n";
        Vec<ZZ_pEX> a;
        return a;
    }
    ZZ q1 = power(p,s);
    ZZ q2 = power(p,r);
    ZZ q3 = power(p,k);

    ZZ_p::init(p);
    ZZ_pX F;
    BuildIrred(F,s);
    cout << "--IrrPoly--: "<< F << endl;
    ZZ_p::init(q3);
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

    if(!VectorInvertible(v1)){
        cout<< "Interpolation points are not invertible!\n";
        Vec<ZZ_pEX> a;
        return a;
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


/*  Random version
    phi of RMFE for GR
    Degree D - GR(p^k,s)^n -> GR(p^k,r)  r/s>=D(n-1)
    a: input vector
    size of a: D*n*s; D vectors, n dim, s coeffs
    b: interpolation points
    size of b: n*s;  n dim, s coeffs
*/
Vec<ZZ_pEX> RMFE_GR_PHI(ZZ p, long k, long n, long r, long s, long D, long a[])
{   
    if(D*(n-1)>r/s){
        cout<< "Note:D(n-1)<=r/s\n";
        Vec<ZZ_pEX> a;
        return a;
    }
    ZZ q1 = power(p,s);
    ZZ q2 = power(p,r);
    ZZ q3 = power(p,k);

    ZZ_p::init(p);
    ZZ_pX F;
    BuildIrred(F,s);
    cout << "--IrrPoly--: "<< F << endl;
    ZZ_p::init(q3);
    ZZ_pE::init(F);

    Vec<ZZ_pEX> V;

    vec_ZZ_pE v1 = GenerateInvertibleDifferenceVector(n*s);

    if(!VectorInvertible(v1)){
        cout<< "Interpolation points are not invertible!\n";
        Vec<ZZ_pEX> a;
        return a;
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


/*  psi of RMFE for GR
    Degree D - GR(p^k,s)^n -> GR(p^k,r)  r/s>=D(n-1)
    c: vector for poly
    size of c: D
    b: interpolation points
    size of b: n*s;  n dim, s coeffs
*/
vec_ZZ_pE RMFE_GR_PSI(const Vec<ZZ_pEX>& c, long n, long r, long s, long D, long b[]){
    
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
    if(!VectorInvertible(v1)){
        cout<< "Interpolation points are not invertible!\n";
        vec_ZZ_pE a;
        return a;
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
    //Degree D - GR(p^k,s)^n -> GR(p^k,r)  r/s>=D(n-1)
    // long p = 2;
    // long k = 2;
    // long n = 2;
    // long r = 9;
    // long s = 1;
    // long D = 4;

    // long a[]={1,2,2,3,3,0,1,3}; // size: D * n * s; D vectors, n dim, s coeffs
    // long b[]={2,3}; // size: n * s; n dim, s coeffs
    // Vec<ZZ_pEX> c = RMFE_GR_PHI(ZZ(p),k,n,r,s,D,b,a);
    // cout<< "PHI result:  "<< c << endl;
    // vec_ZZ_pE d = RMFE_GR_PSI(c,n,r,s,D,b);
    // cout<< "PSI result:  "<< d << endl;


    // long p = 2;
    // long k = 2;
    // long n = 2;
    // long r = 12;
    // long s = 3;
    // long D = 2;

    // long a[]={1,2,2,3,3,0,1,2,3,3,2,1}; // size: D * n * s; D vectors, n dim, s coeffs
    // long b[]={1,1,0,0,2,3}; // size: n * s; n dim, s coeffs
    // Vec<ZZ_pEX> c = RMFE_GR_PHI(ZZ(p),k,n,r,s,D,b,a);
    // cout<< "PHI result:  "<< c << endl;
    // vec_ZZ_pE d = RMFE_GR_PSI(c,n,r,s,D,b);
    // cout<< "PSI result:  "<< d << endl;



    return 0;
}