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
#include <fstream>

void interpolate3(ZZ_pEX& f, const vec_ZZ_pE& a, const vec_ZZ_pE& b, ZZ p, long l, long s)
{
   long m = a.length();
   if (b.length() != m) LogicError("interpolate: vector length mismatch");

   if (m == 0) {
      clear(f);
      return;
   }
   vec_ZZ_pE prod;
   prod = a;
   ZZ_pE t1, t2;

   long k, i;

   vec_ZZ_pE res;
   res.SetLength(m);

   for (k = 0; k < m; k++) {

      const ZZ_pE& aa = a[k];

      set(t1);
      for (i = k-1; i >= 0; i--) {
         mul(t1, t1, aa);
         add(t1, t1, prod[i]);
      }

      clear(t2);
      for (i = k-1; i >= 0; i--) {
         mul(t2, t2, aa);
         add(t2, t2, res[i]);
      }
      t1 = Inv(t1, p, l, s);


    ofstream file3("Inverse.txt", ios::app);
    if (!file3.is_open()) {
        std::cout << "无法打开文件。" << std::endl;
        return;
    }
    file3 << t1 << endl;;
    file3.close();


      sub(t2, b[k], t2);
      mul(t1, t1, t2);

      for (i = 0; i < k; i++) {
         mul(t2, prod[i], t1);
         add(res[i], res[i], t2);
      }

      res[k] = t1;
      if (k < m-1) {
         if (k == 0)
            prod[0]=-prod[0];
         else {
            t1=-a[k];
            add(prod[k], t1, prod[k-1]);
            for (i = k-1; i >= 1; i--) {
               mul(t2, prod[i], t1);
               add(prod[i], t2, prod[i-1]);
            }
            mul(prod[0], prod[0], t1);
         }
      }
   }

   while (m > 0 && IsZero(res[m-1])) m--;
   res.SetLength(m);
   f.rep = res;
}



int main(){
    ZZ p(ZZ(2));
    long pp = 2;
    long k = 16;
    long s = 9;
    long n = 4;
    long r = 4*4*3;
    long D = 4;

    ZZ q1 = power(p,s);
    ZZ q2 = power(p,r);
    ZZ q3 = power(p,k);
    long q4;
    conv(q4,q1-ZZ(1));

    ZZ_p::init(p);
    ZZ_pX F;
    FindPrimitivePoly(F,p,s);


    ofstream file1("Irred.txt", ios::out);
    if (!file1.is_open()) {
        std::cout << "无法打开文件。" << std::endl;
        return 1;
    }
    file1 << F << endl;;
    file1.close();


    ZZ_pX f;
    SetCoeff(f,q4,1);
    SetCoeff(f,0,-1);

    ZZ_pX g_;
    HenselLift(g_,f,F,p,k-1);
    ZZ_pE::init(g_);

    ofstream file2("HenselLift.txt", ios::out);
    if (!file2.is_open()) {
        std::cout << "无法打开文件。" << std::endl;
        return 1;
    }
    file2 << g_ << endl;;
    file2.close();

    ZZ_pX H;
    SetCoeff(H,1,1);
    ZZ_pE a;
    conv(a,H);

    vec_ZZ_pE v1;
    v1.append(a);
    if(IsOne(a) && n==2)
    {
        v1.append(ZZ_pE(0));
    }

    else{
        ZZ_pE a_tmp = a;
        for(int i=1;i<q4;i++){
            a_tmp = a_tmp * a;
            v1.append(a_tmp);
        }
    }
   
    vec_ZZ_pE result;
    result.SetLength(q4);
    random(result,q4);

    ZZ_pEX h;
    interpolate3(h,v1,result,p,k,s);


    return 0;
}