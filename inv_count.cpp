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
#include <exception>
#include <NTL/mat_ZZ_pE.h>

using namespace std;
using namespace NTL;

#include <NTL/ZZ_pE.h>


// 函数定义
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

// vec_ZZ_pE GenerateInvertibleDifferenceVector(ZZ p, long n,long s) {
//     vec_ZZ_pE result;
//     ZZ q3 = power(p,s);
//     long a;
//     conv(a,q3-ZZ(1));
//     vector<long> factors = FindFactor(a);

//     ZZ_pE b;
//     bool flag = 0;
//     while(flag!=1){
//         random(b);
//         cout<<b<<power(b,a)<<endl;
//         flag = 1;
//             if(!IsOne(power(b,a)))
//             {   
//                 flag = 0;
//                 continue;
//             }
//             for(int i=0;i<factors.size()-1;i++){
//                 if(IsOne(power(b,factors[i])))
//                 {
//                     flag = 0;
//                     break;
//                 }
//                 }
//     }

//     for(int i=1;i<n+1;i++){
//         result.append(power(b,i));
//     }
//     cout<<result<<endl;
//     return result;
// }

int main(){
    ZZ_p::init(ZZ(2));
    ZZ_pX F;
    BuildIrred(F,3);
    cout << F << '\n';

    ZZ_p::init(ZZ(4));
    ZZ_pE::init(F);

    long countInvertible = 0;

    for (long a0 = 0; a0 < 4; a0++) {
        for (long a1 = 0; a1 < 4; a1++) {
            for (long a2 = 0; a2 < 4; a2++) {
                // 更多系数可以继续添加
                ZZ_pX aa;
                SetCoeff(aa, 0, a2);
                SetCoeff(aa, 1, a1);
                SetCoeff(aa, 2, a0);
                // 更多 SetCoeff 可以继续添加
                ZZ_pE a;
                conv(a,aa);
                if (Inv(a)) {
                    cout << "Invertible element: " << a <<" inv: "<<inv(a) << endl;
                    countInvertible++;
                }
            }
        }
    }

    cout << "Total Invertible Elements: " << countInvertible << endl;

    long a0 = 0, a1 = 1, a2 = 2;
    ZZ_pX aa;
    SetCoeff(aa, 0, a0);
    SetCoeff(aa, 1, a1);
    SetCoeff(aa, 2, a2);
    ZZ_pE a;
    conv(a,aa);

    a0 = 1;
    a1 = 0;
    a2 = 3;
    ZZ_pX bb;
    SetCoeff(bb, 0, a0);
    SetCoeff(bb, 1, a1);
    SetCoeff(bb, 2, a2);
    ZZ_pE b;
    conv(b,bb);
    cout<<a<<b<<a*b<<endl;
    ZZ_pE x;
    inv(x,b);
    cout<<x<<endl;
    return 0;
}