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
#include <vector>
#include <fstream>


using namespace std;
using namespace NTL;

/*  
    Find all factors of long n
    eg: n = 8,return [1,2,4,8]
*/
vector<long> FindFactor(long n) {
    std::vector<long> factors;

    for (long i = 1; i <= n; ++i) {
        if (n % i == 0) {
            factors.push_back(i);
        }
    }

    return factors;
}

Vec<ZZ> FindFactor(ZZ n) {
    Vec<ZZ> factors;

    ZZ i(ZZ(1));
    for (i; i <= n; ++i) {
        if (n % i == 0) {
            factors.append(i);
        }
    }

    return factors;
}


/*  
    Generate the exceptional set of GR(p^k,s)
    Suppose ZZ_pE has been inited and Irred has been fixed.
*/
vec_ZZ_pE GenerateExceptSet(ZZ p, long s) {
    vec_ZZ_pE result;
    ZZ q3 = power(p,s);
    long a;
    conv(a,q3-ZZ(1));
    vector<long> factors = FindFactor(a);

    ZZ_pE b;
    bool flag = 0;
    while(flag!=1){
        random(b);
        flag = 1;
            if(!IsOne(power(b,a)))
            {   
                flag = 0;
                continue;
            }
            for(long i=0;i<factors.size()-1;i++){
                if(IsOne(power(b,factors[i])))
                {
                    flag = 0;
                    break;
                }
                }
    }

    ZZ_pE m;
    clear(m);
    result.append(m);
    for(long i=0;i<a;i++){
        result.append(power(b,i));
    }
    //cout<<result<<endl;
    return result;
}

/*  
    Find all the invertible elemnts in GR(p^k,s)
    Suppose ZZ_pE has been inited and Irred has been fixed.
*/
void GetUnitGroupHelper(vec_ZZ_pE& result, const vec_ZZ_pE& ExceptSet, ZZ_pE element,ZZ p, long k, long s, long currentLevel) {
    if (currentLevel == k) {
        result.append(element);
    } else if(currentLevel == 0){
        ZZ q3 = power(p,s);
        long c;
        conv(c,q3);
        for (long i = 1; i < c; i++) {
            ZZ_pE a = element;
            a += ExceptSet[i];
            GetUnitGroupHelper(result,ExceptSet, a,p, k, s, currentLevel + 1);
        }
    } else {
        ZZ q3 = power(p,s);
        long c;
        conv(c,q3);
        for (long i = 0; i < c; i++) {
            ZZ a = power(p,currentLevel);
            ZZ_pE b;
            conv(b,a);
            ZZ_pE d = element;
            d += b*ExceptSet[i];
            GetUnitGroupHelper(result,ExceptSet, d,p, k, s, currentLevel + 1);
        }
    }
}

vec_ZZ_pE GetUnitGroup(ZZ p, long k, long s){
    vec_ZZ_pE result;
    ZZ_pE element;
    clear(element);
    vec_ZZ_pE ExceptSet = GenerateExceptSet(p,s);
    GetUnitGroupHelper(result,ExceptSet, element,p, k, s, 0);
    //cout<<result<<endl;
    return result;
}

/*  
    Get the inverse of element in GR(p^k,s)
    Element must be in the unit group, one must verify in advance.
*/
ZZ_pE Inv(const ZZ_pE& element, ZZ p, long k, long s){
    ZZ q1 = power(p,k*s)-power(p,s*(k-1));
    // ZZ qq1 = power(p,s)-ZZ(1);
    // ZZ qq2 = power(p,s*(k-1));

    // long q2;
    // conv(q2,qq1);
    // long q3;
    // conv(q3,qq2);

    // vector<long> factors1 = FindFactor(q2);
    // ZZ_pE a;
    // set(a);
    // long current = 0;
    // for(long i = 0;i<factors1.size();i++){
    //     ZZ_pE tmp = a;
    //     ZZ_pE b = power(element,factors1[i]-current);
    //     a *= b;
    //     if(IsOne(a)){
    //         ZZ_pE c = power(element,factors1[i]-current-1);
    //         return tmp * c;
    //     }
    //     current = i;
    // }

    // vector<long> factors2 = FindFactor(q3);
    // set(a);
    // current = 0;
    // for(long i = 0;i<factors2.size();i++){
    //     ZZ_pE tmp = a;
    //     ZZ_pE b = power(element,factors2[i]-current);
    //     a *= b;
    //     if(IsOne(a)){
    //         ZZ_pE c = power(element,factors2[i]-current-1);
    //         return tmp * c;
    //     }
    //     current = i;
    // }
    ZZ_pE a;
    set(a);
    ZZ current1(ZZ(0));
    ZZ i(ZZ(1));
    // if(qq1 > qq2)
    //     i = qq1;
    // else i = qq2;

    for(i;i<q1;i++){
        if(q1%i!=0)
            continue;
        ZZ_pE tmp = a;
        ZZ_pE b = power(element,i-current1);
        a *= b;
        if(IsOne(a)){
            ZZ_pE c = power(element,i-current1-1);
            return tmp * c;
        }
        current1 = i;
        
    }
    return ZZ_pE(0);
}

/*
    Determine whether element is an invertible element in GR(p^k,s)
*/
bool IsInv(const ZZ_pE& element, ZZ p, long k, long s){
    vec_ZZ_pE UnitGroup = GetUnitGroup(p,k,s);

    ZZ q1 = power(p,k*s)-power(p,s*(k-1));
    long q2;
    conv(q2,q1);

    for(long i = 0; i < q2;i++){
        ZZ_pE a = UnitGroup[i];
        if(a == element){
            return true;
        }
    }

    return false;

}
/*
    Interpolation 
*/
void interpolate2(ZZ_pEX& f, const vec_ZZ_pE& a, const vec_ZZ_pE& b, ZZ p, long l, long s)
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

/*
    Hensel Lift
    Input: f = g*h mod m, sg+th = 1 mod m
    Output: f = g_*h_ mod m^2, s_*g_+t_*h_ = 1 mod m^2
*/

void HenselLift(ZZ_pX& g_, ZZ_pX& h_, ZZ_pX& s_, ZZ_pX& t_, const ZZ_pX& f, 
                const ZZ_pX& g, const ZZ_pX& h, const ZZ_pX& s, const ZZ_pX& t){
    ZZ_pX e = f - g*h;

    ZZ_pX q,r;
    DivRem(q,r,s*e,h);

    g_ = g + t*e + q*g;
    h_ = h + r;

    ZZ_pX b = s*g_ + t*h_ -1;

    ZZ_pX c,d;
    DivRem(c,d,s*b,h_);

    s_ = s - d;
    t_ = t - t*b - c*g_;

    return;
}
/*
    Hensel Lift of g: g|f in Zp[x], g_|f in Zp^(n+1)[x]
*/
void HenselLift(ZZ_pX& g_, const ZZ_pX& f, const ZZ_pX& g, const ZZ p, long n){
    ZZ_p::init(p);

    ZZ_pX h = f / g;
    ZZ_pX d,s,t;
    XGCD(d,s,t,g,h);

    ZZ_pX g_tmp = g;
    ZZ_pX h_tmp = h;
    ZZ_pX f_tmp = f;
    for(int i = 0; i < n; i++){
        ZZ_pX h_, s_, t_;
        ZZ_p::init(power(p,i+2));
        SetCoeff(f_tmp,0,-1);
        HenselLift(g_,h_,s_,t_,f_tmp,g_tmp,h_tmp,s,t);
        g_tmp = g_;
        h_tmp = h_;
        s = s_;
        t = t_;
    }
    return;

}

/*
    Find the primitive polynomial of degree n in Zp[x]
*/

void FindPrimitivePoly(ZZ_pX& g, ZZ p, long n){
    ZZ q2 = power(p,n);
    long q1;
    conv(q1,q2-ZZ(1));

    vector<long> factors = FindFactor(q1);

    bool flag = 0;
    while(!flag){
        ZZ_pX F;
        BuildIrred(F,n);

        ZZ_pX f;
        SetCoeff(f,0,-1);
        flag = 1;
        for(int i=0;i<factors.size()-1;i++){
            if(factors[i]<=n)
                continue;
            
            SetCoeff(f,factors[i],1);
            ZZ_pX q,r;
            DivRem(q,r,f,F);
            if(r == 0){
                flag = 0;
                break;
            }
        }
        if(flag == 1){
            g = F;
        }

    }
    return;
}

/*
    ZZ_pX to ZZ_pEX
*/

void ZZ_pX2ZZ_pEX(ZZ_pEX& F, ZZ_pX& f){
    long n = deg(f);
    for (int i = 0;i<=n;i++){
        ZZ_p a = coeff(f,i);
        long b;
        conv(b,a);
        SetCoeff(F,i,b);
    }
    return;
}

/*
    Find a root of F in GR(p^k,s)
*/
ZZ_pE FindRootHelper(ZZ_pE b, const ZZ_pEX& F, long currentLevel, ZZ p, long k, long s){
    ZZ q1 = power(p,k);
    long q2;
    conv(q2,q1);

    if(currentLevel == 0){
        for(int i = 0; i<q2; i++){
            ZZ_pX f;
            SetCoeff(f,currentLevel,i);
            ZZ_pE a;
            conv(a,f);
            b = b + a;
            if(IsZero(eval(F,b))){
                return b;
            }
        }
    }
    else if(currentLevel < 0){
        if(IsZero(eval(F,b))){
                return b;
            } 
    }
    else {
        for(int i = 0; i<q2; i++){
            ZZ_pX f;
            SetCoeff(f,currentLevel,i);
            ZZ_pE a;
            conv(a,f);
            b = b + a;
            if(IsZero(eval(F,b))){
                return b;
            }
            cout<<b<<" ";
            ZZ_pE c = FindRootHelper(b,F,currentLevel-1,p,k,s);
            if(!IsZero(c)){
                return c;
            }
        }
    }
    return ZZ_pE(0);

}

void FindRoot2(ZZ_pE& a, const ZZ_pEX& F, ZZ p, long k, long s){
    ZZ q1 = power(p,k);
    long q2;
    conv(q2,q1);


    for(int j=s-1;j<s;j++){
        for(int i=1;i<q2;i++){
            ZZ_pX f;
            SetCoeff(f,j,i);
            ZZ_pE b;
            conv(b,f);

            ZZ_pE c = FindRootHelper(b,F,j-1,p,k,s);
            if(!IsZero(c)){
                a = c;
                return;
            }
        }
        
    }
    return;
}

void FindRoot3(ZZ_pE& a, const ZZ_pEX& F, ZZ p, long k, long s){
    ZZ_pE b;
    while(!IsZero(eval(F,b))){
        random(b);
        cout<<b<<" ";
    }
    a = b;

    return;
}


void fillResult(const Vec<ZZ_pEX>& V, vector<long>& result,long D, long n, long s) {
    result.clear(); // 清空 result

    for (long i = 0; i < D; i++) {
        const ZZ_pEX poly = V[i];
        vec_ZZ_pE coeffs;
        coeffs.SetLength(n);

        for (long j = 0; j < n; j++) {
            coeffs[j] = coeff(poly, j); // 获取多项式 poly 的第 j 个系数
        }


        // 将 ZZ_pE 系数转换为 long，并添加到 result 中
        for (long j = 0; j < n; j++) {
        const ZZ_pE element = coeffs[j];
        const ZZ_pX polyRep = rep(element);

        for (long k = 0; k < s; k++) {
            long c;
            conv(c,coeff(polyRep, k));
            result.push_back(c);
        }
    }
    }
}

void Long2VecZZ_pEX(const vector<long>& result, Vec<ZZ_pEX>& V, long D, long n, long s) {
    if (result.size() != D * n * s) {
        // 输入大小不符合预期，抛出错误或进行其他处理
        // 这里假设你想抛出一个异常
        cout<<"Input size error\n";
        return;
    }


    long index = 0; // 用于遍历输入向量的索引

    for (long i = 0; i < D; i++) {
        ZZ_pEX poly;
        for (long j = 0; j < n; j++) {
            ZZ_pX polyRep;
            for (long k = 0; k < s; k++) {
                // 从输入向量中取值并添加到多项式的系数中
                SetCoeff(polyRep, k, to_ZZ_p(result[index++]));
            }
            ZZ_pE a;
            conv(a,polyRep);
            // 将 ZZ_pX 系数添加到多项式中
            SetCoeff(poly, j, a);
        }
        // 将多项式添加到输出 Vec<ZZ_pEX> 中
        V.append(poly);
    }
}

void print(vector<long>& vec){
    cout<<"[";
    for (const auto& element : vec) {
        cout << element << " ";
    }
    cout<<"]"<< endl;
}

ZZ_pE long2ZZpE(const vector<long>& coefficients) {
    ZZ_pX poly;
    for (size_t i = 0; i < coefficients.size(); i++) {
        SetCoeff(poly, i, to_ZZ_p(coefficients[i]));
    }

    ZZ_pE result;
    conv(result, poly);

    return result;
}

ZZ_pX long2ZZpX(const vector<long>& coefficients) {
    ZZ_pX poly;
    for (size_t i = 0; i < coefficients.size(); i++) {
        SetCoeff(poly, i, to_ZZ_p(coefficients[i]));
    }

    return poly;
}


bool VectorInvertible(const vec_ZZ_pE& vector, ZZ p, long k, long s) {
    long n = vector.length();

    vec_ZZ_pE result;
    result.SetLength(n);
    random(result,n);
        try {
            ZZ_pEX h;
            interpolate2(h,vector,result,p,k,s);

        return true;
    } catch (...) {
        return false;
    }
    return true;
}

void fillIrred(const ZZ_pX& F, vector<long>& Irred) {
    Irred.clear(); // 清空 Irred

    for (long i = 0; i <= deg(F); i++) {
        long c;
        conv(c,(coeff(F, i)));
        Irred.push_back(c); // 取多项式 F 的第 i 个系数的值
    }
}

void fillInterpolation(const vec_ZZ_pE& v1, vector<long>& Interpolation, long s) {
    Interpolation.clear(); // 清空原有的数据

    for (long i = 0; i < v1.length(); i++) {
        const ZZ_pE element = v1[i];
        const ZZ_pX polyRep = rep(element);


        for (long j = 0; j < s; j++) {
            long c;
            conv(c,coeff(polyRep, j));
            Interpolation.push_back(c);
        }
    }
}

vector<long> PadVectorToLength(const vector<long>& input, size_t targetLength) {
    vector<long> paddedVector(targetLength, 0);  // 初始化长度为targetLength的全0向量

    // 将input的内容复制到paddedVector中，不超过targetLength
    size_t copyLength = min(input.size(), static_cast<size_t>(targetLength));
    copy(input.begin(), input.begin() + copyLength, paddedVector.begin());

    return paddedVector;
}

vector<long> TrimVector(const std::vector<long>& input) {
    // 找到最后一个非零元素的位置
    size_t lastNonZeroIndex = input.size();

    for (size_t i = input.size(); i > 0; --i) {
        if (input[i - 1] != 0) {
            lastNonZeroIndex = i - 1;
            break;
        }
    }

    // 截取向量
   vector<long> trimmedVector(input.begin(), input.begin() + lastNonZeroIndex + 1);

    return trimmedVector;
}

void RMFE_GR_Init(ZZ p, long k, long n, long r, long s, long D, vector<long>& Interpolation, vector<long>& Irred){
    if(D*(n-1)>r/s){
        cout<< "Note:D(n-1)<=r/s\n";
        return;
    }
    if(power(p,s)<n){
        cout<< "Lenstra Error\n";
        return;
    }
    ZZ q1 = power(p,s);
    ZZ q2 = power(p,r);
    ZZ q3 = power(p,k);
    long q4;
    conv(q4,q1-ZZ(1));

    ZZ_p::init(p);
    ZZ_pX F;
    FindPrimitivePoly(F,p,s);
    cout<<"--Irred--:"<<F;


    ZZ_pX f;
    SetCoeff(f,q4,1);
    SetCoeff(f,0,-1);

    ZZ_pX g_;
    HenselLift(g_,f,F,p,k-1);

    fillIrred(g_,Irred);
    ZZ_pE::init(g_);
    cout<<"  --Hensel Lift--:"<<g_<<endl;

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
        for(int i=1;i<n;i++){
            a_tmp = a_tmp * a;
            v1.append(a_tmp);
        }
    }
    cout<<"--interpolation points--"<< v1 << '\n';
    if(!VectorInvertible(v1,p,k,s)){
        cout<< "Interpolation points are not invertible!\n";
        return;
    }
    fillInterpolation(v1,Interpolation,s);

}

/* 
    phi of RMFE for GR
    Degree D - GR(p^k,s)^n -> GR(p^k,r)  r/s>=D(n-1)
    size of Input: D*n*s; D vectors, n dim, s coeffs
    size of Interpolation: n*s;  n dim, s coeffs
    size of result: D*n*s
    size of Irred: s
*/
void RMFE_GR_PHI(ZZ p, long k, long n, long r, long s, long D, const vector<long>& Interpolation, 
                 const vector<long>& Irred,const vector<long>& Input, vector<long>& result)
{   
    ZZ q1 = power(p,s);
    ZZ q2 = power(p,r);
    ZZ q3 = power(p,k);

    ZZ_p::init(q3);
    ZZ_pX F = long2ZZpX(Irred);
    ZZ_pE::init(F);

    Vec<ZZ_pEX> V;

    long q4;
    conv(q4,q1-ZZ(1));
    vec_ZZ_pE v1 ;
    v1.SetLength(n);
    for(int i = 0; i < n; i++){
        vector<long> coeff(Interpolation.begin()+i*s, Interpolation.begin()+i*s+s);
        v1[i] = long2ZZpE(coeff);
    }

    for(int i = 0; i < D; i++){
        vec_ZZ_pE v2;
        for (int j = 0; j < n; j++){
            ZZ_pX f;
            ZZ_pE g;
            for (int k = 0; k < s; k++){
                SetCoeff(f,k,Input[i*n*s+j*s+k]);
            }
            conv(g,f);
            v2.append(g);
        }
        ZZ_pEX h;
        cout<<v1<<endl;
        cout<<v2<<endl;
        interpolate2(h,v1,v2,p,k,s);
        V.append(h);
    }
    cout<<"--poly--"<<V<<endl;
    fillResult(V,result,D,n,s);
    return;
}

/*  psi of RMFE for GR
    Degree D - GR(p^k,s)^n -> GR(p^k,r)  r/s>=D(n-1)
    c: vector for poly
    size of c: D
    b: interpolation points
    size of b: n*s;  n dim, s coeffs
*/
void RMFE_GR_PSI(ZZ p, long k, long n, long r, long s, long D, const vector<long>& Interpolation, 
                 const vector<long>& Irred,const vector<long>& result, vector<long>& Output)

{

    ZZ q1 = power(p,s);
    ZZ q2 = power(p,r);
    ZZ q3 = power(p,k);

    ZZ_p::init(q3);
    ZZ_pX F = long2ZZpX(Irred);
    ZZ_pE::init(F);

    vec_ZZ_pE v1 ; //interpolation points
    v1.SetLength(n);
    for(int i = 0; i < n; i++){
        vector<long> coeff(Interpolation.begin()+i*s, Interpolation.begin()+i*s+s);
        v1[i] = long2ZZpE(coeff);
    }

    ZZ_pEX m;
    set(m);
    Vec<ZZ_pEX> V;
    Long2VecZZ_pEX(result,V,D,n,s);

    for(int i = 0; i < D; i++){
        m *= V[i];
    }
    cout << "--Poly--" << m << '\n';

    vec_ZZ_pE res;
    for(int i = 0; i < n; i++){
        res.append(eval(m,v1[i]));
    }

    cout<<"--Result--" <<res<<endl;

    fillInterpolation(res,Output,s);
}

void RMFE_GR_com_Init(ZZ p, long k, long n1,long n2, long r1, long r2, long s, long D, vector<long>& Interpolation1, 
                vector<long>& Interpolation2, vector<long>& Irred1, vector<long>& Irred2)
{
    RMFE_GR_Init(p,k,n1,r1,s,D,Interpolation1,Irred1);
    RMFE_GR_Init(p,k,n2,r2,r1,D,Interpolation2,Irred2);
}

/* 
    phi of RMFE_com for GR
    Degree D - GR(p^k,s)^n1 -> GR(p^k,r1)  r1/s>=D(n1-1)
               GR(p^k,r1)^n2 -> GR(p^k,r2)  r2/s>=D(n2-1)

    size of Input: D*n1*n2*s; D vectors, n1*n2 dim, s coeffs
    size of Interpolation: n*s;  n dim, s coeffs
    size of result: D*n*s
    size of Irred: s
*/
void RMFE_GR_com_PHI(ZZ p, long k, long n1,long n2, long r1, long r2, long s, long D,const vector<long>& Interpolation1, 
                const vector<long>& Interpolation2, const vector<long>& Irred1, const vector<long>& Irred2,
                const vector<long>& Input, vector<long>& result)
{   
    vector<long> Input2; //D*n2*n1*s
    for(int i = 0; i < D*n2; i++){
        vector<long> Input1(Input.begin()+i*n1*s,Input.begin()+i*n1*s+n1*s);
        vector<long> result1;
        RMFE_GR_PHI(ZZ(p),k,n1,r1,s,1,Interpolation1,Irred1,Input1,result1);

        vector<long> result1_pad;
        result1_pad = PadVectorToLength(result1,r1);
        Input2.insert(Input2.end(),result1_pad.begin(),result1_pad.end());
    }
    RMFE_GR_PHI(ZZ(p),k,n2,r2,r1,D,Interpolation2,Irred2,Input2,result);
}

// result: D*n2*s
void RMFE_GR_com_PSI(ZZ p, long k, long n1,long n2, long r1, long r2, long s, long D,const vector<long>& Interpolation1, 
                const vector<long>& Interpolation2, const vector<long>& Irred1, const vector<long>& Irred2,
                const vector<long>& result, vector<long>& Output)
{   
    //RMFE_GR_PSI(ZZ(p),k,n2,r2,r1,D,Interpolation2,Irred2,result,output1);

    ZZ q3 = power(p,k);

    ZZ_p::init(q3);
    ZZ_pX F = long2ZZpX(Irred2);
    ZZ_pE::init(F);

    vec_ZZ_pE v1 ; //interpolation points
    v1.SetLength(n2);
    for(int i = 0; i < n2; i++){
        vector<long> coeff(Interpolation2.begin()+i*r1, Interpolation2.begin()+i*r1+r1);
        v1[i] = long2ZZpE(coeff);
    }

    ZZ_pEX m;
    set(m);
    Vec<ZZ_pEX> V;
    Long2VecZZ_pEX(result,V,D,n2,r1);

    for(int i = 0; i < D; i++){
        m *= V[i];
    }
    cout << "--Poly--" << m << '\n';

    vec_ZZ_pE res;
    for(int i = 0; i < n2; i++){
        res.append(eval(m,v1[i]));
    }

    cout<<res<<endl;

    vec_ZZ_pE v2 ; //interpolation points
    v2.SetLength(n1);
    for(int i = 0; i < n1; i++){
        vector<long> coeff(Interpolation1.begin()+i*s, Interpolation1.begin()+i*s+s);
        v2[i] = long2ZZpE(coeff);
    }

    vec_ZZ_pE v3 ;
    for(int i=0;i<n2;i++){
        ZZ_pE a = res[i];
        ZZ_pX a_rep = rep(a);  // 获取 ZZ_pE 的内部表示 ZZ_pX

        ZZ_pEX result2;
        for(long j = 0; j <= deg(a_rep); j++){
            SetCoeff(result2,j,a_rep[j]);
        }

        for(int j=0;j<n1;j++){
            ZZ_pE c = eval(result2,v2[j]);
            v3.append(c);
        }

    }
    fillInterpolation(v3,Output,s);
}


void interpolate_for_cache(ZZ_pEX& f, const vec_ZZ_pE& a, const vec_ZZ_pE& b, ZZ p, long l, long s)
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

    ifstream file2("Inverse.txt",ios::in);

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

      file2>>t1;
      cout<<"t1:"<<t1<<endl;
      streampos currentPos = file2.tellg();

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
      file2.seekg(currentPos);
   }
    file2.close();
   while (m > 0 && IsZero(res[m-1])) m--;
   res.SetLength(m);
   f.rep = res;
}


void RMFE_GR_Init_cache(ZZ p, long k, long n, long r, long s, long D, vector<long>& Interpolation, vector<long>& Irred){
    if(D*(n-1)>r/s){
        cout<< "Note:D(n-1)<=r/s\n";
        return;
    }
    if(power(p,s)<n){
        cout<< "Lenstra Error\n";
        return;
    }
    ZZ q1 = power(p,s);
    ZZ q2 = power(p,r);
    ZZ q3 = power(p,k);
    long q4;
    conv(q4,q1-ZZ(1));

    ZZ_p::init(q3);
    ZZ_pX g_;

    ifstream file1("HenselLift.txt");
    file1 >> g_;
    file1.close();

    fillIrred(g_,Irred);
    ZZ_pE::init(g_);
    cout<<"  --Hensel Lift--:"<<g_<<endl;

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
        for(int i=1;i<n;i++){
            a_tmp = a_tmp * a;
            v1.append(a_tmp);
        }
    }
    fillInterpolation(v1,Interpolation,s);
}

void RMFE_GR_PHI_cache(ZZ p, long k, long n, long r, long s, long D, const vector<long>& Interpolation, 
                 const vector<long>& Irred,const vector<long>& Input, vector<long>& result)
{   
    ZZ q1 = power(p,s);
    ZZ q2 = power(p,r);
    ZZ q3 = power(p,k);

    ZZ_p::init(q3);
    ZZ_pX F = long2ZZpX(Irred);
    ZZ_pE::init(F);

    Vec<ZZ_pEX> V;

    long q4;
    conv(q4,q1-ZZ(1));
    vec_ZZ_pE v1 ;
    v1.SetLength(n);
    for(int i = 0; i < n; i++){
        vector<long> coeff(Interpolation.begin()+i*s, Interpolation.begin()+i*s+s);
        v1[i] = long2ZZpE(coeff);
    }

    for(int i = 0; i < D; i++){
        vec_ZZ_pE v2;
        for (int j = 0; j < n; j++){
            ZZ_pX f;
            ZZ_pE g;
            for (int k = 0; k < s; k++){
                SetCoeff(f,k,Input[i*n*s+j*s+k]);
            }
            conv(g,f);
            v2.append(g);
        }
        ZZ_pEX h;
        interpolate_for_cache(h,v1,v2,p,k,s);
        V.append(h);
    }
    cout<<"--poly--"<<V<<endl;
    fillResult(V,result,D,n,s);
    return;
}

void RMFE_GR_com_Init_cache(ZZ p, long k, long n1,long n2, long r1, long r2, long s, long D, vector<long>& Interpolation1, 
                vector<long>& Interpolation2, vector<long>& Irred1, vector<long>& Irred2)
{
    RMFE_GR_Init(p,k,n1,r1,s,D,Interpolation1,Irred1);
    RMFE_GR_Init_cache(p,k,n2,r2,r1,D,Interpolation2,Irred2);
}

void RMFE_GR_com_PHI_cache(ZZ p, long k, long n1,long n2, long r1, long r2, long s, long D,const vector<long>& Interpolation1, 
                const vector<long>& Interpolation2, const vector<long>& Irred1, const vector<long>& Irred2,
                const vector<long>& Input, vector<long>& result)
{   
    vector<long> Input2; //D*n2*n1*s
    for(int i = 0; i < D*n2; i++){
        vector<long> Input1(Input.begin()+i*n1*s,Input.begin()+i*n1*s+n1*s);
        vector<long> result1;
        RMFE_GR_PHI(ZZ(p),k,n1,r1,s,1,Interpolation1,Irred1,Input1,result1);

        vector<long> result1_pad;
        result1_pad = PadVectorToLength(result1,r1);
        Input2.insert(Input2.end(),result1_pad.begin(),result1_pad.end());
    }
    RMFE_GR_PHI_cache(ZZ(p),k,n2,r2,r1,D,Interpolation2,Irred2,Input2,result);
}
