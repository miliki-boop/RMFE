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
#include <vector>

using namespace std;
using namespace NTL;

#include <NTL/ZZ_pE.h>

void print(vector<long>& vec){
    cout<<"[";
    for (const auto& element : vec) {
        cout << element << " ";
    }
    cout<<"]"<< endl;
}

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

    // 使用循环生成元素并检查差是否可逆
    for (long i = 0; i < length; i++) {
        // 生成一个元素
        ZZ_pE element;
        random(element);
        // 将元素添加到向量中
        result[i] = element;

        // 检查向量中之前的元素与当前元素的差是否可逆
        for (long j = 0; j < i; j++) {
            if (!Inv(result[i] - result[j])) {
                cout << j << " " << i <<endl;
                // 如果差不可逆，重新生成当前元素并重新检查
                i--;
                break;
            }
        }
    }

    return result;
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

void fillIrred(const ZZ_pX& F, vector<long>& Irred) {
    Irred.clear(); // 清空 Irred

    for (long i = 0; i <= deg(F); i++) {
        long c;
        conv(c,(coeff(F, i)));
        Irred.push_back(c); // 取多项式 F 的第 i 个系数的值
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

vector<long> PadVectorToLength(const vector<long>& input, size_t targetLength) {
    vector<long> paddedVector(targetLength, 0);  // 初始化长度为targetLength的全0向量

    // 将input的内容复制到paddedVector中，不超过targetLength
    size_t copyLength = min(input.size(), static_cast<size_t>(targetLength));
    copy(input.begin(), input.begin() + copyLength, paddedVector.begin());

    return paddedVector;
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

    ZZ_p::init(p);
    ZZ_pX F;
    BuildIrred(F,s);

    fillIrred(F,Irred);

    ZZ_p::init(q3);
    ZZ_pE::init(F);

    vec_ZZ_pE v1 = GenerateInvertibleDifferenceVector(n);

    if(!VectorInvertible(v1)){
        cout<< "Interpolation points are not invertible!\n";
        return;
    }
    cout<<"--interpolation points--"<< v1 << '\n';
    
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

    ZZ_p::init(q3);
    ZZ_pX F = long2ZZpX(Irred);
    ZZ_pE::init(F);

    Vec<ZZ_pEX> V;

    vec_ZZ_pE v1 ;
    v1.SetLength(n);
    for(int i = 0; i < n; i++){
        vector<long> coeff(Interpolation.begin()+i*s, Interpolation.begin()+i*s+s);
        v1[i] = long2ZZpE(coeff);
    }
    cout<<v1<<endl;
    if(!VectorInvertible(v1)){
        cout<< "Interpolation points are not invertible!\n";
        return;
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
        cout<<v2<<endl;
        ZZ_pEX h = interpolate(v1,v2);
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
    if(D*(n-1)>r/s){
        cout<< "Note:D(n-1)<=r/s\n";
        return;
    }
    if(power(p,s)<n){
        cout<< "Lenstra Error\n";
        return;
    }


    vec_ZZ_pE v1 ;
    v1.SetLength(n);
    for(int i = 0; i < n; i++){
        vector<long> coeff(Interpolation.begin()+i*s, Interpolation.begin()+i*s+s);
        v1[i] = long2ZZpE(coeff);
    }

    if(!VectorInvertible(v1)){
        cout<< "Interpolation points are not invertible!\n";
        return;
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

    cout<<res<<endl;

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
    print(Input2);
    RMFE_GR_PHI(ZZ(p),k,n2,r2,r1,D,Interpolation2,Irred2,Input2,result);
}

int main(){
    //Degree D - GR(p^k,s)^n -> GR(p^k,r)  r/s>=D(n-1)
    // long p = 2;
    // long k = 2;
    // long n = 2;
    // long r = 9;
    // long s = 1;
    // long D = 4;

    // vector<long> Interpolation;
    // vector<long> Irred;

    // RMFE_GR_Init(ZZ(p),k,n,r,s,D,Interpolation,Irred);
    // print(Interpolation);
    // print(Irred);
    
    // vector<long> result;
    // vector<long> Input = {1,2,2,3,3,0,1,3};
    // RMFE_GR_PHI(ZZ(p),k,n,r,s,D,Interpolation,Irred,Input,result);
    // print(result);

    // vector<long> Output;
    // RMFE_GR_PSI(ZZ(p),k,n,r,s,D,Interpolation,Irred,result,Output);
    // print(Output);


    long p = 2;
    long k = 2;
    long n1 = 2;
    long r1 = 9;
    long s = 1;
    long D = 4;
    long n2 = 5;
    long r2 = 21*9;

    vector<long> Interpolation1,Interpolation2;
    vector<long> Irred1,Irred2;
    RMFE_GR_com_Init(ZZ(p),k,n1,n2,r1,r2,s,D,Interpolation1,Interpolation2,Irred1,Irred2);

    vector<long> result;
    vector<long> Input = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1};
    RMFE_GR_com_PHI(ZZ(p),k,n1,n2,r1,r2,s,D,Interpolation1,Interpolation2,Irred1,Irred2,Input,result);

    // long p = 2;
    // long k = 2;
    // long n = 2;
    // long r = 12;
    // long s = 3;
    // long D = 2;





    return 0;
}