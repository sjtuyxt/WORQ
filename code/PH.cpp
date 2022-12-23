//PH
#include<cstdio>
#include<ctime>
#include<cstdlib>
#include<cmath>
#include <chrono>
#include "iomanip"
#include<iostream>
#include<fstream>
#include<algorithm>
#include<unordered_map>
#define _FILE_OFFSET_BITS 64
using namespace std;
using namespace chrono;

const int D = 3;
const int B = 200;
const int M = 50;
const int P = 1 << D;
const long long n = 60000;
const int CacheSize = 1024;

struct data{
    double x[D];
    long long w, add;
};

unsigned char buffer[B * (D + 2)][8], bufferP[512][8];
data tmp[B];
long long itmp[65][P], tmp2[P], tmp3[P], prefix[D], tmp4[65][D], tmp5[P];
char str[100];
long long DataEnd = P;
long long Cache[CacheSize][512];
long long CacheCnt[CacheSize], CacheStart[CacheSize];
unordered_map <long long, int> Map;


bool Check(data x, data start, data end){
    for (int i = 0; i < D; i++)
        if (x.x[i] < start.x[i] || x.x[i] > end.x[i])
            return false;
    return true;
}

long long F(long long x, long long y){
    return x + y;
}

long long DoubleToLong(double x){
    long long y;
    memcpy(&y, &x, 8);
    return y;
}

double LongToDouble(long long x){
    double y;
    memcpy(&y, &x, 8);
    return y;
}

void LongToChar(long long x, unsigned char *ch){

    for (int i = 7; i >= 0; i--) {
        ch[i] = (int)(x & 255);
        x >>= 8;
    }
}

long long CharToLong(const unsigned char *ch){
    long long x = 0;
    for (int i = 0; i < 8; i++) {
        x = (x << 8) + ch[i];
    }
    return x;
}

void Read(long long start, FILE *fp, data *input){
    fseeko(fp, start * (D + 2) * 8, SEEK_SET);
    fread(buffer, sizeof(char), (D + 2) * B * 8, fp);
    for (int i = 0; i < B; i++){
        for (int j = 0; j < D; j++){
            input[i].x[j] = LongToDouble(CharToLong(buffer[i * (D + 2) + j]));
        }
        input[i].w = CharToLong(buffer[i * (D + 2) + D]);
        input[i].add = CharToLong(buffer[i * (D + 2) + D + 1]);
    }
}

void Write(long long start, FILE *fp, data *output){
    for (int i = 0; i < B; i++){
        for (int j = 0; j < D; j++){
            LongToChar(DoubleToLong(output[i].x[j]), buffer[i * (D + 2) + j]);
        }
        LongToChar(output[i].w, buffer[i * (D + 2) + D]);
        LongToChar(output[i].add, buffer[i * (D + 2) + D + 1]);
    }
    fseeko(fp, start * (D + 2) * 8, SEEK_SET);
    fwrite(buffer, sizeof(char), (D + 2) * B * 8, fp);
    fflush(fp);
}

long long CacheA = 1;
int CacheB = 19;

int Get(long long start, FILE *fp){
    int x = 0;
    auto ans = Map.find(start);
    if (ans == Map.end()) {
        if (Map.size() < CacheSize) {
            for (int i = 0; i < CacheSize; i++)
                if (CacheCnt[i] == 0) {
                    x = i;
                    break;
                }
        }
        else {
            CacheB++;
            if (CacheB == 20){
                CacheB = 0;
                CacheA++;
            }
            for (int i = 0; i < CacheSize; i++)
                if (CacheCnt[i] < CacheCnt[x])
                    x = i;
            for (int i = 0; i < 512; i++)
                LongToChar(Cache[x][i], bufferP[i]);
            fseeko(fp, CacheStart[x] * 8, SEEK_SET);
            fwrite(bufferP, sizeof(char), 512 * 8, fp);
            Map.erase(CacheStart[x]);
        }
        CacheCnt[x] = CacheA;
        Map[start] = x;
        CacheStart[x] = start;
        fseeko(fp, start * 8, SEEK_SET);
        fread(bufferP, sizeof(char), 512 * 8, fp);
        for (int i = 0; i < 512; i++)
            Cache[x][i] = CharToLong(bufferP[i]);
    }
    else {
        x = ans -> second;
//        cout << "+++" << x << endl;
        CacheCnt[x]++;
    }
    return x;
}

void Dump(FILE *fp){
    for (int i = 0; i < CacheSize; i++)
        if (CacheCnt[i]){
            for (int j = 0; j < 512; j++)
                LongToChar(Cache[i][j], bufferP[j]);
            fseeko(fp, CacheStart[i] * 8, SEEK_SET);
            fwrite(bufferP, sizeof(char), 512 * 8, fp);
            fflush(fp);
        }
}

void ReadP(long long start, FILE *fp, long long *input){
    fseeko(fp, start * 8, SEEK_SET);
    fread(bufferP, sizeof(char), P * 8, fp);
    for (int i = 0; i < P; i++)
        input[i] = CharToLong(bufferP[i]);
}

void WriteP(long long start, FILE *fp, long long *output){
    for (int i = 0; i < P; i++)
        LongToChar(output[i], bufferP[i]);
    fseeko(fp, start * 8, SEEK_SET);
    fwrite(bufferP, sizeof(char), P * 8, fp);
    fflush(fp);
}


void RandomDB(long long n){
    int x = time(NULL);
    x = 1617214389;
    srand(x);
    cout << x << endl;
    FILE *fp;
    fp = fopen("DB.txt","w+");

    for (long long p = 0; p < n / B; p++) {
        for (int i = 0; i < B; i++) {
            for (int j = 0; j < D; j++)
                tmp[i].x[j] = rand() / 10000.0;
            tmp[i].w = rand() % 100 + 1;
            tmp[i].add = -1;
        }
        Write(p * B, fp, tmp);
    }
    fclose(fp);
}

void SaveDB(){
    ofstream out;
    out.open("DBvar.txt");

    out << DataEnd << endl;

    out.close();
}

void LoadDB(){
    ifstream in;
    in.open("DBvar.txt");

    in >> DataEnd;

    in.close();
}

int Order(data input, int pos){
    int ans = 0;
    for (int i = 0; i < D; i++){
        ans = ans * 2 + ((DoubleToLong(input.x[i]) >> pos) & 1);
    }
    return ans;
}

void Insert(FILE *fp, data input){
    long long i = 0;
    long long j = -1;
    for (int pos = 62; pos > 0; pos--){
        ReadP(i, fp, tmp2);
        int ord = Order(input, pos);
        if (tmp2[0] >= 0 && tmp2[ord] == 0) {
            bool flag = 1;
            for (int k = 0; k < P; k++)
                if (tmp2[k]) {
                    flag = 0;
                    break;
                }
            if (flag){
                DataEnd -= P - 2;
                tmp2[0] = -DataEnd;
                tmp2[1] = ord;
            }
            else
                tmp2[ord] = DataEnd;
            WriteP(i, fp, tmp2);
            WriteP(DataEnd, fp, tmp3);
            DataEnd += P;
        }
        else if (tmp2[0] < 0 && tmp2[1] != ord){
            tmp3[ord] = DataEnd + P;
            tmp3[tmp2[1]] = -tmp2[0];
            memcpy(tmp2, tmp3, sizeof(long long) * P);
            WriteP(DataEnd, fp, tmp2);
            memset(tmp3, 0, sizeof(long long) * P);
            WriteP(DataEnd + P, fp, tmp3);
            if (j != -1){
                ReadP(j, fp, tmp5);
                if (tmp5[0] < 0)
                    tmp5[0] = -DataEnd;
                else {
                    for (int k = 0; k < P; k++)
                        if (tmp5[k] == i)
                            tmp5[k] = DataEnd;
                }
                WriteP(j, fp, tmp5);
            }
            DataEnd += 2 * P;
        }
        j = i;
        if (tmp2[0] < 0)
            i = -tmp2[0];
        else
            i = tmp2[ord];
    }
    ReadP(i, fp, tmp2);
    int ord = Order(input, 0);
    if (tmp2[0] >= 0){
        bool flag = 1;
        for (int k = 0; k < P; k++)
            if (tmp2[k]) {
                flag = 0;
                break;
            }
        if (flag){
            DataEnd -= P - 2;
            tmp2[0] = -input.w;
            tmp2[1] = ord;
        }
        else
            tmp2[ord] = F(tmp2[ord], input.w);
        WriteP(i, fp, tmp2);
    }
    else {
        if (tmp2[1] == ord) {
            tmp2[0] = -F(-tmp2[0], input.w);
            WriteP(i, fp, tmp2);
        }
        else {
            tmp3[ord] = input.w;
            tmp3[tmp2[1]] = -tmp2[0];
            WriteP(DataEnd, fp, tmp3);
            tmp3[ord] = 0;
            tmp3[tmp2[1]] = 0;
            ReadP(j, fp, tmp5);
            if (tmp5[0] < 0)
                tmp5[0] = -DataEnd;
            else
                for (int k = 0; k < P; k++)
                    if (tmp5[k] == i)
                        tmp5[k] = DataEnd;
            WriteP(j, fp, tmp5);
            DataEnd += P;
        }
    }
}

bool CheckIn(int pos, int ord, data left, data right, long long bac) {
    for (int i = 0; i < D; i++)
        if ((1 << (D - 1 - i)) & ord)
            prefix[i] += (1LL << pos);
    for (int i = 0; i < D; i++)
        if (prefix[i] < (DoubleToLong(left.x[i]) & bac) || prefix[i] > (DoubleToLong(right.x[i]) & bac))
            return false;
    return true;
}

long long Query(FILE *fp, int pos, long long node, data left, data right, long long bac){
    long long ans = 0;
    ReadP(node, fp, itmp[pos]);

    if (itmp[pos][0] < 0){
        memcpy(tmp4[pos], prefix, sizeof(long long) * D);
        if (CheckIn(pos, itmp[pos][1], left, right, bac)) {
            if (pos > 0)
                ans = Query(fp, pos - 1, -itmp[pos][0], left, right, bac | (1LL << (pos - 1)));
            else
                ans = -itmp[pos][0];
        }
        memcpy(prefix, tmp4[pos], sizeof(long long) * D);
        return ans;
    }
    for (int i = 0; i < P; i++)
        if (itmp[pos][i]) {
            memcpy(tmp4[pos], prefix, sizeof(long long) * D);
            if (CheckIn(pos, i, left, right, bac)) {
                if (pos > 0)
                    ans = F(ans, Query(fp, pos - 1, itmp[pos][i], left, right, bac | (1LL << (pos - 1))));
                else
                    ans = F(ans, itmp[pos][i]);
            }
            memcpy(prefix, tmp4[pos], sizeof(long long) * D);
        }
    return ans;
}

int main(){
    FILE *fp;
    fp = fopen("DB.txt","r+");

    FILE *tfp;
    tfp = fopen("1.txt","w+");

    ofstream result;
    result.open("result.txt");

    WriteP(0, tfp, tmp3);

    int cnt = 0;

    auto start = system_clock::now();
    cout << fixed << setprecision(2);
    for (long long i = 0; i < n / B; i++){
        Read(i * B, fp, tmp);
        for (int j = 0; j < B; j++)
            Insert(tfp, tmp[j]);
        while (n / 100 * cnt < (i + 1) * B){
            cnt++;
            cout << cnt << "% insertion finished" << endl;
            auto end  = system_clock::now();
            auto duration = duration_cast<microseconds>(end - start);
            result << fixed << setprecision(10) << double(duration.count()) * microseconds::period::num / microseconds::period::den << endl;
        }
    }
    SaveDB();

//    LoadDB();

    auto end   = system_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << fixed << setprecision(10) << double(duration.count()) * microseconds::period::num / microseconds::period::den << endl;


    fp = fopen("query.txt","r+");

    for (int i = 0; i < 1000; i++) {
        data x, y;
        for (int j = 0; j < D; j++){
            fscanf(fp, "%lf%lf", &x.x[j], &y.x[j]);
        }
        start = system_clock::now();
        long long anss;
        anss = F(anss, Query(tfp, 62, 0, x, y, 1LL << 62));

        end   = system_clock::now();
        duration = duration_cast<microseconds>(end - start);
        result << fixed << setprecision(10) << double(duration.count()) * microseconds::period::num / microseconds::period::den << endl;
    }

    fclose(fp);
    fclose(tfp);
    result.close();

    return 0;
}