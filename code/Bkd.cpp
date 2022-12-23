//Bkd
#include<cstdio>
#include<ctime>
#include<cstdlib>
#include<cmath>
#include <chrono>
#include "iomanip"
#include<iostream>
#include<fstream>
#include<algorithm>
#define _FILE_OFFSET_BITS 64
using namespace std;
using namespace chrono;

const int D = 3;
const int B = 200;
const int M = 50;
const long long n = 60000;

struct data{
    double x[D];
    long long w, add;
};

unsigned char buffer[B * (D + 2)][8];
data tmp[B], tmp2[B], tmp3[B];
data itmp[M][B], tmp1[M][B];
char str[100];
int MCnt = 0, LSMCnt = 0;
long long LSMNum[100], LSMEnd[100], DataEnd[100];
FILE *Fdim[20];

double Sum(const data& x){
    double sum = 0;
    for (int i = 0; i < D; i++)
        sum += x.x[i];
    sum += x.w;
    return sum;
}

bool l0(const data& x, const data& y) {return x.x[0] < y.x[0] || (x.x[0] == y.x[0] && Sum(x) < Sum(y));}
bool l1(const data& x, const data& y) {return x.x[1] < y.x[1] || (x.x[1] == y.x[1] && Sum(x) < Sum(y));}
bool l2(const data& x, const data& y) {return x.x[2] < y.x[2] || (x.x[2] == y.x[2] && Sum(x) < Sum(y));}
bool l3(const data& x, const data& y) {return x.x[3] < y.x[3] || (x.x[3] == y.x[3] && Sum(x) < Sum(y));}
bool l4(const data& x, const data& y) {return x.x[4] < y.x[4] || (x.x[4] == y.x[4] && Sum(x) < Sum(y));}
bool l5(const data& x, const data& y) {return x.x[5] < y.x[5] || (x.x[5] == y.x[5] && Sum(x) < Sum(y));}
bool l6(const data& x, const data& y) {return x.x[6] < y.x[6] || (x.x[6] == y.x[6] && Sum(x) < Sum(y));}
bool l7(const data& x, const data& y) {return x.x[7] < y.x[7] || (x.x[7] == y.x[7] && Sum(x) < Sum(y));}
bool l8(const data& x, const data& y) {return x.x[8] < y.x[8] || (x.x[8] == y.x[8] && Sum(x) < Sum(y));}
bool l9(const data& x, const data& y) {return x.x[9] < y.x[9] || (x.x[9] == y.x[9] && Sum(x) < Sum(y));}
bool l10(const data& x, const data& y) {return x.x[10] < y.x[10]  || (x.x[10] == y.x[10] && Sum(x) < Sum(y));}
bool (*cmp[11])(const data&, const data&) = {l0, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10};

bool g0(const data& x, const data& y) {return x.x[0] > y.x[0] || (x.x[0] == y.x[0] && Sum(x) > Sum(y));}
bool g1(const data& x, const data& y) {return x.x[1] > y.x[1] || (x.x[1] == y.x[1] && Sum(x) > Sum(y));}
bool g2(const data& x, const data& y) {return x.x[2] > y.x[2] || (x.x[2] == y.x[2] && Sum(x) > Sum(y));}
bool g3(const data& x, const data& y) {return x.x[3] > y.x[3] || (x.x[3] == y.x[3] && Sum(x) > Sum(y));}
bool g4(const data& x, const data& y) {return x.x[4] > y.x[4] || (x.x[4] == y.x[4] && Sum(x) > Sum(y));}
bool g5(const data& x, const data& y) {return x.x[5] > y.x[5] || (x.x[5] == y.x[5] && Sum(x) > Sum(y));}
bool g6(const data& x, const data& y) {return x.x[6] > y.x[6] || (x.x[6] == y.x[6] && Sum(x) > Sum(y));}
bool g7(const data& x, const data& y) {return x.x[7] > y.x[7] || (x.x[7] == y.x[7] && Sum(x) > Sum(y));}
bool g8(const data& x, const data& y) {return x.x[8] > y.x[8] || (x.x[8] == y.x[8] && Sum(x) > Sum(y));}
bool g9(const data& x, const data& y) {return x.x[9] > y.x[9] || (x.x[9] == y.x[9] && Sum(x) > Sum(y));}
bool g10(const data& x, const data& y) {return x.x[10] > y.x[10] || (x.x[10] == y.x[10] && Sum(x) > Sum(y));}
bool (*gcmp[11])(const data&, const data&) = {g0, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10};

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


void RandomDB(long long n){
    int x = time(NULL);
    srand(x);
    cout << x << endl;
    FILE *fp;
    fp = fopen("DB.txt","w+");

    for (long long p = 0; p < n / B; p++) {
        for (int i = 0; i < B; i++) {
            for (int j = 0; j < D; j++)
                tmp[i].x[j] = rand() / 10000.0;
            tmp[i].w = rand() % 100;
            tmp[i].add = -1;
        }
        Write(p * B, fp, tmp);
    }
    fclose(fp);
}

void SaveDB(){
    ofstream out;
    out.open("DBvar.txt");

    out << LSMCnt << endl;
    for (int i = 0; i < LSMCnt; i++)
        out << LSMNum[i] << " " << LSMEnd[i] << " " << DataEnd[i] << endl;

    out.close();
}

void LoadDB(){
    ifstream in;
    in.open("DBvar.txt");

    in >> LSMCnt;
    for (int i = 0; i < LSMCnt; i++)
        in >> LSMNum[i] >> LSMEnd[i] >> DataEnd[i];

    in.close();
}

void ExSort(long long start, long long end, int dim, FILE *fp){
    FILE *t1fp, *t2fp;
    t1fp = fp;
    t2fp = fopen("temp.txt", "w+");
    long long offset1 = start;
    long long offset2 = 0;

    int cnt = 0;

    for (long long i = start; i < end; i += B) {
        Read(i, fp, itmp[cnt++]);
        if (cnt == M) {
            sort(itmp[0], itmp[cnt - 1] + B, cmp[dim]);
            for (int j = 0; j < cnt; j++)
                Write(i - B * (cnt - j - 1), t1fp, itmp[j]);
            cnt = 0;
        }
    }

    if (cnt) {
        sort(itmp[0], itmp[cnt - 1] + B, cmp[dim]);
        for (int j = 0; j < cnt; j++)
            Write(end - B * (cnt - j), t1fp, itmp[j]);
    }

    long long length = end - start;
    int head[M];
    long long finished[M];
    data Heap[M];

    for (long long step = B * M; step < length; step *= M) {

        for (long long st = 0; st < length; st += M * step) {
            int cand = 0;
            while (cand < M && st + cand * step < length) {
                head[cand] = 0;
                finished[cand] = 0;
                Read(st + cand * step + offset1, t1fp, itmp[cand]);
                Heap[cand] = itmp[cand][0];
                Heap[cand].w = cand;
                cand++;
            }
            make_heap(Heap, Heap + cand, gcmp[dim]);
            int cnt = 0;
            long long pos = st;
            while (cand > 0) {
                int minp = Heap[0].w;
                pop_heap(Heap, Heap + cand, gcmp[dim]);
                tmp[cnt] = itmp[minp][head[minp]];
                cnt++;
                head[minp]++;

                if (cnt == B) {
                    Write(pos + offset2, t2fp, tmp);
                    cnt = 0;
                    pos += B;
                }

                if (head[minp] == B) {
                    finished[minp]++;
                    if (finished[minp] < step / B && st + minp * step + finished[minp] * B < length) {
                        head[minp] = 0;
                        Read(st + minp * step + finished[minp] * B + offset1, t1fp, itmp[minp]);
                    } else {
                        head[minp] = -1;
                        cand--;
                    }
                }

                if (head[minp] != -1) {
                    Heap[cand - 1] = itmp[minp][head[minp]];
                    Heap[cand - 1].w = minp;
                    push_heap(Heap, Heap + cand, gcmp[dim]);
                }
            }
        }
        swap(t1fp, t2fp);
        swap(offset1, offset2);
    }

    if (t1fp != fp) {
        for (long long i = start; i < end; i += B) {
            Read(i - start, t1fp, tmp);
            Write(i, fp, tmp);
        }
        fclose(t1fp);
    }
    else
        fclose(t2fp);
}

void Merge(int tree, long long start, long long end, int dim){
//    cout << "Merge" << tree << " " << start << " " << end << " " << dim << endl;
    if (end - start <= B)
        return;
    if (start - end <= B * M) {
        FILE *fp;
        fp = Fdim[0];

        ExSort(start, end, dim, fp);
        long long mid = (start + end) / 2;
        dim--;
        if (dim == -1)
            dim = D - 1;
        Merge(tree, start, mid - mid % B, dim);
        Merge(tree, mid - mid % B + B, end, dim);

    }
    else{
        long long mid = (start + end) / 2;
        long long mstart = mid - mid % B;
        Read(mstart, Fdim[dim], tmp);
        data x = tmp[0];
        data y = tmp[B - 1];
        for (int i = 0; i < D; i++)
            if (i != dim){
                long long pos2 = 0;
                long long pos3 = mstart + B;
                int cnt2 = 0;
                int cnt3 = 0;
                Read(mstart, Fdim[dim], tmp);
                Write(mstart, Fdim[D], tmp);
                for (long long j = start; j < end; j += B){
                    Read(j, Fdim[i], tmp);
                    for (int k = 0; k < B; k++){
                        if (cmp[dim](tmp[k], x))
                            tmp2[++cnt2] = tmp[k];
                        else if (gcmp[dim](tmp[k], y))
                            tmp3[++cnt3] = tmp[k];
                        if (cnt2 == B){
                            Write(pos2, Fdim[D], tmp2);
                            pos2 += B;
                            cnt2 = 0;
                        }
                        else if (cnt3 == B){
                            Write(pos3, Fdim[D], tmp3);
                            pos3 += B;
                            cnt3 = 0;
                        }
                    }
                }
                swap(Fdim[D], Fdim[i]);
            }
    }
}

void Insert(data *input){
    memcpy(itmp[MCnt], input, (D + 2) * B * 8);
    MCnt++;

    if (MCnt == M){
        MCnt = 0;
        int dst = 0;
        while (LSMNum[dst]) {
            LSMNum[dst] = 0;
            LSMEnd[dst] = 0;
            DataEnd[dst] = 0;
            dst++;
        }
        if (dst + 1 > LSMCnt)
            LSMCnt = dst + 1;
        LSMNum[dst] = 1;
        FILE *fp;
        long long length = 1;

        sprintf(str, "%d.txt", dst);
        fp = fopen(str, "w+");
        for (int j = 0; j < M; j++)
            Write(j * B, fp, itmp[j]);
        for (int i = 0; i < dst; i++){
            FILE *tfp;
            sprintf(str, "%d.txt", i);
            tfp = fopen(str,"r+");
            for (long long j = 0; j < length * M; j++){
                Read(j * B, tfp, tmp);
                Write(length * M * B + j * B, fp, tmp);
            }
            fclose(tfp);
            tfp = fopen(str,"w+");
            fclose(tfp);
            length *= 2;
        }
        DataEnd[dst] = length * M * B;
        fclose(fp);
        sprintf(str, "%d.txt", dst);
        Fdim[0] = fopen(str, "r+");
        ExSort(0, length * M * B, 0, Fdim[0]);
        for (int i = 1; i <= D; i++) {
            sprintf(str, "D%d.txt", dst);
            Fdim[i] = fopen(str, "w+");
            if (i < D) {
                for (long long j = 0; j < length * M * B; j += B) {
                    Read(j, Fdim[0], tmp);
                    Write(j, Fdim[i], tmp);
                }
                ExSort(0, length * M * B, i, Fdim[i]);
            }
        }
        Merge(dst, 0, length * M * B, D - 1);
        for (int i = 0; i <= D; i++)
            fclose(Fdim[i]);
        for (int i = 1; i <= D; i++) {
            sprintf(str, "D%d.txt", dst);
            Fdim[i] = fopen(str, "w+");
            fclose(Fdim[i]);
        }
    }
}


long long Query(int tree, long long start, long long end, data left, data right, int dim){
    sprintf(str, "%d.txt", tree);
    FILE *fp = fopen(str, "r+");
    long long ans = 0;
    long long mid = (start + end) / 2;
    Read(mid - mid % B, fp, tmp);
    for (int i = 0; i < B; i++) {
        if (Check(tmp[i], left, right))
            ans = F(ans, tmp[i].w);
    }
    int newdim = dim - 1;
    if (newdim == -1)
        newdim = D - 1;
    double line = tmp[mid % B].x[dim];
    if (left.x[dim] <= line && mid - mid % B != start)
        ans = F(ans, Query(tree, start, mid - mid % B, left, right, newdim));
    if (right.x[dim] >= line && mid - mid % B + B != end)
        ans = F(ans, Query(tree, mid - mid % B + B, end, left, right, newdim));
    fclose(fp);
    return ans;
}

int main(){
    FILE *fp;
    fp = fopen("DB.txt","r+");
    ofstream result;
    result.open("result.txt");

    int cnt = 0;

    auto start = system_clock::now();
    for (long long i = 0; i < n / B; i++){
        Read(i * B, fp, tmp);
        Insert(tmp);
        while (n / 100 * cnt < (i + 1) * B){
            cnt++;
            cout << cnt << "% insertion finished" << endl;
            auto end  = system_clock::now();
            auto duration = duration_cast<microseconds>(end - start);
            result << fixed << setprecision(10) << double(duration.count()) * microseconds::period::num / microseconds::period::den << endl;
        }
    }
    SaveDB();
    fclose(fp);

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
        long long anss = 0;
        long long node = M * B;
        for (int j = 0; j < LSMCnt; j++) {
            if (LSMNum[j])
                anss = F(anss, Query(j, 0, node, x, y, D - 1));
            node *= 2;
        }
        end   = system_clock::now();
        duration = duration_cast<microseconds>(end - start);
        result << fixed << setprecision(10) << double(duration.count()) * microseconds::period::num / microseconds::period::den << endl;
    }

    fclose(fp);
    result.close();

    return 0;
}