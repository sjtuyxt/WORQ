//LSMR
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
data tmp[B], tmp2[B];
data itmp[M][B], tmp1[M * 100][B];
char str[100];
int MCnt = 0, LSMCnt = 0;
long long LSMNum[100], LSMEnd[100], DataEnd[100];

bool Check(data x, data start, data end){
    for (int i = 0; i < D; i++)
        if (x.x[i] < start.x[i] || x.x[i] > end.x[i])
            return false;
    return true;
}

int CheckR(data x, data y, data start, data end){
    int ans = 2;
    for (int i = 0; i < D; i++) {
        if (y.x[i] < start.x[i] || x.x[i] > end.x[i])
            return 0;
        if (x.x[i] <= end.x[i] || y.x[i] >= start.x[i])
            ans = 1;
    }
    return ans;
}

long long *HilbertIndexTransposed(long long *hilbertAxes, long long bits)
{
    long long *X = new long long[D];
    memcpy(X, hilbertAxes, sizeof(long long) * D);
    int n = D; // n: Number of dimensions
    long long M = 1 << (bits - 1);
    long long P, Q, t;
    int i;
    // Inverse undo
    for (Q = M; Q > 1; Q >>= 1)
    {
        P = Q - 1;
        for (i = 0; i < n; i++)
            if ((X[i] & Q) != 0)
                X[0] ^= P; // invert
            else {
                t = (X[0] ^ X[i]) & P;
                X[0] ^= t;
                X[i] ^= t;
            }
    } // exchange
    // Gray encode
    for (i = 1; i < n; i++)
        X[i] ^= X[i - 1];
    t = 0;
    for (Q = M; Q > 1; Q >>= 1)
        if ((X[n - 1] & Q)!=0)
            t ^= Q - 1;
    for (i = 0; i < n; i++)
        X[i] ^= t;

    return X;
}

long long F(long long x, long long y){
    return x ^ y;
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

bool lcmp(const data& x, const data& y) {
    long long tmp[D];
    for (int i = 0; i < D; i++)
        tmp[i] = DoubleToLong(x.x[i]);
    long long *X = HilbertIndexTransposed(tmp, 60);
    for (int i = 0; i < D; i++)
        tmp[i] = DoubleToLong(y.x[i]);
    long long *Y = HilbertIndexTransposed(tmp, 60);
    for (int i = 0; i < D; i++) {
        if (X[i] < Y[i]) {
            delete(X);
            delete(Y);
            return true;
        }
        else if (X[i] > Y[i]) {
            delete(X);
            delete(Y);
            return false;
        }
    }
    delete(X);
    delete(Y);
    return false;
}

bool gcmp(const data& x, const data& y) {
    return !lcmp(x, y);
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
//    x = 1616491169;
    srand(x);
    cout << x << endl;
    FILE *fp;
    fp = fopen("DB.txt","w+");

    for (long long p = 0; p < n / B; p++) {
        for (int i = 0; i < B; i++) {
            for (int j = 0; j < D; j++)
                tmp[i].x[j] = rand() / 10000.0;
            tmp[i].w = rand();
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

void ExSort(long long start, long long end, FILE *fp){
    FILE *t1fp, *t2fp;
    t1fp = fp;
    t2fp = fopen("temp.txt", "w+");
    long long offset1 = start;
    long long offset2 = 0;

//    cout << start << "\t" << end << endl;
    int cnt = 0;
    for (long long i = start; i < end; i += B) {
        Read(i, fp, itmp[cnt++]);
        if (cnt == M) {
            sort(itmp[0], itmp[cnt - 1] + B, lcmp);
            for (int j = 0; j < cnt; j++)
                Write(i - B * (cnt - j - 1), t1fp, itmp[j]);
            cnt = 0;
        }
    }

//    cout << start << "\t" << end << endl;

    if (cnt) {
        sort(itmp[0], itmp[cnt - 1] + B, lcmp);
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
            make_heap(Heap, Heap + cand, gcmp);
            int cnt = 0;
            long long pos = st;
            while (cand > 0) {
                int minp = Heap[0].w;
                pop_heap(Heap, Heap + cand, gcmp);
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
                    push_heap(Heap, Heap + cand, gcmp);
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

void Merge(int tree, long long start, long long end){
    FILE *fp;
    sprintf(str, "%d.txt", tree);
    fp = fopen(str, "r+");
    long long head = DataEnd[tree];
    Read(start, fp, tmp);
    int cnt = 0;
    if (tmp[0].add == -1){
        for (long long i = start; i < end; i += B){
            Read(i, fp, tmp);
            data x, y;
            for (int j = 0; j < D; j++) {
                x.x[j] = 1e100;
                y.x[j] = -1e100;
            }
            x.add = i;
            y.add = i;
            y.w = 0;
            for (int j = 0; j < B; j++){
                y.w = F(y.w, tmp[j].w);
                for (int k = 0; k < D; k++){
                    x.x[k] = min(x.x[k], tmp[j].x[k]);
                    y.x[k] = max(y.x[k], tmp[j].x[k]);
                }
            }
            tmp2[cnt] = x;
            tmp2[cnt + 1] = y;
            cnt += 2;
            if (cnt == B){
                cnt = 0;
                Write(DataEnd[tree], fp, tmp2);
                DataEnd[tree] += B;
            }
        }
        if (cnt){
            for (int i = cnt; i < B; i++)
                tmp2[i].w = -1;
            Write(DataEnd[tree], fp, tmp2);
            DataEnd[tree] += B;
        }
        if (DataEnd[tree] - head > B)
            Merge(tree, head, DataEnd[tree]);
    }
    else {
        for (long long i = start; i < end; i += B){
            Read(i, fp, tmp);
            data x, y;
            for (int j = 0; j < D; j++) {
                x.x[j] = 1e100;
                y.x[j] = -1e100;
            }
            x.add = i;
            y.add = i;
            y.w = 0;
            for (int j = 0; j < B; j += 2){
                if (tmp[j + 1].w == -1)
                    break;
                y.w = F(y.w, tmp[j + 1].w);
                for (int k = 0; k < D; k++){
                    x.x[k] = min(x.x[k], tmp[j].x[k]);
                    y.x[k] = max(y.x[k], tmp[j + 1].x[k]);
                }
            }
            tmp2[cnt] = x;
            tmp2[cnt + 1] = y;
            cnt += 2;
            if (cnt == B){
                cnt = 0;
                Write(DataEnd[tree], fp, tmp2);
                DataEnd[tree] += B;
            }
        }
        if (cnt){
            for (int i = cnt; i < B; i++)
                tmp2[i].w = -1;
            Write(DataEnd[tree], fp, tmp2);
            DataEnd[tree] += B;
        }
        if (DataEnd[tree] - head > B)
            Merge(tree, head, DataEnd[tree]);
    }
    fclose(fp);
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
        ExSort(0, length * M * B, fp);
        fclose(fp);
        Merge(dst, 0, length * M * B);
    }
}


long long Query(int tree, long long node, data start, data end, int d){
//    cout << endl << "Query: " << tree << " " << node << " " << d << endl;
    FILE *fp;
    sprintf(str, "%d.txt", tree);
    fp = fopen(str, "r+");
    long long ans = 0;

    Read(node, fp, tmp1[d]);

//    for (int i  = 0; i < B; i++){
//        for (int j = 0; j < D; j++)
//            cout << tmp1[d][i].x[j] << "\t";
//        cout << "|\t" << tmp1[d][i].w << "\t" << tmp1[d][i].add << endl;
//    }

    if (tmp1[d][0].add == -1){
        for (int i = 0; i < B; i++)
            if (Check(tmp1[d][i], start, end))
                ans = F(ans, tmp1[d][i].w);
    }
    else {
        for (int i = 0; i < B; i += 2){
            if (tmp1[d][i].w == -1)
                break;
            int status = CheckR(tmp1[d][i], tmp1[d][i + 1], start, end);
//            cout << status << endl;
            if (status == 1)
                ans = F(ans, Query(tree, tmp1[d][i + 1].add, start, end, d + 1));
            else if (status == 2)
                ans = F(ans, tmp1[d][i + 1].w);
        }
    }
//    cout  << "Query: " << tree << " " << node << " " << d  << endl;
//    cout << "Ans: " << ans << endl;

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
            if (LSMNum[j]) {
                anss = F(anss, Query(j, DataEnd[j] - B, x, y, 0));
            }
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