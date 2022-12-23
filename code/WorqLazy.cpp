//WorqL
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
#pragma GCC optimize(2)
using namespace std;
using namespace chrono;

const int D = 3;
const int P = 200;
const int B = 200;
const int M = 50;
const long long n = 60000;
const long long LazyC = log(n / B) / log(2);

struct data{
    double x[D];
    long long w, add;
};

unsigned char buffer[B * (D + 2)][8];
data tmp[B], tmp2[B], tmp4[B];
data itmp[M][B], tmp1[D][B], tmp3[D][B * 2], tmp5[D][B * 2], tmp6[100][B];
char str[100];
int MCnt = 0, LSMCnt = 0;
long long LSMNum[200], LSMEnd[200], DataEnd[200], RecordEnd[200];
FILE *LSMtree[200], *Segtree[200], *tmpFile;

bool l0(const data& x, const data& y) {return x.x[0] < y.x[0];}
bool l1(const data& x, const data& y) {return x.x[1] < y.x[1];}
bool l2(const data& x, const data& y) {return x.x[2] < y.x[2];}
bool l3(const data& x, const data& y) {return x.x[3] < y.x[3];}
bool l4(const data& x, const data& y) {return x.x[4] < y.x[4];}
bool (*cmp[5])(const data&, const data&) = {l0, l1, l2, l3, l4};

bool g0(const data& x, const data& y) {return x.x[0] > y.x[0];}
bool g1(const data& x, const data& y) {return x.x[1] > y.x[1];}
bool g2(const data& x, const data& y) {return x.x[2] > y.x[2];}
bool g3(const data& x, const data& y) {return x.x[3] > y.x[3];}
bool g4(const data& x, const data& y) {return x.x[4] > y.x[4];}
bool (*gcmp[5])(const data&, const data&) = {g0, g1, g2, g3, g4};

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

long long IOcnt = 0;
long long tmpIO = 0;

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
    IOcnt++;
}


void ReadP(long long start, FILE *fp, data *input) {
    fseeko(fp, start * (D + 2) * 8, SEEK_SET);
    fread(buffer, sizeof(char), (D + 2) * P * 8, fp);
    for (int i = 0; i < P; i++) {
        for (int j = 0; j < D; j++) {
            input[i].x[j] = LongToDouble(CharToLong(buffer[i * (D + 2) + j]));
        }
        input[i].w = CharToLong(buffer[i * (D + 2) + D]);
        input[i].add = CharToLong(buffer[i * (D + 2) + D + 1]);
    }
    IOcnt++;
}

void ReadB(long long start, FILE *fp, data *input) {
    fseeko(fp, start * (D + 2) * 8, SEEK_SET);
    fread(buffer, sizeof(char), (D + 2) * (B - P) * 8, fp);
    for (int i = 0; i < (B - P); i++) {
        for (int j = 0; j < D; j++) {
            input[i].x[j] = LongToDouble(CharToLong(buffer[i * (D + 2) + j]));
        }
        input[i].w = CharToLong(buffer[i * (D + 2) + D]);
        input[i].add = CharToLong(buffer[i * (D + 2) + D + 1]);
    }
    IOcnt++;
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
    IOcnt++;
}


void RandomDB(long long n){
    int x = time(NULL);
    x = 1617640380;
    srand(x);
    cout << x << endl;
    FILE *fp;
    fp = fopen("DB.txt","w+");

    for (long long p = 0; p < n / B; p++) {
        for (int i = 0; i < B; i++) {
            for (int j = 0; j < D; j++)
                tmp[i].x[j] = rand() / 10000.0;
            tmp[i].w = rand() % 4 + 1;
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

    for (int i = 0; i < LSMCnt; i++) {
        sprintf(str, "%d.txt", i);
        LSMtree[i] = fopen(str, "r+");
        sprintf(str, "%dT.txt", i);
        Segtree[i] = fopen(str, "r+");
        if (i != LSMCnt - 1){
            sprintf(str, "%d.txt", i + 100);
            LSMtree[i + 100] = fopen(str, "r+");
            sprintf(str, "%dT.txt", i + 100);
            Segtree[i + 100] = fopen(str, "r+");
        }
    }

    in.close();
}

void ExSort(long long start, long long end, int dim, FILE *fp){
    FILE *t1fp, *t2fp;
    t1fp = fp;
    t2fp = tmpFile;
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

    if (t1fp != fp)
        for (long long i = start; i < end; i += B) {
            Read(i - start, t1fp, tmp);
            Write(i, fp, tmp);
        }
}

long long SegMerge(int tree, long long start, long long end, int dim){
    FILE *fp = LSMtree[tree], * fpT = Segtree[tree];

    ExSort(start, end, dim, fp);

    for (long long i = start; i < end; i += B){
        Read(i, fp, tmp);
        Write(i - start + LSMEnd[tree], fpT, tmp);
    }

    long long i = LSMEnd[tree];
    LSMEnd[tree] += end - start;
    while (i + B < LSMEnd[tree]){
        int cnt = 1;
        Read(i, fpT, tmp);
        tmp2[1] = tmp[B - 1];
        tmp2[0] = tmp[0];
        tmp2[0].add = i;
        for (int j = 1; j < B; j++)
            tmp2[0].w = F(tmp2[0].w, tmp[j].w);
        while (i + cnt * B < LSMEnd[tree] && cnt < B - 1){
            Read(i + cnt * B, fpT, tmp);
            if (tmp[0].x[dim] < tmp2[cnt].x[dim])
                break;
            tmp2[cnt + 1] = tmp[B - 1];
            tmp2[cnt] = tmp[0];
            tmp2[cnt].add = i + cnt * B;
            for (int j = 1; j < B; j++)
                tmp2[cnt].w = F(tmp2[cnt].w, tmp[j].w);
            cnt++;
        }
        for (int j = cnt; j < B; j++) {
            tmp2[j].w = 0;
            tmp2[j].add = -1;
        }
        if (cnt != B - 1){
            for (int j = 0; j < D; j++)
                tmp2[cnt].x[j] = 1e100;
            tmp2[B - 1] = tmp2[cnt];
        }

        Write(LSMEnd[tree], fpT, tmp2);
        LSMEnd[tree] += B;
        i += cnt * B;
    }
    return LSMEnd[tree] - B;
}


long long PartMerge(int tree, long long start, long long end, int dim){
    FILE *fp = LSMtree[tree], * fpT = Segtree[tree];

    if (dim == 1)
        ExSort(start, end, dim, fp);

    long long head = DataEnd[tree];
    int cnt = 0;
    long long length = pow((end - start) / B, 2 / 3.0);
    if (dim == 0)
        length = pow((end - start) / B, 1 / 2.0);
    for (long long i = start; i < end; i += B * length){
        Read(i, fp, tmp2);
        tmp[cnt] = tmp2[0];
        tmp[cnt].w = i;
        cnt++;
        if (cnt == B){
            cnt = 0;
            Write(DataEnd[tree], fp, tmp);
            DataEnd[tree] += B;
        }
    }
    Read(end - B, fp, tmp2);
    tmp[cnt] = tmp2[B - 1];
    tmp[cnt].w = end;
    for (int i = cnt + 1; i < B; i++)
        tmp[i].w = tmp[i].add = -1;
    Write(DataEnd[tree], fp, tmp);
    DataEnd[tree] += B;

    long long ed = DataEnd[tree];
    long long x, y;
    for (long long i = head; i < ed; i += B){
        Read(i, fp, tmp1[dim]);

        for (int j = 0; j < B; j++){
            if (tmp1[dim][j].w == -1)
                break;
            y = tmp1[dim][j].w;
            if (i != head || j != 0){
                tmp1[dim][j].add = SegMerge(tree, x, y, dim ^ 1);
                if (dim == 1)
                    tmp1[dim][j].w = PartMerge(tree, x, y, 0);
            }
            else
                tmp1[dim][j].add = ed;
            x = y;
        }
        Write(i, fp, tmp1[dim]);
    }
    return head;
}

long long WorqMerge(int tree, long long start, long long end, int dim){
    if (dim == 1)
        return PartMerge(tree, start, end ,dim);

    FILE *fp = LSMtree[tree];

    ExSort(start, end, dim, fp);
    long long head = DataEnd[tree];
    int cnt = 0;
    long long length = pow((end - start) / B, (2.0 * dim - 3) / (2.0 * dim - 1)) * 3;

    for (long long i = start; i < end; i += B * length){
        Read(i, fp, tmp2);
        tmp[cnt] = tmp2[0];
        tmp[cnt].w = i;
        cnt++;
        Read(min(i + B * (length - 1), end - B), fp, tmp2);
        tmp[cnt] = tmp2[B - 1];
        tmp[cnt].w = min(i + B * length, end);
        cnt++;
        if (cnt == B){
            cnt = 0;
            Write(DataEnd[tree], fp, tmp);
            DataEnd[tree] += B;
        }
    }
    if (cnt) {
        for (int i = cnt; i < B; i++)
            tmp[i].w = tmp[i].add = -1;
        Write(DataEnd[tree], fp, tmp);
        DataEnd[tree] += B;
    }
    long long ed = DataEnd[tree];
    for (long long i = head; i < ed; i += B){
        Read(i, fp, tmp1[dim]);
        for (int j = 0; j < B; j += 2){
            if (tmp1[dim][j].w == -1)
                break;
            tmp1[dim][j + 1].add = WorqMerge(tree, tmp1[dim][j].w, tmp1[dim][j + 1].w, dim - 1);
            if (i == head && j == 0)
                tmp1[dim][j].add = ed;
        }
        Write(i, fp, tmp1[dim]);
    }

    return head;
}

long long CopySeg(int next, int tree, long long node, int level){
    FILE *fp = Segtree[next], *tfp = Segtree[tree];
    long long ans = LSMEnd[next];
    LSMEnd[next] += B;
    Read(node, tfp, tmp6[level]);
    for (int i = 0; i < B; i++)
        if (tmp6[level][i].add != -1)
            tmp6[level][i].add = CopySeg(next, tree, tmp6[level][i].add, level + 1);
    Write(ans, fp, tmp6[level]);
    return ans;
}

long long CopyMerge(int next, int tree, long long node, int dim){
    long long ans = DataEnd[next];
    FILE *fp = LSMtree[next], *t1fp = LSMtree[tree];
    if (dim == 1){
        long long start = RecordEnd[next];
        Read(node, t1fp, tmp);
        long long ed1 = tmp[0].add;
        long long st1 = DataEnd[next];
        DataEnd[next] += ed1 - node;

        for (long long i = node; i < ed1; i += B){
            Read(i, t1fp, tmp);
            for (int j = 0; j < B; j++) {
                if (tmp[j].w == -1)
                    break;
                if (i != node || j != 0) {
                    long long st2 = tmp[j].w;
                    tmp[j].add = CopySeg(next, tree, tmp[j].add, 0);
                    Read(st2, t1fp, tmp2);
                    long long ed2 = tmp2[0].add;
                    tmp[j].w = DataEnd[next];
                    long long x, y;
                    for (long long k = st2; k < ed2; k += B){
                        Read(k, t1fp, tmp2);
                        for (int p = 0; p < B; p++) {
                            if (tmp2[p].w == -1)
                                break;
                            y = tmp2[p].w;
                            if (k != st2 || p != 0) {
                                tmp2[p].add = CopySeg(next, tree, tmp2[p].add, 0);
                                for (long long q = x; q < y; q += B) {
                                    Read(q, t1fp, tmp4);
                                    Write(RecordEnd[next], fp, tmp4);
                                    RecordEnd[next] += B;
                                }
                            }
                            else
                                tmp2[0].add = DataEnd[next] + ed2 - st2;
                            tmp2[p].w = RecordEnd[next];
                            x = y;
                        }
                        Write(DataEnd[next], fp, tmp2);
                        DataEnd[next] += B;
                    }
                }
                else
                    tmp[0].add = DataEnd[next];
            }
            Write(st1, fp, tmp);
            st1 += B;
        }
    }
    else{
        Read(node, t1fp, tmp);
        long long ed = tmp[0].add;
        for (long long i = node; i < ed; i += B) {
            Read(i, t1fp, tmp);
            if (i == node)
                tmp[0].add = ans + ed - node;
            Write(DataEnd[next], fp, tmp);
            DataEnd[next] += B;
        }
        long long end = DataEnd[next];
        for (long long i = ans; i < end; i += B) {
            Read(i, fp, tmp1[dim]);
            for (int j = 0; j < B; j += 2) {
               if (tmp1[dim][j].w != -1) {
                    tmp1[dim][j].w = RecordEnd[next];
                    tmp1[dim][j + 1].add = CopyMerge(next, tree, tmp1[dim][j + 1].add, dim - 1);
                    tmp1[dim][j + 1].w = RecordEnd[next];
                }
                else
                    break;
            }
            Write(i, fp, tmp1[dim]);
        }
    }
    return ans;
}

long long LazyMerge(int next, int ltree, int rtree, long long lnode, long long rnode, int dim, long long length){
    long long ans = DataEnd[next];
    FILE *fp = LSMtree[next], *t1fp = LSMtree[ltree], *t2fp = LSMtree[rtree];

    if (dim == 1){
        long long start = RecordEnd[next];

        FILE *tfp = t1fp;
        long long node = lnode;
        for (int I = 0; I < 2; I++){
            Read(node, tfp, tmp);
            long long ed1 = tmp[0].add;
            for (long long i = node; i < ed1; i += B){
                Read(i, tfp, tmp);
                for (int j = 0; j < B; j++) {
                    if (tmp[j].w == -1)
                        break;
                    if (i != node || j != 0) {
                        Read(tmp[j].w, tfp, tmp2);
                        long long ed2 = tmp2[0].add;
                        long long x, y;
                        for (long long k = tmp[j].w; k < ed2; k += B){
                            Read(k, tfp, tmp2);
                            for (int p = 0; p < B; p++) {
                                if (tmp2[p].w == -1)
                                    break;
                                y = tmp2[p].w;
                                if (k != tmp[j].w || p != 0)
                                    for (long long q = x; q < y; q += B) {
                                        Read(q, tfp, tmp4);
                                        Write(RecordEnd[next], fp, tmp4);
                                        RecordEnd[next] += B;
                                    }
                                x = y;
                            }
                        }
                    }
                }
            }
            tfp = t2fp;
            node = rnode;
        }

        PartMerge(next, start, RecordEnd[next], 1);
    }
    else {
        length = pow(length, (2.0 * dim - 3) / (2.0 * dim - 1)) * 3;

//        cout << endl << length << endl;

        Read(lnode, t1fp, tmp);
        long long led = tmp[0].add;
        Read(rnode, t2fp, tmp2);
        long long red = tmp2[0].add;
        long long lpos = lnode, rpos = rnode;
        long long i = 0, j = 0;
        int cnt = 0;

        while (lpos < led || rpos < red){
            data cand;
            cand.x[dim] = 1e30;
            int from;
            if (lpos != led && tmp[i].x[dim] < cand.x[dim]){
                cand = tmp[i];
                from = 1;
            }
            if (rpos != red && tmp2[j].x[dim] < cand.x[dim]){
                cand = tmp2[j];
                from = 2;
            }
            if (from == 1){
                tmp4[cnt] = tmp[i];
                tmp4[cnt].add = -1;
                tmp4[cnt + 1] = tmp[i + 1];
                cnt += 2;
                i += 2;
                if (i == B){
                    i = 0;
                    lpos += B;
                    if (lpos < led){
                        Read(lpos, t1fp, tmp);
                    }
                }
                if (tmp[i].w == -1)
                    lpos = led;
            }

            if (from == 2){
                tmp4[cnt] = tmp2[j];
                tmp4[cnt].add = -2;
                tmp4[cnt + 1] = tmp2[j + 1];
                cnt += 2;
                j += 2;
                if (j == B){
                    j = 0;
                    rpos += B;
                    if (rpos < red){
                        Read(rpos, t2fp, tmp2);
                    }
                }
                if (tmp2[j].w == -1)
                    rpos = red;
            }

            if (cnt == B){
                cnt = 0;
                Write(DataEnd[next], fp, tmp4);
                DataEnd[next] += B;
            }
        }
        if (cnt) {
            for (int i = cnt; i < B; i++)
                tmp4[i].w = -1;
            Write(DataEnd[next], fp, tmp4);
            DataEnd[next] += B;
            cnt = 0;
        }

        long long ed = DataEnd[next];
        long long pos = ans;
        cnt = 0;
        j = B;
        Read(ans, fp, tmp3[dim] + B);
        for (long long i = ans; i < ed; i += B){
            memcpy(tmp3[dim], tmp3[dim] + B, sizeof(data) * B);
            j -= B;
            if (i + B < ed)
                Read(i + B, fp, tmp3[dim] + B);
            else
                tmp3[dim][B].w = -1;
            for (; j < B; j += 2){

//                cout << (tmp3[dim][j + 1].w - tmp3[dim][j].w) / B << " ";

                if (tmp3[dim][j].w == -1)
                    break;
                if ((tmp3[dim][j + 1].w - tmp3[dim][j].w) / B > length || tmp3[dim][j + 2].w == -1 || (tmp3[dim][j + 3].w - tmp3[dim][j + 2].w) / B > length){
                    tmp5[dim][cnt] = tmp3[dim][j];
                    tmp5[dim][cnt + 1] = tmp3[dim][j + 1];
                    tmp5[dim][cnt].w = RecordEnd[next];
                    int ntree = tmp5[dim][cnt].add == -1 ? ltree : rtree;
                    tmp5[dim][cnt + 1].add = CopyMerge(next, ntree, tmp5[dim][cnt + 1].add, dim - 1);
                    tmp5[dim][cnt + 1].w = RecordEnd[next];
                    cnt += 2;
                }
                else {
                    tmp5[dim][cnt] = tmp3[dim][j];
                    tmp5[dim][cnt + 1] = tmp3[dim][j + 3];
                    if (tmp3[dim][j+1].x[dim] > tmp3[dim][j + 3].x[dim])
                        tmp5[dim][cnt + 1] = tmp3[dim][j + 1];
                    tmp5[dim][cnt].w = RecordEnd[next];
                    int n1tree = tmp3[dim][j].add == -1 ? ltree : rtree;
                    int n2tree = tmp3[dim][j + 2].add == -1 ? ltree : rtree;
                    tmp5[dim][cnt + 1].add = LazyMerge(next, n1tree, n2tree, tmp3[dim][j + 1].add, tmp3[dim][j + 3].add, dim - 1, (tmp3[dim][j + 1].w - tmp3[dim][j].w + tmp3[dim][j + 3].w - tmp3[dim][j + 2].w) / B);
                    tmp5[dim][cnt + 1].w = RecordEnd[next];
                    j += 2;
                    cnt += 2;
                }

                int lazycnt = 0;
                int lazystart = 0;
                for (int k = cnt - 4; k >= 0 && cnt - k <= 8 * LazyC; k -= 2)
                    if (tmp5[dim][k].x[dim] < tmp5[dim][cnt - 1].x[dim] && tmp5[dim][k + 1].x[dim] > tmp5[dim][cnt - 2].x[dim]){
                        lazycnt++;
                        lazystart = k;
                    }

                if (lazycnt >= LazyC){
                    int size = (cnt - lazystart) / 2;
                    long long st = tmp5[dim][lazystart].w;
                    long long ed = tmp5[dim][cnt - 1].w;
                    ExSort(st, ed, dim, fp);
                    long long newlength = ((ed - st) / B) / size;
                    for (long long k = st; k < ed && lazystart != cnt; k += newlength * B){
                        Read(k, fp, tmp);
                        tmp5[dim][lazystart] = tmp[0];
                        tmp5[dim][lazystart].w = k;
                        long long ned = k + newlength * B - B;
                        if (lazystart == cnt - 2)
                            ned = ed - B;
                        Read(ned, fp, tmp);
                        tmp5[dim][lazystart + 1] = tmp[B - 1];
                        tmp5[dim][lazystart + 1].w = ned + B;
                        tmp5[dim][lazystart + 1].add = WorqMerge(next, k, ned + B, dim - 1);
                        lazystart += 2;
                    }
                }

                if (cnt == B * 2){
                    Write(pos, fp, tmp5[dim]);
                    pos += B;
                    cnt = B;
                    memcpy(tmp5[dim], tmp5[dim] + B, sizeof(data) * B);
                }
            }
        }
        for (int k = cnt; k < 2 * B; k++)
            tmp5[dim][k].w = -1;
        Write(pos, fp, tmp5[dim]);
        pos += B;
        if (cnt > B) {
            Write(pos, fp, tmp5[dim] + B);
            pos += B;
        }
        Read(ans, fp, tmp);
        tmp[0].add = pos;
        Write(ans, fp, tmp);


//        cout << endl;
    }
    return ans;
}

void Clear(int tree){
    if (LSMNum[tree]) {
        fclose(LSMtree[tree]);
        fclose(Segtree[tree]);
    }
    sprintf(str, "%d.txt", tree);
    LSMtree[tree] = fopen(str, "w+");
    sprintf(str, "%dT.txt", tree);
    Segtree[tree] = fopen(str, "w+");
}

void Insert(data *input){
    memcpy(itmp[MCnt], input, (D + 2) * B * 8);
    MCnt++;

    if (MCnt == M){
        MCnt = 0;
        int dst = 0;
        while (LSMNum[dst])
            dst++;
        if (dst + 1 > LSMCnt)
            LSMCnt = dst + 1;
        long long length = 1;

        int next = 100;
        if (dst == 0)
            next = 0;
        Clear(next);
        FILE *fp = LSMtree[next];
        for (int j = 0; j < M; j++)
            Write(j * B, fp, itmp[j]);
        DataEnd[next] = M * B;
        WorqMerge(next, 0, M * B, D - 1);
        LSMNum[next] = 1;

        for (int i = 0; i < dst; i++){
            int next = i + 101;
            if (i == dst  - 1)
                next = dst;
            Clear(next);

            DataEnd[next] = M * B * length * 2;
            LazyMerge(next, i, i + 100, length * B * M, length * B * M, D - 1, length * M * 2);
            LSMNum[next] = 1;

            Clear(i);
            Clear(i + 100);
            length *= 2;
            LSMEnd[i] = 0;
            LSMEnd[i + 100] = 0;
            RecordEnd[i] = 0;
            RecordEnd[i + 100] = 0;
            LSMNum[i] = 0;
            LSMNum[i + 100] = 0;
        }
    }
}

long long SegQuery(FILE *tree, long long node, double start, double end, int dim) {
    ReadP(node, tree, tmp);

    if (end < tmp[0].x[dim])
        return 0;

    if (tmp[0].add != -1) {
        long long ans = 0;
        long long s[3];
        s[0] = 0;
        for (int j = 0; j < B - 1; j++) {
            long long i = j % P;
            if (j != 0 && i == 0)
                ReadP(node + j, tree, tmp);
            if (tmp[i].x[dim] > end)
                break;
            if (start <= tmp[i].x[dim] && end >= tmp[i + 1].x[dim])
                ans = F(ans, tmp[i].w);
            else if (start <= tmp[i + 1].x[dim] && end >= tmp[i].x[dim])
                s[++s[0]] = tmp[i].add;
        }
        for (int i = 1; i <= s[0]; i++)
            ans = F(ans, SegQuery(tree, s[i], start, end, dim));
        return ans;
    } else {
        long long ans = 0;
        ReadB(node + P, tree, tmp + P);
        for (int i = 0; i < B; i++) {
            if (tmp[i].x[dim] > end)
                break;
            if (tmp[i].x[dim] >= start && tmp[i].x[dim] <= end)
                ans = F(ans, tmp[i].w);
        }
        return ans;
    }
}

long long PartQuery(int tree, long long node, double start_x, double end_x, double start_y, double end_y){
    FILE *fp = LSMtree[tree], * fpT = Segtree[tree];

    Read(node, fp, tmp2);
    long long ed = tmp2[0].add;

    double y1, y2;
    long long Part1[100];
    Part1[0] = 0;
    long long ans = 0;
    for (long long i = node; i < ed; i += B){
        if (i != node)
            Read(i, fp, tmp2);
        for (int j = 0; j < B; j++){
            if (tmp2[j].w == -1)
                break;
            y2 = tmp2[j].x[1];
            if (i != node || j != 0){
                if (y1 > end_y)
                    break;
                if (start_y <= y1 && end_y >= y2)
                    ans = F(ans, SegQuery(fpT, tmp2[j].add, start_x, end_x, 0));
                else if (start_y <= y2 && end_y >= y1) {
                    Part1[++Part1[0]] = tmp2[j].w;
                }
            }
            y1 = y2;
        }
    }

    for (int p = 1; p <= Part1[0]; p++){
        Read(Part1[p], fp, tmp2);
        ed = tmp2[0].add;
        double x1, x2;
        long long start, end;
        for (long long i = Part1[p]; i < ed; i += B){
            if (i != Part1[p])
                Read(i, fp, tmp2);
            for (int j = 0; j < B; j++){
                if (tmp2[j].w == -1)
                    break;
                x2 = tmp2[j].x[0];
                end = tmp2[j].w;
                if (i != Part1[p] || j != 0){
                    if (x1 > end_x)
                        break;
                    if (start_x <= x1 && end_x >= x2) {
                        ans = F(ans, SegQuery(fpT, tmp2[j].add, start_y, end_y, 1));
                    }
                    else if (start_x <= x2 && end_x >= x1)
                        for (long long k = start; k < end; k += B){
                            Read(k, fp, tmp);
                            for (int q = 0; q < B; q++)
                                if (tmp[q].x[0] >= start_x && tmp[q].x[0] <= end_x && tmp[q].x[1] >= start_y && tmp[q].x[1] <= end_y)
                                    ans = F(ans, tmp[q].w);
                        }
                }
                x1 = x2;
                start = end;
            }
        }
    }
    return ans;
}

long long Search(int tree, long long node, data start, data end, int dim){
    FILE *fp = LSMtree[tree];

    Read(node, fp, itmp[dim]);
    long long ed = itmp[dim][0].add;
    long long ans = 0;

    if (dim <= 1) {
        long long w1, w2;
        double x1, x2;
        for (long long i = node; i < ed; i += B) {
            if (i != node)
                Read(i, fp, itmp[dim]);
            for (int j = 0; j < B; j++) {
                if (itmp[dim][j].w == -1)
                    break;
                x2 = itmp[dim][j].x[dim];
                w2 = itmp[dim][j].w;
                if (i != node || j != 0) {
                    if (start.x[dim] <= x2 && end.x[dim] >= x1) {
                        if (dim == 1)
                            ans = F(ans, Search(tree, itmp[dim][j].w, start, end, dim - 1));
                        else {
                            for (long long k = w1; k < w2; k += B) {
                                Read(k, fp, tmp);
                                for (int q = 0; q < B; q++)
                                    if (Check(tmp[q], start, end))
                                        ans = F(ans, tmp[q].w);
                            }
                        }
                    }
                }
                x1 = x2;
                w1 = w2;
            }
        }
    }
    else{
        for (long long i = node; i < ed; i += B) {
            if (i != node)
                Read(i, fp, itmp[dim]);
            for (int j = 0; j < B; j += 2) {
                if (itmp[dim][j].w == -1)
                    break;
                if (start.x[dim] <= itmp[dim][j+1].x[dim] && end.x[dim] >= itmp[dim][j].x[dim]) {
                        ans = F(ans, Search(tree, itmp[dim][j + 1].add, start, end, dim - 1));
                }
            }
        }
    }
    return ans;
}

long long WorqQuery(int tree, long long node, data start, data end, int dim){
    if (dim == 1)
        return PartQuery(tree, node, start.x[0], end.x[0], start.x[1], end.x[1]);

    FILE *fp = LSMtree[tree];

    Read(node, fp, itmp[dim]);
    long long ed = itmp[dim][0].add;

    long long ans = 0;
    for (long long i = node; i < ed; i += B){
        if (i != node)
            Read(i, fp, itmp[dim]);
        for (int j = 0; j < B; j += 2){
            if (itmp[dim][j].w == -1)
                break;
            if (start.x[dim] <= itmp[dim][j].x[dim] && end.x[dim] >= itmp[dim][j + 1].x[dim])
                ans = F(ans, WorqQuery(tree, itmp[dim][j + 1].add, start, end, dim - 1));
            else if (start.x[dim] <= itmp[dim][j + 1].x[dim] && end.x[dim] >= itmp[dim][j].x[dim])
                ans = F(ans, Search(tree, itmp[dim][j + 1].add, start, end, dim - 1));
        }
    }
    return ans;
}

int main(){
    FILE *fp;
    fp = fopen("DB.txt","r+");
    tmpFile = fopen("temp.txt", "w+");
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
                anss = F(anss, WorqQuery(j, node, x, y, D - 1));
            node *= 2;
        }
        end   = system_clock::now();
        duration = duration_cast<microseconds>(end - start);
        result << fixed << setprecision(10) << double(duration.count()) * microseconds::period::num / microseconds::period::den << endl;
    }

    fclose(fp);
    fclose(tmpFile);
    for (int i = 0; i < LSMCnt; i++) {
        fclose(LSMtree[i]);
        fclose(Segtree[i]);
    }
    for (int i = 0; i < LSMCnt - 1; i++) {
        fclose(LSMtree[i + 100]);
        fclose(Segtree[i + 100]);
    }
    result.close();

    return 0;
}