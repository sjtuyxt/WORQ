//Worq
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
const int P = 200;
const int M = 50;
const long long n = 60000;

struct data {
    double x[D];
    long long w, add;
};

unsigned char buffer[B * (D + 2)][8];
data tmp[B], tmp2[B], tmp3[B], tmp4[B], tmp5[B], tmp6[B];
data itmp[M][B], tmp1[M][B];
char str[100];
int MCnt = 0, LSMCnt = 0;
long long LSMNum[100], LSMEnd[100], DataEnd[100], DataStart[100];
FILE *LSMtree[100], *Segtree[100], *tmpFile;

bool l0(const data &x, const data &y) { return x.x[0] < y.x[0]; }

bool l1(const data &x, const data &y) { return x.x[1] < y.x[1]; }

bool l2(const data &x, const data &y) { return x.x[2] < y.x[2]; }

bool l3(const data &x, const data &y) { return x.x[3] < y.x[3]; }

bool l4(const data &x, const data &y) { return x.x[4] < y.x[4]; }

bool (*cmp[5])(const data &, const data &) = {l0, l1, l2, l3, l4};

bool g0(const data &x, const data &y) { return x.x[0] > y.x[0]; }

bool g1(const data &x, const data &y) { return x.x[1] > y.x[1]; }

bool g2(const data &x, const data &y) { return x.x[2] > y.x[2]; }

bool g3(const data &x, const data &y) { return x.x[3] > y.x[3]; }

bool g4(const data &x, const data &y) { return x.x[4] > y.x[4]; }

bool (*gcmp[5])(const data &, const data &) = {g0, g1, g2, g3, g4};

bool Check(data x, data start, data end) {
    for (int i = 0; i < D; i++)
        if (x.x[i] < start.x[i] || x.x[i] > end.x[i])
            return false;
    return true;
}

long long F(long long x, long long y) {
    return x ^ y;
}

long long DoubleToLong(double x) {
    long long y;
    memcpy(&y, &x, 8);
    return y;
}

double LongToDouble(long long x) {
    double y;
    memcpy(&y, &x, 8);
    return y;
}

void LongToChar(long long x, unsigned char *ch) {

    for (int i = 7; i >= 0; i--) {
        ch[i] = (int) (x & 255);
        x >>= 8;
    }
}

long long CharToLong(const unsigned char *ch) {
    long long x = 0;
    for (int i = 0; i < 8; i++) {
        x = (x << 8) + ch[i];
    }
    return x;
}

long long IOcnt = 0, SegIO = 0;

void Read(long long start, FILE *fp, data *input) {
    fseeko(fp, start * (D + 2) * 8, SEEK_SET);
    fread(buffer, sizeof(char), (D + 2) * B * 8, fp);
    for (int i = 0; i < B; i++) {
        for (int j = 0; j < D; j++) {
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

void Write(long long start, FILE *fp, data *output) {
    for (int i = 0; i < B; i++) {
        for (int j = 0; j < D; j++) {
            LongToChar(DoubleToLong(output[i].x[j]), buffer[i * (D + 2) + j]);
        }
        LongToChar(output[i].w, buffer[i * (D + 2) + D]);
        LongToChar(output[i].add, buffer[i * (D + 2) + D + 1]);
    }
    fseeko(fp, start * (D + 2) * 8, SEEK_SET);
    fwrite(buffer, sizeof(char), (D + 2) * B * 8, fp);
    IOcnt++;
}


void RandomDB(long long n) {
    int x = time(NULL);
//    x = 1617640380;
    srand(x);
    cout << x << endl;
    FILE *fp;
    fp = fopen("DB.txt", "w+");

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

void SaveDB() {
    ofstream out;
    out.open("DBvar.txt");

    out << LSMCnt << endl;
    for (int i = 0; i < LSMCnt; i++)
        out << LSMNum[i] << " " << LSMEnd[i] << " " << DataEnd[i] << endl;

    out.close();
}

void LoadDB() {
    ifstream in;
    in.open("DBvar.txt");

    in >> LSMCnt;
    for (int i = 0; i < LSMCnt; i++)
        in >> LSMNum[i] >> LSMEnd[i] >> DataEnd[i];

    in.close();

    for (int i = 0; i < LSMCnt; i++) {
        sprintf(str, "%d.txt", i);
        LSMtree[i] = fopen(str, "r+");
        sprintf(str, "%dT.txt", i);
        Segtree[i] = fopen(str, "r+");
    }
}

void ExSort(long long start, long long end, int dim, FILE *fp) {
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

long long SegMerge(int tree, long long start, long long end, int dim) {
    FILE *fp = LSMtree[tree], *fpT = Segtree[tree];

//    if (dim == 0)
//        cout << "seg: " << start << " " << end << endl;

    ExSort(start, end, dim, fp);

    for (long long i = start; i < end; i += B) {
        Read(i, fp, tmp);
        Write(i - start + LSMEnd[tree], fpT, tmp);
    }

    long long i = LSMEnd[tree];
    LSMEnd[tree] += end - start;
    while (i + B < LSMEnd[tree]) {
        int cnt = 1;
        Read(i, fpT, tmp);
        tmp2[1] = tmp[B - 1];
        tmp2[0] = tmp[0];
        tmp2[0].add = i;
        for (int j = 1; j < B; j++)
            tmp2[0].w = F(tmp2[0].w, tmp[j].w);
        while (i + cnt * B < LSMEnd[tree] && cnt < B - 1) {
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
        for (int j = cnt; j < B; j++)
            tmp2[j].w = 0;
        if (cnt != B - 1) {
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

long long PartMerge(int tree, long long start, long long end, int dim) {
    FILE *fp = LSMtree[tree];

    if (dim == 1)
        ExSort(start, end, dim, fp);

    long long head = DataEnd[tree];
    int cnt = 0;
    long long length = pow((end - start) / B, 2 / 3.0);

    if (dim == 0)
        length = pow((end - start) / B, 1 / 2.0);

//    cout << dim << " " << (end - start) / B << " " << (end - start) / B / length << endl;
    for (long long i = start; i < end; i += B * length) {
        Read(i, fp, tmp2);
        tmp[cnt] = tmp2[0];
        tmp[cnt].w = i;
        cnt++;
        if (cnt == B) {
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
    for (long long i = head; i < ed; i += B) {
        Read(i, fp, tmp1[dim]);

        for (int j = 0; j < B; j++) {
            if (tmp1[dim][j].w == -1)
                break;
            y = tmp1[dim][j].w;
            if (i != head || j != 0) {
                tmp1[dim][j].add = SegMerge(tree, x, y, dim ^ 1);
                if (dim == 1)
                    tmp1[dim][j].w = PartMerge(tree, x, y, 0);
            } else
                tmp1[dim][j].add = ed;
            x = y;
        }
        Write(i, fp, tmp1[dim]);
    }

    return head;
}

long long WorqMerge(int tree, long long start, long long end, int dim) {
    if (dim == 1)
        return PartMerge(tree, start, end, dim);

    FILE *fp = LSMtree[tree];

    ExSort(start, end, dim, fp);

    long long head = DataEnd[tree];
    int cnt = 0;
    long long length = pow((end - start) / B, (2.0 * dim - 3) / (2.0 * dim - 1)) * 3;

    for (long long i = start; i < end; i += B * length) {
        Read(i, fp, tmp2);
        tmp[cnt] = tmp2[0];
        tmp[cnt].w = i;
        cnt++;
        if (cnt == B) {
            cnt = 0;
            Write(DataEnd[tree], fp, tmp);
            DataEnd[tree] += B;
        }
    }
    Read(end - B, fp, tmp2);
    tmp[cnt] = tmp2[B - 1];
    tmp[cnt].w = end;
    cnt++;
    long long ed = DataEnd[tree] + cnt;
    if (DataEnd[tree] == head) {
        tmp[0].add = ed;
        Write(DataEnd[tree], fp, tmp);
        DataEnd[tree] += B;
    } else {
        Write(DataEnd[tree], fp, tmp);
        Read(head, fp, tmp);
        tmp[0].add = ed;
        Write(head, fp, tmp);
        DataEnd[tree] += B;
    }

//    for (int i = 0; i < cnt; i++)
//        cout << tmp[i].w << endl;

    long long st = head;
    long long nextst = head;
    long long nexted = ed;
    long long x = 0;
    long long y = 0;
    long long beg1 = 0;
    long long beg2 = 0;

    for (int i = dim - 1; i >= 0; i--){
        cnt = 0;
        st = nextst;
        nextst = DataEnd[tree];
        ed = nexted;
//        cout << "i " << i << " " << st << " " << ed  << endl;
        for (long long j = st; j < ed; j += B){
            Read(j, fp, tmp3);
//            for (int k = 0; k < B && j + k < ed; k++)
//                cout << tmp3[k].w << " ";
//            cout << endl;
            for (int k = 0; k < B && j + k < ed; k++){
                x = y;
                y = tmp3[k].w;
                if (tmp3[k].add == -1){
//                    cout << "xy " << x << " " << y << endl;
                    if (i == 0){
                        tmp3[k].add = SegMerge(tree, x, y, 0);
                        tmp3[k].w = DataEnd[tree] + cnt;
                    }
                    else {
                        tmp3[k].add = DataEnd[tree] + cnt;
                        ExSort(x, y, i, fp);
                    }
                    if (i > 1) {
                        length = pow((y - x) / B, (2.0 * i - 3) / (2.0 * i - 1)) * 3;
                        if (length > (y - x) / B)
                            length = (y - x) / B;
                    }
                    else if (i == 1)
                        length = pow((y - x) / B, 2 / 3.0);
                    else
                        length = pow((y - x) / B, 1 / 2.0);
                    long long p = x;
                    beg1 = DataEnd[tree];
                    beg2 = cnt;
                    for (; p < y; p += B * length) {
                        Read(p, fp, tmp4);
                        tmp5[cnt] = tmp4[0];
                        tmp5[cnt].w = p;
                        if (i == 0 && p != x)
                            tmp5[cnt].add = SegMerge(tree, p - B * length, p, 1);
                        cnt++;
                        if (cnt == B) {
                            cnt = 0;
                            Write(DataEnd[tree], fp, tmp5);
                            DataEnd[tree] += B;
                        }
                    }
                    Read(y - B, fp, tmp4);
                    tmp5[cnt] = tmp4[B - 1];
                    tmp5[cnt].w = y;
                    if (i == 0)
                        tmp5[cnt].add = SegMerge(tree, p - B * length, y, 1);
                    cnt++;
                    if (beg1 == DataEnd[tree])
                        tmp5[beg2].add = DataEnd[tree] + cnt;
                    else {
                        Read(beg1, fp, tmp6);
                        tmp6[beg2].add = DataEnd[tree] + cnt;
                        Write(beg1, fp, tmp6);
                    }
                    if (cnt == B) {
                        cnt = 0;
                        Write(DataEnd[tree], fp, tmp5);
                        DataEnd[tree] += B;
                        nexted = DataEnd[tree];
                    }
                }
            }
            Write(j, fp, tmp3);
        }
        if (cnt) {
            nexted = DataEnd[tree] + cnt;
            Write(DataEnd[tree], fp, tmp5);
            DataEnd[tree] += B;
        }
    }

    return head;
}

void Clear(int tree) {

    if (LSMNum[tree]) {
        fclose(LSMtree[tree]);
        fclose(Segtree[tree]);
    }
    sprintf(str, "%d.txt", tree);
    LSMtree[tree] = fopen(str, "w+");
    sprintf(str, "%dT.txt", tree);
    Segtree[tree] = fopen(str, "w+");
}

void Insert(data *input) {
    memcpy(itmp[MCnt], input, (D + 2) * B * 8);
    MCnt++;

    if (MCnt == M) {
        MCnt = 0;
        int dst = 0;
        while (LSMNum[dst])
            dst++;
        if (dst + 1 > LSMCnt)
            LSMCnt = dst + 1;
        Clear(dst);
        LSMNum[dst] = 1;
        FILE *fp = LSMtree[dst];
        long long length = 1;
        for (int j = 0; j < M; j++)
            Write(j * B, fp, itmp[j]);
        for (int i = 0; i < dst; i++) {
            FILE *tfp = LSMtree[i];
            for (long long j = 0; j < length * M; j++) {
                Read(j * B, tfp, tmp);
                Write(length * M * B + j * B, fp, tmp);
            }
            Clear(i);
            length *= 2;
            LSMNum[i] = 0;
            LSMEnd[i] = 0;
            DataEnd[i] = 0;
        }
        DataEnd[dst] = length * M * B;
//        cout << endl << "worqM: " << dst << endl;
        WorqMerge(dst, 0, length * M * B, D - 1);
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

long long PartQuery(int tree, long long node, double start_x, double end_x, double start_y, double end_y) {
    FILE *fp = LSMtree[tree], *fpT = Segtree[tree];

    if (node < DataStart[1] || node >= DataStart[1] + B) {
        Read(node, fp, itmp[1]);
        DataStart[1] = node;
    }
    long long ed = itmp[1][node - DataStart[1]].add;

    double y1 = 0;
    double y2 = 0;

    long long ans = 0;
    bool loop;
    for (long long i = node; i < ed; i++) {
        if (i >= DataStart[1] + B){
            Read(i, fp, itmp[1]);
            DataStart[1] = i;
        }
        long long j = i - DataStart[1];
        y2 = itmp[1][j].x[1];
        if (i != node) {
            if (y1 > end_y)
                break;
            if (start_y <= y1 && end_y >= y2) {
                ans = F(ans, SegQuery(fpT, itmp[1][j].add, start_x, end_x, 0));
            }
            else if (start_y <= y2 && end_y >= y1) {

                long long x_st = itmp[1][j].w;

                if (x_st < DataStart[0] || x_st >= DataStart[0] + B) {
                    Read(x_st, fp, itmp[0]);
                    DataStart[0] = x_st;
                }
                long long ed = itmp[0][x_st - DataStart[0]].add;

                double x1, x2;
                long long start, end;
                for (long long k = x_st; k < ed; k++) {
                    if (k >= DataStart[0] + B){
                        Read(k, fp, itmp[0]);
                        DataStart[0] = k;
                    }
                    long long j = k - DataStart[0];
                    x2 = itmp[0][j].x[0];
                    end = itmp[0][j].w;
                    if (k != x_st) {
                        if (x1 > end_x)
                            break;
                        if (start_x <= x1 && end_x >= x2)
                            ans = F(ans, SegQuery(fpT, itmp[0][j].add, start_y, end_y, 1));
                        else if (start_x <= x2 && end_x >= x1) {
                            loop = true;
                            for (long long k = start; k < end && loop; k += B) {
                                Read(k, fp, tmp);
                                for (int q = 0; q < B; q++) {
                                    if (tmp[q].x[1] > end_y){
                                        loop = false;
                                        break;
                                    }
                                    if (tmp[q].x[0] >= start_x && tmp[q].x[0] <= end_x && tmp[q].x[1] >= start_y)
                                        ans = F(ans, tmp[q].w);
                                }
                            }
                        }
                    }
                    x1 = x2;
                    start = end;
                }

            }
        }
        y1 = y2;
    }

//    cout << fixed << setprecision(1) << "part: " << node << " " << 1 << " " << ans << endl;

    return ans;
}

//long long SegSearch(FILE *tree, long long node, data start, data end, int dim) {
//    Read(node, tree, tmp);
//
//    if (end.x[dim] < tmp[0].x[dim])
//        return 0;
//
//    if (tmp[0].add != -1) {
//        long long ans = 0;
//        long long s[B + 1];
//        s[0] = 0;
//        for (int i = 0; i < B - 1; i++) {
//            if (tmp[i].x[dim] > end.x[dim])
//                break;
//            if (start.x[dim] <= tmp[i + 1].x[dim] && end.x[dim] >= tmp[i].x[dim])
//                s[++s[0]] = tmp[i].add;
//        }
//        for (int i = 1; i <= s[0]; i++)
//            ans = F(ans, SegSearch(tree, s[i], start, end, dim));
//        return ans;
//    } else {
//        long long ans = 0;
//        for (int i = 0; i < B; i++) {
//            if (tmp[i].x[dim] > end.x[dim])
//                break;
//            if (Check(tmp[i], start, end))
//                ans = F(ans, tmp[i].w);
//        }
//        return ans;
//    }
//}

long long Search(int tree, long long node, data start, data end, int dim) {
    FILE *fp = LSMtree[tree];
    if (node < DataStart[dim] || node >= DataStart[dim] + B) {
        Read(node, fp, itmp[dim]);
        DataStart[dim] = node;
    }
    long long ed = itmp[dim][node - DataStart[dim]].add;
    long long ans = 0;
    long long w1, w2;


    double x1, x2;
    bool loop;
    for (long long i = node; i < ed; i++) {
        if (i >= DataStart[dim] + B){
            Read(i, fp, itmp[dim]);
            DataStart[dim] = i;
        }
        long long j = i - DataStart[dim];
        x2 = itmp[dim][j].x[dim];
        w2 = itmp[dim][j].w;
        if (i != node) {
            if (x1 > end.x[dim])
                break;
            if (start.x[dim] <= x2 && end.x[dim] >= x1) {
                if (dim > 1)
                    ans = F(ans, Search(tree, itmp[dim][j].add, start, end, dim - 1));
                else if (dim == 1)
                    ans = F(ans, Search(tree, itmp[dim][j].w, start, end, dim - 1));
                else {
                    loop = true;
                    for (long long k = w1; k < w2 && loop; k += B){
                        Read(k, fp, tmp);
                        for (int q = 0; q < B; q++) {
                            if (tmp[q].x[1] > end.x[1]) {
                                loop = false;
                                break;
                            }
                            if (Check(tmp[q], start, end)) {
                                ans = F(ans, tmp[q].w);
                            }
                        }
                    }
                }
            }
        }
        x1 = x2;
        w1 = w2;
    }

    return ans;
}

long long WorqQuery(int tree, long long node, data start, data end, int dim) {
    if (dim == 1)
        return PartQuery(tree, node, start.x[0], end.x[0], start.x[1], end.x[1]);

    FILE *fp = LSMtree[tree];

    if (node < DataStart[dim] || node >= DataStart[dim] + B) {
        Read(node, fp, itmp[dim]);
        DataStart[dim] = node;
    }
    long long ed = itmp[dim][node - DataStart[dim]].add;

    double x1, x2;
    long long ans = 0;
    for (long long i = node; i < ed; i++) {
        if (i >= DataStart[dim] + B){
            Read(i, fp, itmp[dim]);
            DataStart[dim] = i;
        }
        long long j = i - DataStart[dim];
        x2 = itmp[dim][j].x[dim];
        if (i != node) {
            if (x1 > end.x[dim])
                break;
            if (start.x[dim] <= x1 && end.x[dim] >= x2)
                ans = F(ans, WorqQuery(tree, itmp[dim][j].add, start, end, dim - 1));
            else if (start.x[dim] <= x2 && end.x[dim] >= x1)
                ans = F(ans, Search(tree, itmp[dim][j].add, start, end, dim - 1));
        }
        x1 = x2;
    }

    return ans;
}

int main() {
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
        IOcnt = 0;
        for (int j = 0; j < LSMCnt; j++) {
            if (LSMNum[j]) {
                anss = F(anss, WorqQuery(j, node, x, y, D - 1));
            }
            node *= 2;
        }
        end = system_clock::now();
        duration = duration_cast<microseconds>(end - start);
        result << fixed << setprecision(10)
             << double(duration.count()) * microseconds::period::num / microseconds::period::den << "\t" << endl;
    }

    fclose(fp);
    fclose(tmpFile);
    for (int i = 0; i < LSMCnt; i++) {
        fclose(LSMtree[i]);
        fclose(Segtree[i]);
    }
    result.close();

    return 0;
}