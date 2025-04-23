#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <chrono>
#include <cmath>
#include <numeric>
#include <climits>
#include <functional>  
#include<assert.h>
#include <stack>

using namespace std;

//==== Data structures ====  
struct HardBlock { string name; int width, height; };
struct Pad { string name; int x, y; };
struct Net { string name; vector<string> pins; };

struct Data {
    vector<HardBlock>                blocks;
    unordered_map<string, int>        blockIndex;
    vector<Pad>                      pads;
    unordered_map<string, pair<int, int>> padPos;
    vector<Net>                      nets;
    long long                        totalBlockArea = 0;
};

struct Shape {
    int w, h;
    vector<pair<int, int>> offsets;  // offsets[i] = (x,y) of block i
};

//==== Read input .txt ====  
bool readInput(const string& filename, Data& in) {
    ifstream fin(filename);
    if (!fin) return false;
    string token;
    int numHard;
    fin >> token >> numHard;  // NumHardBlocks
    in.blocks.reserve(numHard);
    for (int i = 0; i < numHard; i++) {
        string kind, name; int w, h;
        fin >> kind >> name >> w >> h;
        in.blockIndex[name] = in.blocks.size();
        in.blocks.push_back({ name, w, h });
        in.totalBlockArea += 1LL * w * h;
    }
    int numPads;
    fin >> token >> numPads;  // NumPads
    in.pads.reserve(numPads);
    for (int i = 0; i < numPads; i++) {
        string kind, name; int x, y;
        fin >> kind >> name >> x >> y;
        in.pads.push_back({ name, x, y });
        in.padPos[name] = { x, y };
    }
    int numNets;
    fin >> token >> numNets;  // NumNets
    in.nets.reserve(numNets);
    for (int i = 0; i < numNets; i++) {
        string kind, netName; int pinCnt;
        fin >> kind >> netName >> pinCnt;
        Net net; net.name = netName;
        for (int j = 0; j < pinCnt; j++) {
            fin >> kind; // Pin
            string pinName;
            fin >> pinName;
            net.pins.push_back(pinName);
        }
        in.nets.push_back(move(net));
    }
    return true;
}

//==== Compute HPWL ====  
double computeHPWL(const Data& in, const vector<pair<int, int>>& pos)
{
    double total = 0;
    for (auto& net : in.nets) {
        int xmin = INT_MAX, xmax = INT_MIN;
        int ymin = INT_MAX, ymax = INT_MIN;
        for (auto& p : net.pins) {
            int cx, cy;
            auto itb = in.blockIndex.find(p);
            if (itb != in.blockIndex.end()) {
                // Block pin
                int i = itb->second;
                int x = pos[i].first, y = pos[i].second;
                int w = in.blocks[i].width, h = in.blocks[i].height;
                // 不需要再進行 swap，因為 rot 已經在其他地方處理過
                cx = x + w / 2;
                cy = y + h / 2;
            }
            else {
                // Pad pin
                auto pp = in.padPos.at(p);
                cx = pp.first;
                cy = pp.second;
            }
            xmin = min(xmin, cx);
            xmax = max(xmax, cx);
            ymin = min(ymin, cy);
            ymax = max(ymax, cy);
        }
        total += (xmax - xmin) + (ymax - ymin);
    }
    return total;
}

//==== Write output ====  
void writeOutput(const string& filename,
    long long wirelength,
    const Data& in,
    const vector<pair<int, int>>& pos,
    const vector<int>& rot)
{
    ofstream fout(filename);
    fout << "Wirelength " << wirelength << "\n\n";
    fout << "NumHardBlocks " << in.blocks.size() << "\n";
    for (size_t i = 0; i < in.blocks.size(); ++i) {
        fout
            << in.blocks[i].name << " "
            << pos[i].first << " "
            << pos[i].second << " "
            << rot[i] << "\n";
    }
}

//==== Preprocess: rotate so height >= width ====  
void normalizeTall(Data& in, vector<int>& rot) {
    int n = in.blocks.size();
    rot.resize(n);
    for (int i = 0; i < n; ++i) {
        auto& b = in.blocks[i];
        if (b.width > b.height) {
            swap(b.width, b.height);
            rot[i] = 1;
        }
        else {
            rot[i] = 0;
        }
    }
}

vector<int> initialLB2(Data& in, int W, int H, vector<pair<int, int>>& pos, vector<int>& rot) {
    int n = in.blocks.size();
    pos.assign(n, { 0, 0 });
    vector<bool> used(n, false);
    int placedCount = 0;

    // Prepare sorted lists by height and by area
    vector<int> heightOrder(n), areaOrder(n);
    iota(heightOrder.begin(), heightOrder.end(), 0);
    iota(areaOrder.begin(), areaOrder.end(), 0);
    sort(heightOrder.begin(), heightOrder.end(), [&](int a, int b) {
        return in.blocks[a].height > in.blocks[b].height;
        });
    sort(areaOrder.begin(), areaOrder.end(), [&](int a, int b) {
        long long A = 1LL * in.blocks[a].width * in.blocks[a].height;
        long long B = 1LL * in.blocks[b].width * in.blocks[b].height;
        return A > B;
        });

    vector<int> polishExpr;
    int heightPtr = 0;
    int y = 0;

    while (placedCount < n && y < H) {
        int x = 0;
        // Step 2: pick next tallest unused block
        while (heightPtr < n && used[heightOrder[heightPtr]]) ++heightPtr;
        if (heightPtr >= n) break;
        int b1 = heightOrder[heightPtr];
        int h1 = in.blocks[b1].height;
        used[b1] = true;
        pos[b1] = { 0, y };
        placedCount++;
        x += in.blocks[b1].width;

        // Add block index to Polish expression
        polishExpr.push_back(b1);

        if (y + h1 > H) {
            cout << "Height overflow, break!" << endl;
            break;
        }

        // Step 3: place b2 and try vertical stack
        while (true) {
            while (heightPtr < n && used[heightOrder[heightPtr]]) ++heightPtr;
            if (heightPtr >= n) break;
            int b2 = heightOrder[heightPtr];
            int w2 = in.blocks[b2].width;
            int h2 = in.blocks[b2].height;
            if (x + w2 > W) break;

            // find best vertical partner b3
            int best = -1;
            for (int c : areaOrder) {
                if (used[c] || c == b2) continue;
                int w3 = in.blocks[c].width;
                int h3 = in.blocks[c].height;
                if (w3 <= w2 && h3 <= h1 - h2) {
                    best = c;
                    break;
                }
                else if (h3 <= w2 && w3 <= h1 - h2) {
                    best = c;
                    swap(in.blocks[c].width, in.blocks[c].height);
                    rot[c] = 1; // rotate
                    break;
                }
            }
            if (best >= 0) {
                // place vertical stack
                pos[b2] = { x, y };
                pos[best] = { x, y + h2 };
                used[b2] = used[best] = true;
                placedCount += 2;

                // Add block indices and vertical operator to Polish expression
                polishExpr.push_back(b2);
                polishExpr.push_back(best);
                polishExpr.push_back(-2);
                //cout << "Placed vertical stack: " << in.blocks[b2].name << " and " << in.blocks[best].name << endl;
                x += w2;
            }
            else {
                // place b2 alone
                pos[b2] = { x, y };
                used[b2] = true;
                placedCount++;

                // Add block index to Polish expression
                polishExpr.push_back(b2);

                x += w2;
            }
            polishExpr.push_back(-1);
        }

        // Step 4: fill remaining width in row
        int remW = W - x;
        if (remW > 0) {
            int best = -1;
            for (int c : areaOrder) {
                if (used[c]) continue;
                int w3 = in.blocks[c].width;
                int h3 = in.blocks[c].height;
                if (w3 <= remW && h3 <= h1) {
                    best = c;
                    break;
                }
                else if (h3 <= remW && w3 <= h1) {
                    best = c;
                    swap(in.blocks[c].width, in.blocks[c].height);
                    rot[c] = 1; // rotate
                    break;
                }
            }
            if (best >= 0) {
                pos[best] = { x, y };
                used[best] = true;
                placedCount++;

                // Add block index to Polish expression
                polishExpr.push_back(best);
                polishExpr.push_back(-1);
            }
        }

        // Add horizontal operator to Polish expression
        if (y > 0) {
            polishExpr.push_back(-2); // Horizontal operator
        }

        // Step 5: move to next row
        y += h1;
    }

    if (placedCount != n || polishExpr.size() != 2 * n - 1)
        cout << "Not all blocks placed, fuck!" << endl;

    return polishExpr;
}


bool evalExpr(const vector<int>& expr, const Data& in, int W, int H, vector<pair<int, int>>& pos) {
    int n = in.blocks.size();
    struct Frame {
        int w, h; // Frame dimensions
        int x, y; // Bottom-left corner of the frame
        vector<pair<int, int>> offsets; // Offsets for blocks
    };
    stack<Frame> stk;
    pos.assign(n, { -1, -1 }); // Initialize positions with invalid values

    for (int tok : expr) {
        if (tok >= 0) {
            // Block token
            Frame f;
            f.w = in.blocks[tok].width;
            f.h = in.blocks[tok].height;
            f.x = 0; // Default x-coordinate
            f.y = 0; // Default y-coordinate
            f.offsets.push_back({ 0, 0 }); // Block's own position
            stk.push(move(f));
        }
        else {
            // Operator token (H or V)
            if (stk.size() < 2) return false; // Invalid expression
            Frame R = move(stk.top()); stk.pop();
            Frame L = move(stk.top()); stk.pop();
            Frame C;

            if (tok == -2) {
                // Horizontal stacking (H): L below, R above
                C.w = max(L.w, R.w);
                C.h = L.h + R.h;
                C.x = L.x;
                C.y = L.y;

                // Adjust offsets for L and R
                for (auto& offset : L.offsets) {
                    C.offsets.push_back(offset); // L's offsets remain unchanged
                }
                for (auto& offset : R.offsets) {
                    C.offsets.push_back({ offset.first, offset.second + L.h }); // R is above L
                }
            }
            else if (tok == -1) {
                // Vertical stacking (V): L to the left, R to the right
                C.w = L.w + R.w;
                C.h = max(L.h, R.h);
                C.x = L.x;
                C.y = L.y;

                // Adjust offsets for L and R
                for (auto& offset : L.offsets) {
                    C.offsets.push_back(offset); // L's offsets remain unchanged
                }
                for (auto& offset : R.offsets) {
                    C.offsets.push_back({ offset.first + L.w, offset.second }); // R is to the right of L
                }
            }
            else {
                return false; // Invalid operator
            }
            stk.push(move(C));
        }
    }

    if (stk.size() != 1) return false; // Invalid expression

    Frame result = move(stk.top());
    if (result.w > W || result.h > H) return false; // Exceeds bounding box

    // Assign positions to blocks
    int offsetIdx = 0;
    for (size_t i = 0; i < expr.size(); ++i) {
        if (expr[i] >= 0) { // 如果是 block token
            int blockIdx = expr[i]; // 使用 Polish 表達式中的 block 索引
            pos[blockIdx] = { result.offsets[offsetIdx].first, result.offsets[offsetIdx].second };
            offsetIdx++;
        }
    }

    return true;
}


vector<int> simulatedAnneal(const Data& in,  vector<int> expr, const vector<int>& rot, int W, int H, vector<pair<int, int>> pos){
    int nBlocks = in.blocks.size();
    vector<int> bestExpr = expr;

    double bestCost = computeHPWL(in, pos);

    mt19937_64 rng(0);
    uniform_real_distribution<double> ur(0.0, 1.0);

    double T = 1000.0;
    const double cooling = 0.85;

    while (T > 0.001) {
        for (int it = 0; it < 5 * nBlocks; ++it) {
            vector<int> cand = expr;
            double r = ur(rng);

            if (r < 0.33) {
                // —— M1: ?机挑一?相?的 operand 互?
                vector<int> idx;
                for (int i = 0; i + 1 < (int)cand.size(); ++i)
                    if (cand[i] >= 0 && cand[i + 1] >= 0)
                        idx.push_back(i);
                if (idx.empty()) continue;
                int i = idx[rng() % idx.size()];
                swap(cand[i], cand[i + 1]);
            }
            else if (r < 0.66) {
                // —— M2: ?机挑一? operator ?反?
                vector<int> idx;
                for (int i = 0; i < (int)cand.size(); ++i)
                    if (cand[i] < 0)
                        idx.push_back(i);
                if (idx.empty()) continue;
                int i = idx[rng() % idx.size()];
                int j = i;
                while (j + 1 < (int)cand.size() && cand[j + 1] < 0 && cand[j + 1] != cand[j])
                    ++j;
                for (int k = i; k <= j; ++k)
                    cand[k] = (cand[k] == -1 ? -2 : -1);
            }
            else {
                // —— M3: ?机挑一? operator 与后面 operand 互?
                vector<int> idx;
                for (int i = 0; i + 1 < (int)cand.size(); ++i)
                    if (cand[i] < 0 && cand[i + 1] >= 0)
                        idx.push_back(i);
                if (idx.empty()) continue;
                int i = idx[rng() % idx.size()];
                // 交? operator 和后面的???
                swap(cand[i], cand[i + 1]);
            }

            // 用 evalExpr ??并?生新的 pos，?便?掉越界或?构不合法
            vector<pair<int, int>> posC;
            if (!evalExpr(cand, in, W, H, posC))
                continue;

            double costC = computeHPWL(in, posC);
            vector<pair<int, int>> posCur;
            evalExpr(expr, in, W, H, posCur);
            double costCur = computeHPWL(in, posCur);
            double dC = costC - costCur;

            if (dC < 0 || exp(-dC / T) > ur(rng)) {
                expr = move(cand);
                if (costC < bestCost) {
                    bestCost = costC;
                    bestExpr = expr;
                }
            }
        }
        T *= cooling;
    }

    return bestExpr;
}


//==== Main ====  
int main() {
    Data in;
    if (!readInput("public3.txt", in)) {
        cerr << "Failed to read input" << endl;
        return 1;
    }

    long long A = in.totalBlockArea;
    double deadRatio = 0.1;
    int W = (int)floor(sqrt((double)A * (1 + deadRatio)));
    int H = W;

    vector<pair<int, int>> pos;
    vector<int> rot;
    normalizeTall(in, rot);
    vector<int> polishExpr = initialLB2(in, W, H, pos, rot);

    vector<int> bestExpr = simulatedAnneal(in, polishExpr, rot, W, H, pos);


    vector<pair<int, int>> bestPos;
    bool ok = evalExpr(bestExpr, in, W, H, bestPos);
    if (!ok) {
        long long wl = computeHPWL(in, pos);
        writeOutput("output.out", wl, in, bestPos, rot);
        return 0;
    }

    long long bestWL = (long long)computeHPWL(in, bestPos);
    writeOutput("output.out", bestWL, in, bestPos, rot);

    // Output block positions
    cout << "Block positions:" << endl;
    for (size_t i = 0; i < bestPos.size(); ++i) {
        cout << "Block " << in.blocks[i].name << ": (" << bestPos[i].first << ", " << bestPos[i].second << ")" << endl;
    }

    return 0;
}


