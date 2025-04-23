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

//==== New Left-Bottom Initial Floorplan function ====  
void initialLB2(Data& in, int W, int H, vector<pair<int, int>>& pos, vector<int>& rot) {
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
        cout << "Placed block " << in.blocks[b1].name
            << " (width: " << in.blocks[b1].width
            << ", height: " << in.blocks[b1].height
            << ") at (0, " << y << ")" << endl;

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
                cout << "Placed block " << in.blocks[b2].name
                    << " (width: " << in.blocks[b2].width
                    << ", height: " << in.blocks[b2].height
                    << ") at (" << x << ", " << y << ")" << endl;
                cout << "Placed block " << in.blocks[best].name
                    << " (width: " << in.blocks[best].width
                    << ", height: " << in.blocks[best].height
                    << ") at (" << x << ", " << y + h2 << ")" << endl;
                x += w2;
            }
            else {
                // place b2 alone
                pos[b2] = { x, y };
                used[b2] = true;
                placedCount++;
                cout << "Placed block " << in.blocks[b2].name
                    << " (width: " << in.blocks[b2].width
                    << ", height: " << in.blocks[b2].height
                    << ") at (" << x << ", " << y << ")" << endl;
                x += w2;
            }
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
                cout << "Placed block " << in.blocks[best].name
                    << " (width: " << in.blocks[best].width
                    << ", height: " << in.blocks[best].height
                    << ") at (" << x << ", " << y << ")" << endl;
            }
        }

        // Step 5: move to next row
        y += h1;
    }

    // fallback: place any remaining blocks at origin
    for (int i = 0; i < n; ++i) {
        if (!used[i]) {
            pos[i] = { 0, 0 };
            rot[i] = 0;
            cout << "Fallback: Placed block " << in.blocks[i].name
                << " (width: " << in.blocks[i].width
                << ", height: " << in.blocks[i].height
                << ") at (0, 0)" << endl;
        }
    }
}


//==== Main ====  
int main() {
    Data in;
    if (!readInput("public1.txt", in)) {
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
    initialLB2(in, W, H, pos, rot);

    double wl = computeHPWL(in, pos);
    writeOutput("output.out", (long long)wl, in, pos, rot);
    return 0;
}

