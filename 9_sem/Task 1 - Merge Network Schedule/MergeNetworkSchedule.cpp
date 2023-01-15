#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <bitset>
#include <algorithm>

using namespace std;

vector<pair<int, int>> vComparators;

void Join(int iFirst1, int iFirst2, int iStep, int iCount1, int iCount2) {
    int iCountOdd1, iCountEven2, i;

    if (iCount1 * iCount2 < 1) return;
    if (iCount1 == 1 && iCount2 == 1) {
        vComparators.push_back(make_pair(iFirst1, iFirst2));
        return;
    }
    
    iCountOdd1 = iCount1 - (iCount1 / 2);
    iCountEven2 = iCount2 - (iCount2 / 2);

    Join(iFirst1, iFirst2, 2 * iStep, iCountOdd1, iCountEven2);
    Join(iFirst1 + iStep, iFirst2 + iStep, 2 * iStep, iCount1 - iCountOdd1, iCount2 - iCountEven2);
    
    for (i = 1; i < iCount1 - 1; i += 2) {
        vComparators.push_back(make_pair(iFirst1 + iStep * i, iFirst1 + iStep * (i + 1)));
    }
    if (iCount1 % 2 == 0) {
        vComparators.push_back(make_pair(iFirst1 + iStep * (iCount1 - 1), iFirst2));
        i = 1;
    }  
    else i = 0;

    for (; i < iCount2 - 1; i += 2) {
        vComparators.push_back(make_pair(iFirst2 + iStep * i, iFirst2 + iStep * (i + 1)));
    }

}

bool bOK = true;
void check() {
    for (int iSize = 1; iSize <= 24; iSize++) {
        for (int i = 0; i <= pow(2, iSize) - 1; i++) {
            string curSeq = bitset<24>(i).to_string();
            curSeq = curSeq.substr(24 - iSize);

            for (int p1 = 1; p1 <= iSize; p1++) {
                int p2 = iSize - p1;
                bool bOKFirst = is_sorted(begin(curSeq), begin(curSeq) + p1 - 1);
                bool bOKSecond = is_sorted(begin(curSeq) + p1, end(curSeq));

                if (bOKFirst && bOKSecond) {
                    vComparators.clear();

                    Join(0, p1, 1, p1, p2);

                    for (int j = 0; j < vComparators.size(); j++) {
                        if (((int)curSeq[vComparators[j].first]) > ((int)curSeq[vComparators[j].second]))
                            swap(curSeq[vComparators[j].first], curSeq[vComparators[j].second]);
                    }
                    if (!is_sorted(begin(curSeq), end(curSeq))) bOK = false;
                    //cout << curSeq << " p1 = " << p1 << " p2 = " << p2 << endl;
                }
            }
        }
    }
    if (bOK) cout << "Good job, bro:)!" << endl;
    else cout << "Oops, tests failed:c" << endl;
}

int CountTackts(vector<pair<int, int>> vCur, int iSize) {
    int iMax1, iMax2;
    vector<int>vTackts(iSize);

    for (int i = 0; i < vCur.size(); i++) {
        iMax1 = max(vTackts[vComparators[i].first], vTackts[vComparators[i].second]);
        vTackts[vComparators[i].first] = iMax1 + 1;
        vTackts[vComparators[i].second] = iMax1 + 1;
    }

    iMax2 = vTackts[0];
    for (int i = 1; i < iSize; i++) {
        if (vTackts[i] > iMax2) iMax2 = vTackts[i];
    }

    return iMax2;
 };


int main(int argc, char** argv)
{
    int p1, p2;
    vector<int> v1, v2;

    cout << argv[1] << " " <<  argv[2] << " " << "0" << endl;
    p1 = stoi(argv[1]);
    p2 = stoi(argv[2]);

    Join(0, p1, 1, p1, p2);
    
    for (int i = 0; i < vComparators.size(); i++) {
        cout << vComparators[i].first << " " << vComparators[i].second << endl;
    }
    cout << vComparators.size() << endl;
    cout << CountTackts(vComparators, p1 + p2) << endl;
    //check();


    return 0;
}
