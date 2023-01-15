#include <mpi.h>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

using namespace std;

vector<pair<int, int>> vComparators; // вектор компаратаров
typedef vector<pair<int, int>> typeVecPair;
vector<typeVecPair> vTackts; // вектор из тактов

int iCountProc, my_rank; // количество процессоров и ранг текущего процесса
int idxSort = 0; // индекс сортировки, по умолчанию 0 - координата x

int iCountPointsOnProc; // количество точек на процесс

struct Point {
    float coord[2];
    int index;
};

// определяем новый тип MPI для структуры Point
void CreateMPITypePoint(MPI_Datatype& MPI_Type_Point) {
    int CountFields = 2;
    int LengthField[2] = { 2, 1 };
    MPI_Datatype TypeField[2] = { MPI_FLOAT, MPI_INT };
    MPI_Aint OffsetField[2];
    OffsetField[0] = offsetof(Point, coord);
    OffsetField[1] = offsetof(Point, index);
    MPI_Type_create_struct(CountFields, LengthField, OffsetField, TypeField, &MPI_Type_Point);
    MPI_Type_commit(&MPI_Type_Point);
}

float x(int i, int j) {
    return (float)i * i + j * j;
}

float y(int i, int j) {
    return (float)i * i - j * j;
}

// генерация точек, n1, n2 - размер сетки, iCountProc -  число процессоров
void InitPoints(vector<Point>& vPoints, int n1, int n2, int iCountProc) {
    Point curPoint;
    
    // количество точек, которые надо добавить на последний процессор
    int iAdd = ((n1 * n2) % iCountProc == 0) ? 0 : iCountProc - (n1 * n2) % iCountProc;

    // генерируем точки для основной сетки
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            curPoint.index = n2 * i + j;
            curPoint.coord[0] = x(i, j);
            curPoint.coord[1] = y(i, j);
            vPoints.push_back(curPoint);
        }
    }
    // добавляем фиктивные точки для последнего процессора
    for (int i = 0; i < iAdd; i++) {
        curPoint.index = -1;
        curPoint.coord[0] = 0;
        curPoint.coord[1] = 0;
        vPoints.push_back(curPoint);
    }
};

// печать точек
void PrintPoints(const vector<Point>& vPoints) {
    //cout << "Points (size: " << vPoints.size() << "): " << "rank " << my_rank << endl;
    for (int i = 0; i < vPoints.size(); i++) {
        cout << vPoints[i].coord[0] << " " << vPoints[i].coord[1] << " " << vPoints[i].index << " - rank " << my_rank << endl;
    }
}

// печать компараторов
void PrintComparators(const typeVecPair& vComparators) {
    cout << "Comparators (size: " << vComparators.size() << "):" << endl;
    for (int i = 0; i < vComparators.size(); i++) {
        cout << vComparators[i].first << " " << vComparators[i].second << endl;
    }

}

// печать тактов
void PrintTackts(const vector<typeVecPair>& vTackts) {
    cout << "Tackts (size: " << vTackts.size() << "):" << endl;
    for (int i = 0; i < vTackts.size(); i++) {
        for (int j = 0; j < vTackts[i].size(); j++) cout << "(" << vTackts[i][j].first << " " << vTackts[i][j].second << ") ";
        cout << endl;
    }
}

// сравнение точек для стандартной функции сортировки
int ComparePoints(const void* First, const void* Second)
{
    Point FPoint = *((const Point*)First);
    Point SPoint = *((const Point*)Second);

    if (FPoint.coord[idxSort] < SPoint.coord[idxSort]) return -1;
    if (FPoint.coord[idxSort] > SPoint.coord[idxSort]) return 1;
    return 0;
}

// возвращает true, если первая точка больше второй; idxSort - индекс сортировки (0, если x; 1, если y) 
bool bGreater(Point First, Point Second, int idxSort) {
    return First.coord[idxSort] > Second.coord[idxSort];
}

// возвращает true, если первая точка меньше либо равна второй; idxSort - индекс сортировки (0, если x; 1, если y) 
bool bLess(Point First, Point Second, int idxSort) {
    return First.coord[idxSort] <= Second.coord[idxSort];
}

// Join - S - слияние массивов Бэтчер
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

// Sort - B - сортировка массива Бэтчер
void Sort(int iFirst, int iStep, int iCount) {
    if (iCount < 2) return;
    if (iCount == 2) {
        vComparators.push_back(make_pair(iFirst, iFirst + iStep));
        return;
    }
    
    int iCount1;
    if (iCount % 2 != 0) iCount1 = iCount / 2 + 1;
    else iCount1 = iCount / 2;

    Sort(iFirst, iStep, iCount1);
    Sort(iFirst + iStep * iCount1, iStep, iCount - iCount1);
    Join(iFirst, iFirst + iStep * iCount1, iStep, iCount1, iCount - iCount1);
 }

// проверям можно ли добавить компаратор в текущий такт
bool bAddComp(pair<int, int> Comp, typeVecPair vTact) {
    bool bOK = true;
    int n = vTact.size();
    int i = 0;
    while (bOK && i < n) {
        bOK = (Comp.first != vTact[i].first && Comp.first != vTact[i].second && Comp.second != vTact[i].first && Comp.second != vTact[i].second);
        i++;
    }
    return bOK;
}

// раскидываем компараторы по тактам
void CreateTackts(typeVecPair& vComps, vector<typeVecPair>& vTackts) {
    reverse(vComps.begin(), vComps.end()); // инвертируем вектор компараторов 
    
    vTackts.resize(1);
    vTackts[0].push_back(vComps[vComps.size()-1]);
    vComps.pop_back();
    
    while (!vComps.empty()) { // пока вектор компараторов не пуст
        int n = vTackts.size();
        int idxCurComps = vComps.size()-1;
        int i = n - 1;
        // определяем номер такта в который можно вставить текущий компаратор
        while (i >= 0 && bAddComp(vComps[idxCurComps], vTackts[i])) i--;
        i++; // номер такта в который нужно добавить текущий компаратор

        if (i >= n) { // если тактов не хватает, то создаем новый такт
            vTackts.resize(++n);
        }
        vTackts[i].push_back(vComps[idxCurComps]); // добавляем компаратор в нужный такт
        vComps.pop_back(); // удаляем текущий компаратор из вектора компараторов
    }
}

// Merge - слияние массивов между процессорами, книга стр. 152
void MergeOld(vector<Point>& vFirst, vector<Point>& vSecond, int rank1, int rank2) {
    int iSize = vFirst.size();
    vector<Point> vResult(iSize);
    
    if (my_rank == rank1) { // формирование массива, содержащего меньшие элементы массивов iFirst и iSecond
        for (int iF = 0, iS = 0, k = 0; k < iSize;) {
            if (bLess(vFirst[iF], vSecond[iS], idxSort)) vResult[k++] = vFirst[iF++];
            else vResult[k++] = vSecond[iS++];
        }
        vFirst = vResult;
    }
    else if (my_rank == rank2) { // формирование массива, содержащего большие элементы массивов iFirst и iSecond
        for (int iF = iSize - 1, iS = iSize - 1, k = iSize - 1; k >= 0;) {
            if (bGreater(vFirst[iF], vSecond[iS], idxSort)) vResult[k--] = vFirst[iF--];
            else vResult[k--] = vSecond[iS--];
        }
        vSecond = vResult;
    }
}

void Merge(vector<Point>& vFirst, vector<Point>& vSecond, int rank1, int rank2) {
    int iSize = vFirst.size();
    vector<Point> vResult(iSize);

    if (my_rank == rank1) { // формирование массива, содержащего меньшие элементы массивов iFirst и iSecond
        for (int iF = 0, iS = 0, k = 0; k < iSize;) {
            if (bLess(vFirst[iF], vSecond[iS], idxSort)) vResult[k++] = vFirst[iF++];
            else vResult[k++] = vSecond[iS++];
        }
        vFirst = vResult;
    }
    else if (my_rank == rank2) { // формирование массива, содержащего большие элементы массивов iFirst и iSecond
        int iF = 0, iS = 0, k = 0;
        // пропускаем младшие элементы, которые ушли в массив vFirst
        for (; k < iSize;) {
            if (bLess(vFirst[iF], vSecond[iS], idxSort)) {
                k++; iF++;
            }
            else {
                k++; iS++;
            }
        }
        k = 0;
        // формируем массив vSecond из оставшихся элементов
        for (; k < iSize;) {
            if (iF < iSize && iS < iSize) {
                if (bLess(vFirst[iF], vSecond[iS], idxSort)) vResult[k++] = vFirst[iF++];
                else vResult[k++] = vSecond[iS++];
            } 
            else 
                if (iF < iSize)  vResult[k++] = vFirst[iF++];
                else vResult[k++] = vSecond[iS++];
        }
        vSecond = vResult;
    }
}

// параллельная сортировка Бэтчера
void ParallelSort(int n1, int n2) {
    if (my_rank == 0) cout << n1 << ":" << n2 << endl;

    MPI_Status status;

    vector<Point> vPoints; // все точки
    InitPoints(vPoints, n1, n2, iCountProc); // инициализируем точки
    iCountPointsOnProc = vPoints.size() / iCountProc; // определяем количество точек на процесс

    vector<Point> vPointsOnProc(iCountPointsOnProc); // основные точки на процессе
    vector<Point> vBufOnProc(iCountPointsOnProc); // буфер для приема точек

    // определяем новый тип MPI для структуры Point
    MPI_Datatype MPI_Type_Point;
    CreateMPITypePoint(MPI_Type_Point);

    //if (my_rank == 0) cout << "iCountPointsOnProc:" << iCountPointsOnProc << endl;

    // рассылаем точки по процессам
    MPI_Scatter(vPoints.data(), iCountPointsOnProc, MPI_Type_Point, vPointsOnProc.data(), iCountPointsOnProc, MPI_Type_Point, 0, MPI_COMM_WORLD);

    double StartTime = MPI_Wtime(); // время старта сортировки

    // сортируем массивы на каждом процессоре
    qsort(vPointsOnProc.data(), vPointsOnProc.size(), sizeof(Point), ComparePoints);
    MPI_Barrier(MPI_COMM_WORLD);

    //PrintPoints(vPointsOnProc);
     //if (my_rank == 0) PrintPoints(vPoints); // печать точек

    Sort(0, 1, iCountProc); // формируем расписание сети слияния
    // if (my_rank == 0) PrintComparators(vComparators); // печать компараторов

    CreateTackts(vComparators, vTackts); // формируем такты
    reverse(vTackts.begin(), vTackts.end());

    int iCountTackts = vTackts.size(); // сохраняем число тактов

    //if (my_rank == 0) PrintTackts(vTackts); // печать тактов

    MPI_Barrier(MPI_COMM_WORLD); // все процессы готовы к параллельной сортировке

    while (!vTackts.empty()) { // пока не исчерпали все такты
        int iCurTackt = vTackts.size() - 1;
        int iSize = vTackts[iCurTackt].size(); // размер текущего такта

        for (int i = iSize - 1; i >= 0; i--) { // бежим по всем компараторам текущего такта
            int rank1 = vTackts[iCurTackt][i].first;
            int rank2 = vTackts[iCurTackt][i].second;
            vTackts[iCurTackt].pop_back(); // удаляем текущий компаратор

            // пересылка/прием массивов если my_rank == rank1 
            if (my_rank == rank1) {
                //cout << "(" << rank1 << ", " << rank2 << ")" << endl;
                MPI_Send(vPointsOnProc.data(), vPointsOnProc.size(), MPI_Type_Point, rank2, 0, MPI_COMM_WORLD);
                MPI_Recv(vBufOnProc.data(), vBufOnProc.size(), MPI_Type_Point, rank2, 0, MPI_COMM_WORLD, &status);
                Merge(vPointsOnProc, vBufOnProc, rank1, rank2); // слияние массивов
            }
            // пересылка/прием массивов если my_rank == rank2 
            if (my_rank == rank2) {
                MPI_Recv(vBufOnProc.data(), vBufOnProc.size(), MPI_Type_Point, rank1, 0, MPI_COMM_WORLD, &status);
                MPI_Send(vPointsOnProc.data(), vPointsOnProc.size(), MPI_Type_Point, rank1, 0, MPI_COMM_WORLD);
                Merge(vBufOnProc, vPointsOnProc, rank1, rank2); // слияние массивов
            }
        }
        MPI_Barrier(MPI_COMM_WORLD); // окончания работы такта
        //cout << "iCurTackt: " << iCurTackt << endl;
        vTackts.resize(iCurTackt);
    }

    //собираем отсортированный массив в vPoints на my_rank == 0
    //MPI_Gather(vPointsOnProc.data(), iCountPointsOnProc, MPI_Type_Point, vPoints.data(), iCountPointsOnProc, MPI_Type_Point, 0, MPI_COMM_WORLD); 
    MPI_Barrier(MPI_COMM_WORLD);

    double EndTime = MPI_Wtime(); // время окончания сортировки
    double SortTime = EndTime - StartTime; // время сортировки

    if (my_rank == 0) {
        cout << "Sort time: " << SortTime << endl;
        cout << "iCountTackts: " << iCountTackts << endl;
    }

    //PrintPoints(vPointsOnProc);
    //if (my_rank == 0) PrintPoints(vPoints);
}

// последовательная сортировка
void SerialSort(int n1, int n2) {
    cout << n1 << ":" << n2 << endl;

    vector<Point> vPoints; // все точки
    InitPoints(vPoints, n1, n2, iCountProc); // инициализируем точки
    
    double StartTime = MPI_Wtime(); // время старта сортировки
    // сортируем массивы на каждом процессоре
    qsort(vPoints.data(), vPoints.size(), sizeof(Point), ComparePoints);
    double EndTime = MPI_Wtime(); // время окончания сортировки
    double SortTime = EndTime - StartTime; // время сортировки
    cout << "Sort time: " << SortTime << endl;
}

int main(int argc, char* argv[]) {
   
   MPI_Init(&argc, &argv);  // инициализация MPI
 
    int n1 = 3, n2 = 2; // размер сетки по умолчанию
    
    if (argc > 2) {
        n1 = atoi(argv[1]);
        n2 = atoi(argv[2]);
    }
    
    // определяем количество процессоров и ранг текущего процесса
    MPI_Comm_size(MPI_COMM_WORLD, &iCountProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if (iCountProc > 1) ParallelSort(n1, n2);
    else SerialSort(n1, n2);

    MPI_Finalize(); // закрытие MPI
} 