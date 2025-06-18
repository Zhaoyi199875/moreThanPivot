#include <algorithm>
#include <vector>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <utility>
#include <list>
#include <ctime>
#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>
#include <unordered_set>


using namespace std;

class RMCEdegenMTP
{
public:
    RMCEdegenMTP(string const &fileName);
    ~RMCEdegenMTP();

    void computeDegeneracyOrder();
    long computeDegeneracyWithGlobalReduction(int *offset);

    void P2R(int vertex, 
             int* posX, int* posP, int* posR, 
             int* newPosX, int* newPosP, int* newPosR);

    void P2RandFindPivot(int vertex, 
                int* posX, int* posP, int* posR, 
                int* newPosX, int* newPosP, int* newPosR, 
                int* pivot, int* pivotIntersectionSize, int boundSize, long* cliqueCount, int currentDepth);

    void R2X(int vertex, 
             int* posX, int* posP, int* posR);

    void setCandidates(int** candidates, int* numCandidates,
                                   int posX, int posP, int posR, int pivot, int pivotIntersectionSize, int* vertexDegree,
                                   long* cliqueCount, int currentDepth);

    void bpRecur(long* cliqueCount,
                 list<int> &partialClique, 
                 int posX, int posP, int posR, int pivot, int pivotIntersectionSize, int* vertexDegree, int currentDepth);

    long run();

private:
    int n;
    int m;
    int degeneracy=-1;
    int* vals;
    int* inds;
    int* seps;
    int* seps2;
    int* isEdgeDeleted;
    // int* laterStart;
    int* orders;
    int* o2v;
    int* orderedVals;
    int* vertexArray;
    int* vertexPos;
    int** neighborsInP;
    // vector<vector<int>> neighborsInP;
    int* numNeighbors;
    int* tmpNumNeighbors; // record the instant number of neighbors for pruning
    bool* bfsParent;
    int* Xpruned;
    int* vertexDegree;
};

