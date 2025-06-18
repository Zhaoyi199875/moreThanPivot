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

struct NeighborList
{
    NeighborList()
    : vertex(-1)
    , earlier()
    , later()
    , orderNumber(-1) {}

    int vertex; //!< the vertex that owns this neighbor list
    std::list<int> earlier; //!< a linked list of neighbors that come before this vertex in the ordering
    std::list<int> later; //!< a linked list of neighbors that come after this vertex in the ordering
    int orderNumber; //!< the position of this vertex in the ordering
};

typedef struct NeighborList NeighborList;

class NeighborListArray
{
public:
    NeighborListArray()
    : vertex(-1)
    , earlier()
    , earlierDegree(-1)
    , later()
    , laterDegree(-1)
    , orderNumber(-1) {}

    int vertex; //!< the vertex that owns this neighbor list
    std::vector<int> earlier; //!< an array of neighbors that come before this vertex in an ordering
    int earlierDegree; //!< the number of neighbors in earlier
    std::vector<int> later; //!< an array of neighbors that come after this vertex in an ordering
    int laterDegree; //!< an array of neighbors that come after this vertex in an ordering
    int orderNumber; //!< the position of this vertex in the ordering
};

class DegenMTP
{
public:
    DegenMTP(string const &fileName);
    ~DegenMTP();

    NeighborListArray** computeDegeneracyOrderArray();

    void fillInPandXForRecursiveCall( int vertex, int orderNumber,
                                                   int* vertexSets, int* vertexLookup, 
                                                   NeighborListArray** orderingArray,
                                                   int** neighborsInP, int* numNeighbors,
                                                   int* pBeginX, int *pBeginP, int *pBeginR, 
                                                   int* pNewBeginX, int* pNewBeginP, int *pNewBeginR, int& pivot, int& pivot_IntersectionSize,
                                                   int* degreeInPArray);

    void listAllMaximalCliquesRecursive(long* cliqueCount,
                                               list<int> &partialClique, 
                                               int* vertexSets, int* vertexLookup,
                                               int** neighborsInP, int* numNeighbors,
                                               int beginX, int beginP, int beginR,int currentDepth, 
                                               int &pivot, int &pivot_IntersectionSize, int* degreeInPArray);

    void findBestPivotNonNeighbors( int** pivotNonNeighbors, int* numNonNeighbors,
                                                int* vertexSets, int* vertexLookup,
                                                int** neighborsInP, int* numNeighbors,
                                                int beginX, int beginP, int beginR,int &currentDepth, 
                                                int &pivot, int &pivot_IntersectionSize, int* degreeInPArray);

    void moveToR( int vertex, 
                               int* vertexSets, int* vertexLookup, 
                               int** neighborsInP, int* numNeighbors,
                               int* pBeginX, int *pBeginP, int *pBeginR, 
                               int* pNewBeginX, int* pNewBeginP, int *pNewBeginR, int &newPivot,
                               int &newBoundsize,int boundSize, 
                               int currentDepth, int* degreeInPArray);

    void moveFromRToX( int vertex, 
                                        int* vertexSets, int* vertexLookup, 
                                        int* pBeginX, int* pBeginP, int* pBeginR);

    long run();

private:
    int n;
    int m;
    NeighborListArray** orderingArray;
    int* vertexSets;
    int* vertexLookup;
    int** neighborsInP;
    int* numNeighbors;
    vector<list<int>> adjList;
    int* degreeInPArray;
};
