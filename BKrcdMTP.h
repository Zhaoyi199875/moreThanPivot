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

class BKrcdMTP
{
public:
    BKrcdMTP(string const &fileName);
    ~BKrcdMTP();

    NeighborListArray** computeDegeneracyOrderArray();

    void fillInPandX( int vertex, int orderNumber,
                              int* vertexSets, int* vertexLookup, 
                              NeighborListArray** orderingArray,
                              int** neighborsInP, int* numNeighbors,
                              int* pBeginX, int *pBeginP, int *pBeginR, 
                              int* pNewBeginX, int* pNewBeginP, int *pNewBeginR,
                              long& scount, long& kcount, long& noise, long& dive, 
                              int& pivot, int& boundSize, int* degreeInPArray);

    void BKReverseRecursive(long* cliqueCount,
                                    list<int> &partialClique, 
                                    int* vertexSets, int* vertexLookup,
                                    int** neighborsInP, int* numNeighbors,
                                    int beginX, int beginP, int beginR);

    bool findMostNonNeighbors ( int* vertexSets, int *vertexLookup, 
                                   int **neighborsInP, int *numNeighbors,
                                   int beginX, int beginP, int beginR,
                                   int *pivot);
    
    bool existCommonNeighborOfPinX( int* vertexSets, int *vertexLookup, 
                                       int** neighborsInP, int* numNeighbors,
                                       int beginX, int beginP, int beginR);

    void moveToR( int vertex, 
                     int* vertexSets, int* vertexLookup, 
                     int **neighborsInP, int *numNeighbors,
                     int* pBeginX, int *pBeginP, int *pBeginR, 
                     int* pNewBeginX, int* pNewBeginP, int *pNewBeginR);

    void moveFromRToX( int vertex, 
                          int* vertexSets, int* vertexLookup,               
                          int**  neighborsInP, int* numNeighbors,
                          int* pBeginX, int* pBeginP, int* pBeginR );

    void listAllMaximalCliquesDgcyRecursive(long* cliqueCount,
                                               list<int> &partialClique, 
                                               int* vertexSets, int* vertexLookup,
                                               int** neighborsInP, int* numNeighbors,
                                               int beginX, int beginP, int beginR, 
                                               int pivot, int boundSize, int currentDepth, int* degreeInPArray);

    void moveToRDgcy( int vertex, 
                         int* vertexSets, int* vertexLookup, 
                         int** neighborsInP, int* numNeighbors,
                         int* pBeginX, int *pBeginP, int *pBeginR, 
                         int* pNewBeginX, int* pNewBeginP, int *pNewBeginR,
                         int& newPivot, int& newBoundSize, int boundSize, int* degreeInPArray);

    void findBestPivotNonNeighborsDgcy( int** pivotNonNeighbors, int* numNonNeighbors,
                                                int* vertexSets, int* vertexLookup,
                                                int** neighborsInP, int* numNeighbors,
                                                int beginX, int beginP, int beginR, int pivot, int boundSize, int currentDepth, int* degreeInPArray);

    void moveFromRToXDgcy( int vertex, 
                                    int* vertexSets, int* vertexLookup, 
                                    int* pBeginX, int* pBeginP, int* pBeginR );

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
