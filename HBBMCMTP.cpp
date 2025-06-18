// local includes
#include "HBBMCMTP.h"
#include <queue>
#include <functional>

using namespace std;

// static unsigned long numberOfRecurrciveCall=0;

HBBMCMTP::HBBMCMTP(string const &fileName){
    ifstream instream(fileName.c_str());
    if (instream.good() && !instream.eof())
        instream >> n;
    else {
        fprintf(stderr, "problem with line 1 in input file\n");
        exit(1);
    }

    if (instream.good() && !instream.eof())
        instream >> m;
    else {

        fprintf(stderr, "problem with line 2 in input file\n");
        exit(1);
    }
    
    adjList.resize(n);
    int u, v;
    int i = 0;
    while(i < m)
    {
        char comma;
        if (instream.good() && !instream.eof()) {
            instream >> u >> comma >> v;
        } else {
            fprintf(stderr, "problem with line %d in input file\n", i+2);
            exit(1);
        }
        assert(u < n && u > -1);
        assert(v < n && v > -1);
        if(u==v)
            fprintf(stderr, "Detected loop %d->%d\n", u, v);
        assert(u != v);

        adjList[u].push_back(v);

        i++;
    }

    vertexSets = (int*)malloc(n * sizeof(int));
    vertexLookup = (int*)malloc(n * sizeof(int));
    neighborsInP = (int**)malloc(n * sizeof(int*));
    numNeighbors = (int*)malloc(n * sizeof(int));
    inverseNeighborsInP = (int**)calloc(n, sizeof(int*));
    tempNeighbors = (int*)calloc(n, sizeof(int));
    tempInverseNeighbors = (int*)calloc(n, sizeof(int));
    degreeInPArray = (int*)malloc(n * sizeof(int));
    i = 0;
    while(i<n)
    {
        vertexSets[i] = i;
        vertexLookup[i] = i;
        neighborsInP[i] = (int*)malloc(1 * sizeof(int));
        degreeInPArray[i] = 0;
        tempNeighbors[i] = 0;
        tempInverseNeighbors[i] = 0;
        numNeighbors[i] = 1;
        i++;
    }
}

HBBMCMTP::~HBBMCMTP() {
    free(vertexSets);
    free(vertexLookup);
    for(int i=0; i<n; i++)
    {
        free(neighborsInP[i]);
        free(inverseNeighborsInP[i]);
        delete orderingArray[i];
    }
    free(neighborsInP);
    free(numNeighbors);
    free(inverseNeighborsInP);
    free(tempInverseNeighbors);
    free(degreeInPArray);
    delete[] orderingArray;

}

NeighborListArray** HBBMCMTP::computeDegeneracyOrderArray(){
    vector<NeighborList> vOrdering(n);

    int i = 0;

    int degeneracy = 0;
    
    // array of lists of vertices, indexed by degree
    vector<list<int>> verticesByDegree(n);

    // array of lists of vertices, indexed by degree
    vector<list<int>::iterator> vertexLocator(n);

    vector<int> degree(n);

    // fill each cell of degree lookup table
    // then use that degree to populate the 
    // lists of vertices indexed by degree

    for(i=0; i<n; i++)
    {
        degree[i] = adjList[i].size();
        verticesByDegree[degree[i]].push_front(i);
        vertexLocator[i] = verticesByDegree[degree[i]].begin();
    }

    int currentDegree = 0;

    int numVerticesRemoved = 0;

    while(numVerticesRemoved < n)
    {
        if(!verticesByDegree[currentDegree].empty())
        {
            degeneracy = max(degeneracy,currentDegree);
            
            int const vertex = verticesByDegree[currentDegree].front();
            verticesByDegree[currentDegree].pop_front();

            vOrdering[vertex].vertex = vertex;
            vOrdering[vertex].orderNumber = numVerticesRemoved;

            degree[vertex] = -1;

            list<int> const &neighborList = adjList[vertex];

            for(int const neighbor : neighborList)
            {
                if(degree[neighbor]!=-1)
                {
                    verticesByDegree[degree[neighbor]].erase(vertexLocator[neighbor]);
                    (vOrdering[vertex].later).push_back(neighbor);

                    degree[neighbor]--;

                    if(degree[neighbor] != -1)
                    {
                        verticesByDegree[degree[neighbor]].push_front(neighbor);
                        vertexLocator[neighbor] = verticesByDegree[degree[neighbor]].begin();
                    }
                }
                else
                {
                    vOrdering[vertex].earlier.push_back(neighbor);
                }
            }

            numVerticesRemoved++;
            currentDegree = 0;
        }
        else
        {
            currentDegree++;
        }

    }

    // orderingArray = (NeighborListArray**)calloc(n, sizeof(NeighborListArray*));
    orderingArray = new NeighborListArray*[n];
    for(i = 0; i<n;i++)
    {
        orderingArray[i] = new NeighborListArray();
        orderingArray[i]->vertex = vOrdering[i].vertex;
        orderingArray[i]->orderNumber = vOrdering[i].orderNumber;

        orderingArray[i]->laterDegree = vOrdering[i].later.size();
        orderingArray[i]->later.resize(orderingArray[i]->laterDegree);

        int j=0;
        for(int const laterNeighbor : vOrdering[i].later)
        {
            orderingArray[i]->later[j++] = laterNeighbor;
        }

        orderingArray[i]->earlierDegree = vOrdering[i].earlier.size();
        orderingArray[i]->earlier.resize(orderingArray[i]->earlierDegree);

        j=0;
        for (int const earlierNeighbor : vOrdering[i].earlier)
        {
            orderingArray[i]->earlier[j++] = earlierNeighbor;
        }
    }

    return orderingArray;
}

void HBBMCMTP::fillInPandXForRecursiveCall( int vertex, int orderNumber,
                                                   int* vertexSets, int* vertexLookup, 
                                                   NeighborListArray** orderingArray,
                                                   int** neighborsInP, int* numNeighbors,
                                                   int* pBeginX, int *pBeginP, int *pBeginR, 
                                                   int* pNewBeginX, int* pNewBeginP, int *pNewBeginR, int& pivot, int &boundSize,
                                                   int* degreeInPArray, int* tempNeighbors,int** inverseNeighborsInP,int* tempInverseNeighbors,int& isPlex)
{
        int vertexLocation = vertexLookup[vertex];

        (*pBeginR)--;
        vertexSets[vertexLocation] = vertexSets[*pBeginR];
        vertexLookup[vertexSets[*pBeginR]] = vertexLocation;
        vertexSets[*pBeginR] = vertex;
        vertexLookup[vertex] = *pBeginR;

        *pNewBeginR = *pBeginR;
        *pNewBeginP = *pBeginR;

        int j = 0;
        while(j<orderingArray[orderNumber]->laterDegree)
        {
            int neighbor = orderingArray[orderNumber]->later[j];
            int neighborLocation = vertexLookup[neighbor];

            (*pNewBeginP)--;

            vertexSets[neighborLocation] = vertexSets[*pNewBeginP];
            vertexLookup[vertexSets[*pNewBeginP]] = neighborLocation;
            vertexSets[*pNewBeginP] = neighbor;
            vertexLookup[neighbor] = *pNewBeginP;

            j++; 
        }

        *pNewBeginX = *pNewBeginP;

        j = 0;
        while(j<orderingArray[orderNumber]->earlierDegree)
        {
            int neighbor = orderingArray[orderNumber]->earlier[j];
            int neighborLocation = vertexLookup[neighbor];

            (*pNewBeginX)--;
            vertexSets[neighborLocation] = vertexSets[*pNewBeginX];
            vertexLookup[vertexSets[*pNewBeginX]] = neighborLocation;
            vertexSets[*pNewBeginX] = neighbor;
            vertexLookup[neighbor] = *pNewBeginX;

            free(neighborsInP[neighbor]);
            neighborsInP[neighbor] = (int*)calloc(min(*pNewBeginR-*pNewBeginP,orderingArray[neighbor]->laterDegree), sizeof(int));
            numNeighbors[neighbor] = 0;
            
            int k = 0;
            while(k<orderingArray[neighbor]->laterDegree)
            {
                int laterNeighbor = orderingArray[neighbor]->later[k];
                int laterNeighborLocation = vertexLookup[laterNeighbor];
                if(laterNeighborLocation >= *pNewBeginP && laterNeighborLocation < *pNewBeginR)
                {
                    neighborsInP[neighbor][numNeighbors[neighbor]] = laterNeighbor;
                    numNeighbors[neighbor]++;
                }

                k++;
            }

            if(numNeighbors[neighbor]>boundSize)
            {
                boundSize = numNeighbors[neighbor];
                pivot = neighbor;
            }

            j++; 

        }

        j = *pNewBeginP;
        while(j<*pNewBeginR)
        {
            int vertexInP = vertexSets[j];
            numNeighbors[vertexInP] = 0;
            free(neighborsInP[vertexInP]);
            neighborsInP[vertexInP]=(int*)calloc( min( *pNewBeginR-*pNewBeginP, 
                                                 orderingArray[vertexInP]->laterDegree 
                                               + orderingArray[vertexInP]->earlierDegree), sizeof(int));

            free(inverseNeighborsInP[vertexInP]);
            inverseNeighborsInP[vertexInP]=(int*)calloc( min( *pNewBeginR-*pNewBeginP, 
                                                 orderingArray[vertexInP]->laterDegree 
                                               + orderingArray[vertexInP]->earlierDegree), sizeof(int));

            j++;
        }

        j = *pNewBeginP;
        while(j<*pNewBeginR)
        {
            int vertexInP = vertexSets[j];

            int k = 0;
            while(k<orderingArray[vertexInP]->laterDegree)
            {
                int laterNeighbor = orderingArray[vertexInP]->later[k];
                int laterNeighborLocation = vertexLookup[laterNeighbor];

                if(laterNeighborLocation >= *pNewBeginP && laterNeighborLocation < *pNewBeginR)
                {
                    neighborsInP[vertexInP][numNeighbors[vertexInP]] = laterNeighbor;
                    numNeighbors[vertexInP]++;
                    neighborsInP[laterNeighbor][numNeighbors[laterNeighbor]] = vertexInP;
                    numNeighbors[laterNeighbor]++;
                }

                k++;
            }

            j++;
        }

        bool valid3Plex = true;
        bool valid2Plex = true;

        if((*pNewBeginP - *pNewBeginX) != 0)
        {
            valid3Plex = false;
            valid2Plex = false;
        }

        j = *pNewBeginP;
        while(j<*pNewBeginR)
        {
            int vertex = vertexSets[j];
            degreeInPArray[vertex] = numNeighbors[vertex];
            tempNeighbors[vertex] = numNeighbors[vertex];
            int sizeP = *pNewBeginR - *pNewBeginP;
            int missing = sizeP - 1 - numNeighbors[vertex];

            if (missing > 2) {
                valid3Plex = false;
                valid2Plex = false;
            } else if (missing > 1) {
                valid2Plex = false;
            }
            
            if(numNeighbors[vertex]>boundSize)
            {
                boundSize = numNeighbors[vertex];
                pivot = vertex;
            }
            j++;
        }
        
        if (!valid3Plex) {
            isPlex = 0;
        } else if (valid2Plex) {
            isPlex = 2;
        } else {
            isPlex = 3;
        }
        
        if(isPlex == 3)
        {
            int size = *pNewBeginR - *pNewBeginP;
            int* isInverseNeighbors = (int*)malloc(size * sizeof(int));
            for (int i = 0; i < size; i++)
            {
                isInverseNeighbors[i] = -1;
            }
            for (int i = *pNewBeginP; i < *pNewBeginR; i++)
            {
                int v = vertexSets[i];
                int loc = vertexLookup[v];
                tempInverseNeighbors[v] = 0;
                isInverseNeighbors[loc-*pNewBeginP] = v;
                for (int j = 0; j < tempNeighbors[v]; j++)
                {
                    
                    int neighbor = neighborsInP[v][j];
                    int neighborLocation = vertexLookup[neighbor];
                    if (neighborLocation >= *pNewBeginP && neighborLocation < *pNewBeginR)
                    {
                        isInverseNeighbors[neighborLocation - *pNewBeginP] = v;
                    }
                }
                
                for (int j = 0; j < size; j++)
                {
                    if (isInverseNeighbors[j] != v)
                    {
                        int inverseNeighbor = vertexSets[j + *pNewBeginP];
                        inverseNeighborsInP[v][tempInverseNeighbors[v]] = inverseNeighbor;
                        tempInverseNeighbors[v]++;
                    }
                }
            }
            free(isInverseNeighbors);
        }

}

void HBBMCMTP::findBestPivotNonNeighbors( int** pivotNonNeighbors, int* numNonNeighbors,
                                                int* vertexSets, int* vertexLookup,
                                                int** neighborsInP, int* numNeighbors,
                                                int beginX, int beginP, int beginR,int &currentDepth, 
                                                int &pivot, int &boundSize, int* degreeInPArray)
{
//    clock_t clockStart = clock();
    if(currentDepth == 1)
    {
        *pivotNonNeighbors = (int*)malloc((beginR-beginP) *sizeof(int));
        memcpy(*pivotNonNeighbors, &vertexSets[beginP], (beginR-beginP)*sizeof(int));
        *numNonNeighbors = 0;

        queue<int> coverQueue;
        std::queue<int> candidateQueue;
        queue<int> pivotQueue;
        int *vertexFlagSet = (int*)calloc(beginR-beginX, sizeof(int));

        int initialPivot = pivot;
        int initialpivotIndex = vertexLookup[initialPivot];

        int unprocessedVertices = beginR - beginP;
        int threshold = 4;

        pivotQueue.push(initialPivot);

        vertexFlagSet[initialpivotIndex-beginX] = 3; // 3 means pivot
        unprocessedVertices--;

        int N = 2;
        for(int i = 0; i < N; i++)
        // while (!pivotQueue.empty())
        {
            int new_pivot = pivotQueue.front();
            pivotQueue.pop();

            int numPivotNeighbors = min(boundSize, numNeighbors[new_pivot]);

            int k = 0;
            while(k<numPivotNeighbors)
            {
                int neighbor = neighborsInP[new_pivot][k];
                int neighborLocation = vertexLookup[neighbor];

                if(neighborLocation >= beginP && neighborLocation < beginR)
                {
                    if(vertexFlagSet[neighborLocation-beginX] == 0)
                    {
                        coverQueue.push(neighbor);
                        vertexFlagSet[neighborLocation-beginX] = 1; // 1 means covered
                        (*pivotNonNeighbors)[neighborLocation-beginP] = -1;
                        unprocessedVertices--;
                    }
                }
                else
                {
                    break;
                }
                k++;
            }

            if (unprocessedVertices > threshold && i < N-1)
            {            
                while(!coverQueue.empty())
                {
                    int coverVertex = coverQueue.front();
                    coverQueue.pop();
                    int numCoverNeighbors = min(boundSize, numNeighbors[coverVertex]);

                    int j = 0;
                    while(j<numCoverNeighbors)
                    {
                        int candidate = neighborsInP[coverVertex][j];
                        int candidateLocation = vertexLookup[candidate];

                        if(candidateLocation >= beginP && candidateLocation < beginR)
                        {
                            int flag = vertexFlagSet[candidateLocation-beginX];
                            if(flag == 0)
                            {
                                candidateQueue.push(candidate);
                                vertexFlagSet[candidateLocation-beginX] = 2; // 2 means pillar vertex
                                degreeInPArray[candidate]--;
                                unprocessedVertices--;
                            }
                            else if(flag == 2)
                            {
                                degreeInPArray[candidate]--;
                            }
                        }
                        else
                        {
                            break;
                        }
                        j++;
                    }
                }
            }
            else
                break;

            int candidatePivot;
            bool candidatePivotFound = false;
            bool bestPivotFound = false;
            int maxIntersectionSize = 0;

            while(!candidateQueue.empty())
            {
                int vertex = candidateQueue.front();
                candidateQueue.pop();
                int numPotentialNeighbors = min(boundSize, numNeighbors[vertex]);

                int neighborCountLCS = degreeInPArray[vertex];
                int neighborCountLPS = 0;

                int n = 0;
                while(n<numPotentialNeighbors)
                {
                    int neighbor = neighborsInP[vertex][n];
                    int neighborLocation = vertexLookup[neighbor];

                    if(neighborLocation >= beginP && neighborLocation < beginR)
                    {
                        if(vertexFlagSet[neighborLocation-beginX] == 0)
                        {
                            coverQueue.push(neighbor);
                        }
                    }
                    else
                    {
                        break;
                    }

                    n++;
                }

                while(!coverQueue.empty())
                {
                    int coverVertex = coverQueue.front();
                    coverQueue.pop();
                    int numCoverNeighbors = min(boundSize, numNeighbors[coverVertex]);

                    int j = 0;
                    while(j<numCoverNeighbors)
                    {
                        int candidate = neighborsInP[coverVertex][j];
                        int candidateLocation = vertexLookup[candidate];

                        if(candidateLocation >= beginP && candidateLocation < beginR)
                        {
                            int flag = vertexFlagSet[candidateLocation-beginX];
                            if(flag == 0)
                            {
                                neighborCountLPS++;
                            }
                        }
                        else
                        {
                            break;
                        }
                        j++;
                    }
                }

                int neighborCount = neighborCountLCS - neighborCountLPS;
                if( neighborCount > maxIntersectionSize)
                {
                    maxIntersectionSize = neighborCount;
                    candidatePivot = vertex;
                    candidatePivotFound = true;
                    if (neighborCountLCS > unprocessedVertices - 1||neighborCountLPS<1) 
                    {
                        bestPivotFound = true;
                        break;
                    }
                }

            }

            // if(!candidatePivotFound)
            // {
            // int j = beginX;
            // while(j < beginP)
            // {
            //     if (vertexFlagSet[j-beginX] == 3)
            //     {
            //         j++;
            //         continue;
            //     }
                
            //     int vertex = vertexSets[j];
            //     int numPotentialNeighbors = min(boundSize, numNeighbors[vertex]);
            //     int neighborCount = 0;

            //     int n = 0;
            //     while(n<numPotentialNeighbors)
            //     {
            //         int neighbor = neighborsInP[vertex][n];
            //         int neighborLocation = vertexLookup[neighbor];

            //         if(neighborLocation >= beginP && neighborLocation < beginR)
            //         {
            //             if(vertexFlagSet[neighborLocation-beginX] == 0)
            //             {
            //                 neighborCount++;
            //             }
            //         }
            //         else
            //         {
            //             break;
            //         }

            //         n++;
            //     }

            //     if(neighborCount > maxIntersectionSize)
            //     {
            //         maxIntersectionSize = neighborCount;
            //         candidatePivot = vertex;
            //         candidatePivotFound = true;

            //         if (maxIntersectionSize > unprocessedVertices - 1) 
            //         {
            //             bestPivotFound = true;
            //             break;
            //         }
            //     }

            //     j++;
            // }
            // }

            // if (!bestPivotFound)
            // {
            // int j = beginP;
            // while(j < beginR)
            // {
            //     int flag = vertexFlagSet[j-beginX];
            //     if(flag == 1||flag == 3)
            //     {
            //         j++;
            //         continue;
            //     }
            //     int vertex = vertexSets[j];
            //     int neighborCount = degreeInPArray[vertex];

            //     if(neighborCount > maxIntersectionSize)
            //     {
            //         maxIntersectionSize = neighborCount;
            //         candidatePivot = vertex;
            //         candidatePivotFound = true;
            //     }

            //     j++;
            // }
            // }

            if(candidatePivotFound)
            {
                pivotQueue.push(candidatePivot);
                int loc = vertexLookup[candidatePivot];
                if(loc>=beginX && loc < beginP)
                {
                    vertexFlagSet[loc-beginX] = 3;
                }
                else if(loc>=beginP && loc < beginR)
                {
                    if (vertexFlagSet[loc-beginX] == 0)
                    {
                        vertexFlagSet[loc-beginX] = 3;
                        unprocessedVertices--;
                    }
                    else
                    {
                        vertexFlagSet[loc-beginX] = 3;
                    }
                }
            }
            else
                break;
        }

        int tmpV = (*pivotNonNeighbors)[0];
        int j = 1;
        while(j<beginR-beginP)
        {
            int vertex = (*pivotNonNeighbors)[j];
            if(vertex < 0)
            {
                j++;
                continue;
            }
            (*pivotNonNeighbors)[*numNonNeighbors] = vertex;
            (*numNonNeighbors)++;
            j++;
        }

        free(vertexFlagSet);

        if(tmpV < 0)
        {
            return;
        }
        (*pivotNonNeighbors)[*numNonNeighbors] = tmpV;
        (*numNonNeighbors)++;
        return;

        // for(int i = beginP; i < beginR; i++)
        // {
        //     int flag = vertexFlagSet[i-beginX];
        //     if(flag == 3 ||flag == 2 || flag == 0)
        //     {
        //         (*pivotNonNeighbors)[*numNonNeighbors] = vertexSets[i];
        //         (*numNonNeighbors)++;
        //     }
        // }

        // clock_t clockEnd = clock();
        // timeComputingPivot += (clockEnd - clockStart);
        // Free(degreeInPArray);
    }
    else
    {
        *pivotNonNeighbors = (int*)malloc((beginR-beginP)*sizeof(int));
        memcpy(*pivotNonNeighbors, &vertexSets[beginP], (beginR-beginP)*sizeof(int));

        // we will decrement numNonNeighbors as we find neighbors
        // *numNonNeighbors = beginR-beginP;
        *numNonNeighbors = 0;

        int numPivotNeighbors = min(boundSize, numNeighbors[pivot]);

        // mark the neighbors of pivot that are in P.
        int j = 0;
        while(j<numPivotNeighbors)
        {
            int neighbor = neighborsInP[pivot][j];
            int neighborLocation = vertexLookup[neighbor];

            if(neighborLocation >= beginP && neighborLocation < beginR)
            {
                (*pivotNonNeighbors)[neighborLocation-beginP] = -1;
            }
            else
            {
                break;
            }

            j++;
        }

        // j = 0;
        // while(j<*numNonNeighbors)
        // {
        //     int vertex = (*pivotNonNeighbors)[j];

        //     if(vertex == -1)
        //     {
        //         (*numNonNeighbors)--;
        //         (*pivotNonNeighbors)[j] = (*pivotNonNeighbors)[*numNonNeighbors];
        //         continue;
        //     }

        //     j++;
        // }

        int tmpV = (*pivotNonNeighbors)[0];
        j = 1;
        while(j<beginR-beginP)
        {
            int vertex = (*pivotNonNeighbors)[j];
            if(vertex < 0)
            {
                j++;
                continue;
            }
            (*pivotNonNeighbors)[*numNonNeighbors] = vertex;
            (*numNonNeighbors)++;
            j++;
        }
        if(tmpV < 0)
        {
            return;
        }
        (*pivotNonNeighbors)[*numNonNeighbors] = tmpV;
        (*numNonNeighbors)++;
        // clock_t clockEnd = clock();
    
        // timeComputingPivot += (clockEnd - clockStart);
        return;

    }
}

void HBBMCMTP::moveToR( int vertex, 
                               int* vertexSets, int* vertexLookup, 
                               int** neighborsInP, int* numNeighbors,
                               int* pBeginX, int *pBeginP, int *pBeginR, 
                               int* pNewBeginX, int* pNewBeginP, int *pNewBeginR, int &newPivot,
                               int &newBoundSize,int boundSize, 
                               int currentDepth, int* tempNeighbors, int &isPlex, int** inverseNeighborsInP, int* tempInverseNeighbors)
{
//    clock_t clockStart = clock();
        int vertexLocation = vertexLookup[vertex];

        (*pBeginR)--;
        vertexSets[vertexLocation] = vertexSets[*pBeginR];
        vertexLookup[vertexSets[*pBeginR]] = vertexLocation;
        vertexSets[*pBeginR] = vertex;
        vertexLookup[vertex] = *pBeginR;

        *pNewBeginX = *pBeginP;
        *pNewBeginP = *pBeginP;
        *pNewBeginR = *pBeginP;

        int j = *pBeginX;
        while(j<*pNewBeginX)
        {
            int neighbor = vertexSets[j];
            int neighborLocation = j;

            int incrementJ = 1;

            int numPotentialNeighbors = min(boundSize, numNeighbors[neighbor]);
            int k = 0;
            while(k<numPotentialNeighbors)
            {
                if(neighborsInP[neighbor][k] == vertex)
                {
                    (*pNewBeginX)--;
                    vertexSets[neighborLocation] = vertexSets[(*pNewBeginX)];
                    vertexLookup[vertexSets[(*pNewBeginX)]] = neighborLocation;
                    vertexSets[(*pNewBeginX)] = neighbor;
                    vertexLookup[neighbor] = (*pNewBeginX);
                    incrementJ=0;
                    break;
                }
                k++;
            }
            if(incrementJ) j++;
        }

        int numPotentialNeighbors = min(boundSize, numNeighbors[vertex]);
        j = 0;
        while(j <numPotentialNeighbors)
        {
            int neighbor = neighborsInP[vertex][j];
            int neighborLocation = vertexLookup[neighbor];
            if(neighborLocation >= *pBeginP && neighborLocation < *pBeginR)
            {
                vertexSets[neighborLocation] = vertexSets[(*pNewBeginR)];
                vertexLookup[vertexSets[(*pNewBeginR)]] = neighborLocation;
                vertexSets[(*pNewBeginR)] = neighbor;
                vertexLookup[neighbor] = (*pNewBeginR);
                (*pNewBeginR)++;
            }
            j++;
        }

        bool valid3Plex = true;
        bool valid2Plex = true;

        if((*pNewBeginP - *pNewBeginX) != 0)
        {
            valid3Plex = false;
            valid2Plex = false;
        }

        j = (*pNewBeginX);
        while(j < *pNewBeginR)
        {
            int thisVertex = vertexSets[j];

            int numPotentialNeighbors = min(boundSize, numNeighbors[thisVertex]);
            tempNeighbors[thisVertex] = 0;

            int k = 0;
            while(k < numPotentialNeighbors)
            {
                int neighbor = neighborsInP[thisVertex][k];
                int neighborLocation = vertexLookup[neighbor];
                if(neighborLocation >= *pNewBeginP && neighborLocation < *pNewBeginR)
                {
                    neighborsInP[thisVertex][k] = neighborsInP[thisVertex][tempNeighbors[thisVertex]];
                    neighborsInP[thisVertex][tempNeighbors[thisVertex]] = neighbor;
                    tempNeighbors[thisVertex]++;
                }
                k++;
            }
            
            if(tempNeighbors[thisVertex]>newBoundSize)
            {
                newBoundSize = tempNeighbors[thisVertex];
                newPivot = thisVertex;
            }

            if(j>=*pNewBeginP)
            {
                int degreeInP = tempNeighbors[thisVertex];
                int sizeP = *pNewBeginR - *pNewBeginP;
                int missing = sizeP - 1 - degreeInP; 
                if (missing > 2) {
                    valid3Plex = false;
                    valid2Plex = false;
                } else if (missing > 1) {
                    valid2Plex = false;
                }
            }

            j++;;
        }

        if (!valid3Plex) {
            isPlex = 0;
        } else if (valid2Plex) {
            isPlex = 2;
        } else {
            isPlex = 3;
        }

        if(isPlex == 3)
        {
            int size = *pNewBeginR - *pNewBeginP;
            int* isInverseNeighbors = (int*)malloc(size * sizeof(int));
            for (int i = 0; i < size; i++)
            {
                isInverseNeighbors[i] = -1;
            }
            for (int i = *pNewBeginP; i < *pNewBeginR; i++)
            {
                int v = vertexSets[i];
                int loc = vertexLookup[v];
                tempInverseNeighbors[v] = 0;
                isInverseNeighbors[loc-*pNewBeginP] = v;
                for (int j = 0; j < tempNeighbors[v]; j++)
                {
                    int neighbor = neighborsInP[v][j];
                    int neighborLocation = vertexLookup[neighbor];
                    if (neighborLocation >= *pNewBeginP && neighborLocation < *pNewBeginR)
                    {
                        isInverseNeighbors[neighborLocation - *pNewBeginP] = v;
                    }
                }
                
                for (int j = 0; j < size; j++)
                {
                    if (isInverseNeighbors[j] != v)
                    {
                        int inverseNeighbor = vertexSets[j + *pNewBeginP];
                        inverseNeighborsInP[v][tempInverseNeighbors[v]] = inverseNeighbor;
                        tempInverseNeighbors[v]++;
                    }
                }
                // cout<<"sizeofP: "<<size<<" inverseNeighbor: "<< tempInverseNeighbors[v]<<" Neighbors: "<<tempNeighbors[v]<<endl;
            }
            free(isInverseNeighbors);
        }
        if (newBoundSize == (*pNewBeginR - *pNewBeginP))
        {
            newPivot = -2;
            return;
        }
        return;
}

void HBBMCMTP::moveFromRToX( int vertex, 
                                    int* vertexSets, int* vertexLookup, 
                                    int* pBeginX, int* pBeginP, int* pBeginR )
{
    int vertexLocation = vertexLookup[vertex];

    vertexSets[vertexLocation] = vertexSets[*pBeginP];
    vertexLookup[vertexSets[*pBeginP]] = vertexLocation;
    vertexSets[*pBeginP] = vertex;
    vertexLookup[vertex] = *pBeginP;

    *pBeginP = *pBeginP + 1;
    *pBeginR = *pBeginR + 1;
}

void HBBMCMTP::listAllMaximalCliquesRecursive(long* cliqueCount,
                                               list<int> &partialClique, 
                                               int* vertexSets, int* vertexLookup,
                                               int** neighborsInP, int* numNeighbors,
                                               int beginX, int beginP, int beginR,int currentDepth, 
                                               int &pivot, int &boundSize, int* degreeInPArray, 
                                               int* tempNeighbors, int** inverseNeighborsInP, int* tempInverseNeighbors)
{
    // numberOfRecurrciveCall++;
    if(beginX >= beginP && beginP >= beginR)
    {
        (*cliqueCount)++;
        return;
    }

    if(beginP >= beginR)
        return;

    int* myCandidatesToIterateThrough;
    int numCandidatesToIterateThrough = beginR - beginP;

    // get the candidates to add to R to make a maximal clique
    bool findPivot = true;
    if (numCandidatesToIterateThrough<=2)
    {
        findPivot = false;
    }

    if(findPivot)
    {
        findBestPivotNonNeighbors( &myCandidatesToIterateThrough,
                                    &numCandidatesToIterateThrough,
                                    vertexSets, vertexLookup,
                                    neighborsInP, numNeighbors,
                                    beginX, beginP, beginR, currentDepth, 
                                    pivot, boundSize, degreeInPArray);    
    }
    else
    {
        myCandidatesToIterateThrough = (int*)malloc((beginR-beginP)*sizeof(int));
        memcpy(myCandidatesToIterateThrough, &vertexSets[beginP], (beginR-beginP)*sizeof(int));
    }

    if(numCandidatesToIterateThrough != 0)
    {
    int iterator = 0;
    while(iterator < numCandidatesToIterateThrough)
    {
        // vertex to be added to the partial clique
        int vertex = myCandidatesToIterateThrough[iterator];

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("%d ", vertex);
        #endif

        int newBeginX, newBeginP, newBeginR;

        int newBoundSize = -1;
        int newPivot = -1;
        int isPlex;

        moveToR( vertex, 
                vertexSets, vertexLookup, 
                neighborsInP, numNeighbors,
                &beginX, &beginP, &beginR, 
                &newBeginX, &newBeginP, &newBeginR, newPivot, 
                newBoundSize, boundSize, 
                currentDepth,tempNeighbors, isPlex, inverseNeighborsInP, tempInverseNeighbors);

        if(newPivot == -2)
        {
            if(iterator == numCandidatesToIterateThrough-1)
            {
                beginR = beginR + 1;
                break;
            }
            moveFromRToX( vertex, 
                        vertexSets, vertexLookup,
                        &beginX, &beginP, &beginR);
    
            iterator++;
            continue;
        }
        if(isPlex == 2)
        {
            // cout<<"isPlex: "<<isPlex<<endl;
            int sizeofP = newBeginR - newBeginP;
            vector<int> F, disconnectPairs;
            for (int i = newBeginP; i < newBeginR; i++) {
                int v = vertexSets[i];
                if (tempNeighbors[v] == sizeofP - 1) F.push_back(v);
                else disconnectPairs.push_back(v);
            }

            int Lsize = disconnectPairs.size() / 2;
            for (int mask = 0; mask < (1 << Lsize); ++mask) {
                list<int> clique(partialClique.begin(), partialClique.end());
                clique.push_back(vertex);
                for (int f : F) clique.push_back(f);
                for (int i = 0; i < Lsize; ++i) {
                    int idx = i * 2;
                    int u = disconnectPairs[idx], v = disconnectPairs[idx + 1];
                    clique.push_back((mask & (1 << i)) ? u : v);
                }
                (*cliqueCount)++;
            }
            // (*cliqueCount) = (*cliqueCount) + (1 << Lsize);

            if(iterator == numCandidatesToIterateThrough-1)
            {
                beginR = beginR + 1;
                break;
            }

            moveFromRToX( vertex, 
                        vertexSets, vertexLookup,
                        &beginX, &beginP, &beginR);

            iterator++;
            continue;
        }
        else if(isPlex == 3)
        {
            int size = newBeginR - newBeginP;
            vector<int> F;
            for (int i = newBeginP; i < newBeginR; i++) {
                int v = vertexSets[i];
                if (tempInverseNeighbors[v] == 0) F.push_back(v);
            }
            
            bool isCycle = false;
            bool isPath = false;
            vector<int> visited(size, 0);
            vector<vector<int>> paths;
            vector<vector<int>> cycles;
            
            for (int i = newBeginP; i < newBeginR; i++) {
                int v = vertexSets[i];
                if (visited[i - newBeginP] == 1) continue;
            
                if (tempInverseNeighbors[v] == 1) {
                    isPath = true;
                    vector<int> path;
                    int current = v, prev = -1;
            
                    while (true) {
                        path.push_back(current);
                        visited[vertexLookup[current] - newBeginP] = 1;
            
                        int next = -1;
                        for (int j = 0; j < tempInverseNeighbors[current]; j++) {
                            int neighbor = inverseNeighborsInP[current][j];
                            if (visited[vertexLookup[neighbor] - newBeginP] == 0 && neighbor != prev) {
                                next = neighbor;
                                break;
                            }
                        }
            
                        if (next == -1 || visited[vertexLookup[next] - newBeginP] == 1) break;
                        prev = current;
                        current = next;
                    }
                    paths.push_back(path);
                }
            }
            
            for (int i = newBeginP; i < newBeginR; i++) {
                int v = vertexSets[i];
                int loc = vertexLookup[v] - newBeginP;
                if (visited[loc] == 1 || tempInverseNeighbors[v] < 2) continue;
            
                isCycle = true;
                vector<int> cycle;
                int current = v, prev = -1;
            
                while (true) {
                    cycle.push_back(current);
                    visited[vertexLookup[current] - newBeginP] = 1;
            
                    int next = -1;
                    for (int j = 0; j < tempInverseNeighbors[current]; j++) {
                        int neighbor = inverseNeighborsInP[current][j];
                        if (neighbor != prev && visited[vertexLookup[neighbor] - newBeginP] == 0) {
                            next = neighbor;
                            break;
                        }
                    }
            
                    if (next == -1 || visited[vertexLookup[next] - newBeginP] == 1) break;
            
                    prev = current;
                    current = next;
                }
                cycles.push_back(cycle);
            }
            
            vector<vector<vector<int>>> allCliqueSets;
            for (auto& p : paths) {
                vector<vector<int>> cliquesInPath;
                vector<int> p0 = {p[0]};
                vector<int> p1 = {p[1]};
                EnumRec(cliquesInPath, p0, p);
                EnumRec(cliquesInPath, p1, p);
                allCliqueSets.push_back(cliquesInPath);
            }
            
            for (auto& c : cycles) {
                vector<vector<int>> cliquesInCycle = enumFromCycle(c);
                allCliqueSets.push_back(cliquesInCycle);
            }
            
            list<int> baseClique(partialClique.begin(), partialClique.end());
            for (int f : F) baseClique.push_back(f);
            
            function<void(int, list<int>&)> combineCliques = [&](int depth, list<int>& current) {
                if (depth == allCliqueSets.size()) {
                    list<int> fullClique = current;
                    (*cliqueCount)++;
                    return;
                }
            
                for (const auto& option : allCliqueSets[depth]) {
                    for (int x : option) current.push_back(x);
                    combineCliques(depth + 1, current);
                    for (size_t i = 0; i < option.size(); ++i) current.pop_back();
                }
            };
            
            if (!allCliqueSets.empty()) {
                combineCliques(0, baseClique);
            } else {
                (*cliqueCount)++;
            }
            

            if(iterator == numCandidatesToIterateThrough-1)
            {
                beginR = beginR + 1;
                break;
            }

            moveFromRToX( vertex, 
                        vertexSets, vertexLookup,
                        &beginX, &beginP, &beginR);
            iterator++;
            continue;
        }

        partialClique.push_back(vertex);

        listAllMaximalCliquesRecursive(cliqueCount,
                                    partialClique, 
                                    vertexSets, vertexLookup,
                                    neighborsInP, numNeighbors,
                                    newBeginX, newBeginP, newBeginR, currentDepth+1, 
                                    newPivot,newBoundSize, degreeInPArray, 
                                    tempNeighbors, inverseNeighborsInP, tempInverseNeighbors);

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("b ");
        #endif

        partialClique.pop_back();

        if(iterator == numCandidatesToIterateThrough-1)
        {
            beginR = beginR + 1;
            break;
        }

        moveFromRToX( vertex, 
                    vertexSets, vertexLookup,
                    &beginX, &beginP, &beginR);

        iterator++;
    }

    iterator = 0;

////    clock_t clockStart = clock();
    while(iterator < numCandidatesToIterateThrough-1)
    {
        int vertex = myCandidatesToIterateThrough[iterator];
        int vertexLocation = vertexLookup[vertex];

        beginP--;
        vertexSets[vertexLocation] = vertexSets[beginP];
        vertexSets[beginP] = vertex;
        vertexLookup[vertex] = beginP;
        vertexLookup[vertexSets[vertexLocation]] = vertexLocation;

        iterator++;
    }
////    clock_t clockEnd = clock();
////    timeMovingXToP += (clockEnd - clockStart);
    }

    free(myCandidatesToIterateThrough);
}

void HBBMCMTP::EnumRec(vector<vector<int>>& S, vector<int>& Sp, vector<int>& p) {
    int i = Sp.back();  // i 是当前的顶点编号
    int idx = -1;
    // 找到 i 在 p 中的位置
    for (int j = 0; j < p.size(); ++j) {
        if (p[j] == i) {
            idx = j;
            break;
        }
    }
    if (idx == -1) return;

    if (idx + 2 >= p.size()) {
        S.push_back(Sp);  // 已经不能扩展，加入结果集
        return;
    }

    // Expand vi+2
    vector<int> Sp1 = Sp;
    Sp1.push_back(p[idx + 2]);
    EnumRec(S, Sp1, p);

    // Expand vi+3 if possible
    if (idx + 3 < p.size()) {
        vector<int> Sp2 = Sp;
        Sp2.push_back(p[idx + 3]);
        EnumRec(S, Sp2, p);
    }
}

vector<vector<int>> HBBMCMTP::enumFromCycle(const vector<int>& c) {
    int n = c.size();
    vector<vector<int>> result;
    if (n == 3) {
        for (int i = 0; i < 3; ++i)
            result.push_back({c[i]});
        return result;
    }
    if (n == 4) {
        result.push_back({c[0], c[2]});
        result.push_back({c[1], c[3]});
        return result;
    }
    if (n == 5) {
        result.push_back({c[0], c[2]});
        result.push_back({c[0], c[3]});
        result.push_back({c[1], c[3]});
        result.push_back({c[1], c[4]});
        result.push_back({c[2], c[4]});
        return result;
    }
    // |c| >= 6
    vector<vector<int>> S;

    // 从 v1 开始 path1: {v1, ..., v_{|c|-1}}
    vector<int> path1(c.begin(), c.end() - 1);
    vector<int> c0;
    c0.push_back(c[0]);
    vector<int> c1;
    c1.push_back(c[1]);

    EnumRec(S, c0, path1);

    // 从 v2 开始 path2: {v2, ..., v_{|c|}}
    vector<int> path2(c.begin() + 1, c.end());
    EnumRec(S, c1, path2);

    // 从 {v|c|, v3} 开始 path3: {v3, ..., v_{|c|-2}}
    vector<int> path3;
    for (int i = 2; i < c.size() - 2; ++i)
        path3.push_back(c[i]);
    if (!path3.empty()) {
        vector<int> start = {c.back(), c[2]};
        EnumRec(S, start, path3);
    }

    return S;
}

long HBBMCMTP::run(){
    clock_t start, mid, end;
    long cliqueCount = 0;
    list<int> partialClique;
    int* candidates;
    int numCandidates;

    start  = clock();
    int beginX = 0;
    int beginP = 0;
    int beginR = n;
    int currentDepth = 1;
    orderingArray = computeDegeneracyOrderArray();
    for(int i = 0; i<n ;i++)
    {
        int vertex = (int)orderingArray[i]->vertex;
        int newBeginX, newBeginP, newBeginR;

        int isPlex;
        int pivot=-1;
        int boundSize=-1;
        fillInPandXForRecursiveCall( i, vertex, 
                                    vertexSets, vertexLookup, 
                                    orderingArray,
                                    neighborsInP, numNeighbors,
                                    &beginX, &beginP, &beginR, 
                                    &newBeginX, &newBeginP, &newBeginR, pivot, boundSize, degreeInPArray,
                                    tempNeighbors, inverseNeighborsInP,tempInverseNeighbors, isPlex);

        if(isPlex == 2)
        {
            int sizeofP = newBeginR - newBeginP;
            vector<int> F, disconnectPairs;
            for (int j = newBeginP; j < newBeginR; j++) {
                int v = vertexSets[j];
                if (tempNeighbors[v] == sizeofP - 1) F.push_back(v);
                else disconnectPairs.push_back(v);
            }

            int Lsize = disconnectPairs.size() / 2;
            for (int mask = 0; mask < (1 << Lsize); ++mask) {
                list<int> clique(partialClique.begin(), partialClique.end());
                clique.push_back(vertex);
                for (int f : F) clique.push_back(f);
                for (int i = 0; i < Lsize; ++i) {
                    int idx = i * 2;
                    int u = disconnectPairs[idx], v = disconnectPairs[idx + 1];
                    clique.push_back((mask & (1 << i)) ? u : v);
                }
                (cliqueCount)++;
            }

            beginR = beginR + 1;
            continue;
        }
        else if(isPlex == 3)
        {
            int size = newBeginR - newBeginP;
            vector<int> F;
            for (int i = newBeginP; i < newBeginR; i++) {
                int v = vertexSets[i];
                if (tempInverseNeighbors[v] == 0) F.push_back(v);
            }
            
            bool isCycle = false;
            bool isPath = false;
            vector<int> visited(size, 0);
            vector<vector<int>> paths;
            vector<vector<int>> cycles;
            
            for (int i = newBeginP; i < newBeginR; i++) {
                int v = vertexSets[i];
                if (visited[i - newBeginP] == 1) continue;
            
                if (tempInverseNeighbors[v] == 1) {
                    isPath = true;
                    vector<int> path;
                    int current = v, prev = -1;
            
                    while (true) {
                        path.push_back(current);
                        visited[vertexLookup[current] - newBeginP] = 1;
            
                        int next = -1;
                        for (int j = 0; j < tempInverseNeighbors[current]; j++) {
                            int neighbor = inverseNeighborsInP[current][j];
                            if (visited[vertexLookup[neighbor] - newBeginP] == 0 && neighbor != prev) {
                                next = neighbor;
                                break;
                            }
                        }
            
                        if (next == -1 || visited[vertexLookup[next] - newBeginP] == 1) break;
                        prev = current;
                        current = next;
                    }
                    paths.push_back(path);
                }
            }
            
            for (int i = newBeginP; i < newBeginR; i++) {
                int v = vertexSets[i];
                int loc = vertexLookup[v] - newBeginP;
                if (visited[loc] == 1 || tempInverseNeighbors[v] < 2) continue;
            
                isCycle = true;
                vector<int> cycle;
                int current = v, prev = -1;
            
                while (true) {
                    cycle.push_back(current);
                    visited[vertexLookup[current] - newBeginP] = 1;
            
                    int next = -1;
                    for (int j = 0; j < tempInverseNeighbors[current]; j++) {
                        int neighbor = inverseNeighborsInP[current][j];
                        if (neighbor != prev && visited[vertexLookup[neighbor] - newBeginP] == 0) {
                            next = neighbor;
                            break;
                        }
                    }
            
                    if (next == -1 || visited[vertexLookup[next] - newBeginP] == 1) break;
            
                    prev = current;
                    current = next;
                }
                cycles.push_back(cycle);
            }
            
            vector<vector<vector<int>>> allCliqueSets;
            
            for (auto& p : paths) {
                vector<vector<int>> cliquesInPath;
                vector<int> p0 = {p[0]};
                vector<int> p1 = {p[1]};
                EnumRec(cliquesInPath, p0, p);
                EnumRec(cliquesInPath, p1, p);
                allCliqueSets.push_back(cliquesInPath);
            }
            
            for (auto& c : cycles) {
                vector<vector<int>> cliquesInCycle = enumFromCycle(c);
                allCliqueSets.push_back(cliquesInCycle);
            }
            
            list<int> baseClique(partialClique.begin(), partialClique.end());
            for (int f : F) baseClique.push_back(f);
            
            function<void(int, list<int>&)> combineCliques = [&](int depth, list<int>& current) {
                if (depth == allCliqueSets.size()) {
                    list<int> fullClique = current;
                    (cliqueCount)++;
                    return;
                }
            
                for (const auto& option : allCliqueSets[depth]) {
                    for (int x : option) current.push_back(x);
                    combineCliques(depth + 1, current);
                    for (size_t i = 0; i < option.size(); ++i) current.pop_back();
                }
            };
            
            if (!allCliqueSets.empty()) {
                combineCliques(0, baseClique);
            } else {
                (cliqueCount)++;
            }

            beginR = beginR + 1;
            
            continue;
        }
        partialClique.push_back(vertex);
        listAllMaximalCliquesRecursive(&cliqueCount,
                                    partialClique, 
                                    vertexSets, vertexLookup,
                                    neighborsInP, numNeighbors,
                                    newBeginX, newBeginP, newBeginR, currentDepth, 
                                    pivot,boundSize,degreeInPArray,
                                    tempNeighbors,inverseNeighborsInP, tempInverseNeighbors);   
        beginR = beginR + 1;
        partialClique.pop_back();
    }

    partialClique.clear();

    end = clock();
    cout << "******************total time count: " << (double)(end-start)/CLOCKS_PER_SEC << "s" << endl;
    cout << "******************Number of maximal cliques: " << cliqueCount << endl;
    // cout << "******************Number of recursive calls: " << numberOfRecurrciveCall << endl;

    return cliqueCount;
}

int main(int argc, char** argv) {
    string const filename = argv[1];
    HBBMCMTP alg(filename);
    alg.run();
}