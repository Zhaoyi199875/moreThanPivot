#include "BKrcdMTP.h"
#include <queue>

using namespace std;

// static unsigned long numberOfRecurrciveCall=0;

BKrcdMTP::BKrcdMTP(string const &fileName){
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
    degreeInPArray = (int*)malloc(n * sizeof(int));

    i = 0;
    while(i<n)
    {
        vertexSets[i] = i;
        vertexLookup[i] = i;
        degreeInPArray[i] = 0;
        neighborsInP[i] = (int*)malloc(1 * sizeof(int));
        numNeighbors[i] = 1;
        i++;
    }
}

BKrcdMTP::~BKrcdMTP() {
    free(vertexSets);
    free(vertexLookup);
    for(int i=0; i<n; i++)
    {
        free(neighborsInP[i]);
    }
    for(int i=0; i<n; i++)
    {
        delete orderingArray[i];
        // free(orderingArray[i]);
    }
    free(numNeighbors);
    // free(orderingArray);
    free(degreeInPArray);
    delete[] orderingArray;
}

NeighborListArray** BKrcdMTP::computeDegeneracyOrderArray(){
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

void BKrcdMTP::fillInPandX( int vertex, int orderNumber,
                              int* vertexSets, int* vertexLookup, 
                              NeighborListArray** orderingArray,
                              int** neighborsInP, int* numNeighbors,
                              int* pBeginX, int *pBeginP, int *pBeginR, 
                              int* pNewBeginX, int* pNewBeginP, int *pNewBeginR,
                              long& scount, long& kcount, long& noise, long& dive, int& pivot, int& boundSize, int* degreeInPArray)
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

        // reset numNeighbors and neighborsInP for this vertex
        j = *pNewBeginP;
        while(j<*pNewBeginR)
        {
            int vertexInP = vertexSets[j];
            numNeighbors[vertexInP] = 0;
            free(neighborsInP[vertexInP]);
            neighborsInP[vertexInP]=(int*)calloc( min( *pNewBeginR-*pNewBeginP, 
                                                 orderingArray[vertexInP]->laterDegree 
                                               + orderingArray[vertexInP]->earlierDegree), sizeof(int));
            j++;
        }

        // count neighbors in P, and fill in array of neighbors
        // in P
        j = *pNewBeginP;
        long long edges = 0;
        while(j<*pNewBeginR)
        {
            int vertexInP = vertexSets[j];

            int k = 0;
            while(k < orderingArray[vertexInP]->laterDegree)
            {
                int laterNeighbor = orderingArray[vertexInP]->later[k];
                int laterNeighborLocation = vertexLookup[laterNeighbor];

                if(laterNeighborLocation >= *pNewBeginP && laterNeighborLocation < *pNewBeginR)
                {
                    neighborsInP[vertexInP][numNeighbors[vertexInP]] = laterNeighbor;
                    numNeighbors[vertexInP]++;
                    neighborsInP[laterNeighbor][numNeighbors[laterNeighbor]] = vertexInP;
                    numNeighbors[laterNeighbor]++;
                    edges++;
                }

                k++;
            }

            j++;
        }

        int oldBeginR = *pNewBeginR;
        int sizeofP = *pNewBeginR - *pNewBeginP;
        int vers_max_degree = 0;
        int max_degree = sizeofP - 1;

        scount = 0;
        kcount = 0;
        noise = 0;
        dive = 0;

        j = *pNewBeginP;
        while(j < *pNewBeginR)
        {
            int vertexInP = vertexSets[j];

            if (numNeighbors[vertexInP] == max_degree) {
                scount++;
                vers_max_degree += 1;
                int vertexLocation = vertexLookup[vertexInP];

                (*pNewBeginR)--;

                vertexSets[vertexLocation] = vertexSets[*pNewBeginR];
                vertexLookup[vertexSets[*pNewBeginR]] = vertexLocation;
                vertexSets[*pNewBeginR] = vertexInP;
                vertexLookup[vertexInP] = *pNewBeginR;

                continue;
            }

            ++j;
        }

        sizeofP = *pNewBeginR - *pNewBeginP;
        kcount = sizeofP;

        j = *pNewBeginP;
        while(j < *pNewBeginR) 
        {
            int vertexInP = vertexSets[j];
            int k = 0;
            while (k < numNeighbors[vertexInP]) 
            {
                if ( vertexLookup[neighborsInP[vertexInP][k]] >= *pNewBeginR ) 
                {
                    swap(neighborsInP[vertexInP][k], neighborsInP[vertexInP][--numNeighbors[vertexInP]]);
                    continue;
                }
                ++k;
            }

            noise += numNeighbors[vertexInP];
            dive = max((int)dive, numNeighbors[vertexInP]);

            ++j;
        }

        // swap earlier neighbors of vertex into X section of vertexSets
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

            // fill in NeighborsInP
            int k = 0;
            int degreeInR = 0;
            while(k<orderingArray[neighbor]->laterDegree)
            {
                int laterNeighbor = orderingArray[neighbor]->later[k];
                int laterNeighborLocation = vertexLookup[laterNeighbor];
                if(laterNeighborLocation >= *pNewBeginP && laterNeighborLocation < *pNewBeginR)
                {
                    neighborsInP[neighbor][numNeighbors[neighbor]] = laterNeighbor;
                    numNeighbors[neighbor]++;
                } else if (laterNeighborLocation >= *pNewBeginR && laterNeighborLocation < oldBeginR) {
                    degreeInR += 1;
                }

                k++;
            }

            if ( degreeInR < vers_max_degree ) {
                (*pNewBeginX)++;
            }

            j++; 

        }

        j = *pNewBeginX;
        while (j < *pNewBeginR) {
            int vertex = vertexSets[j];
            degreeInPArray[vertex] = numNeighbors[vertex];
            if (numNeighbors[vertex] > boundSize) {
                pivot = vertex;
                boundSize = numNeighbors[vertex];
            }
            j++;
        }
}

bool BKrcdMTP::findMostNonNeighbors ( int* vertexSets, int *vertexLookup, 
                                   int **neighborsInP, int *numNeighbors,
                                   int beginX, int beginP, int beginR,
                                   int *pivot)
{
    int j       = beginP;
    int degree  = beginR-beginP;
    int sizeofP = beginR-beginP;

    while (j < beginR) {
        int numNeighborsInP = 0;
        int vertex = vertexSets[j];
        //int numPotentialNeighbors = min(sizeofP, numNeighbors[vertex]);
        int numPotentialNeighbors = numNeighbors[vertex];

        int k = 0;
        while(k<numPotentialNeighbors)
        {
            int neighbor = neighborsInP[vertex][k];
            int neighborLocation = vertexLookup[neighbor];

            if(neighborLocation >= beginP && neighborLocation < beginR)
                numNeighborsInP++;
            k++;
        }

        if (numNeighborsInP < degree) {
            *pivot = vertexSets[j];
            degree = numNeighborsInP;
        } 

        ++j;
    }

    if (degree == beginR-beginP-1)
        return true;

    return false;
}

bool BKrcdMTP::existCommonNeighborOfPinX( int* vertexSets, int *vertexLookup, 
                                       int** neighborsInP, int* numNeighbors,
                                       int beginX, int beginP, int beginR)
{
    int j = beginX, sizeofP = beginR-beginP;
    while (j < beginP) {
        int numNeighborsInP = 0;
        int vertex = vertexSets[j];
        int numPotentialNeighbors = numNeighbors[vertex];

        int k = 0;
        while(k < numPotentialNeighbors)
        {
            int neighbor = neighborsInP[vertex][k];
            int neighborLocation = vertexLookup[neighbor];
            if(neighborLocation >= beginP && neighborLocation < beginR)
                numNeighborsInP++;
            k++;
        }
        if (numNeighborsInP == sizeofP) 
            return true;
        ++j;
    }
    return false;
}

void BKrcdMTP::moveToR( int vertex, 
                     int* vertexSets, int* vertexLookup, 
                     int **neighborsInP, int *numNeighbors,
                     int* pBeginX, int *pBeginP, int *pBeginR, 
                     int* pNewBeginX, int* pNewBeginP, int *pNewBeginR)
{
        int vertexLocation = vertexLookup[vertex];

        (*pBeginR)--;
        swap(vertexLookup[vertexSets[*pBeginR]], vertexLookup[vertex]);
        swap(vertexSets[vertexLocation], vertexSets[*pBeginR]);

        *pNewBeginX = *pBeginP;
        *pNewBeginP = *pBeginP;
        *pNewBeginR = *pBeginP;

        int sizeOfP = *pBeginR - *pBeginP;
        int j = *pBeginX;

        j = (*pBeginX);
        while(j < *pNewBeginX)
        {
            int neighbor = vertexSets[j];
            int neighborLocation = j;

            int incrementJ = 1;
            int numPotentialNeighbors = numNeighbors[neighbor];

            int k = 0;
            while(k<numPotentialNeighbors)
            {
                if(neighborsInP[neighbor][k] == vertex)
                {
                    (*pNewBeginX)--;
                    swap(vertexLookup[neighbor], vertexLookup[vertexSets[*pNewBeginX]]);
                    swap(vertexSets[neighborLocation], vertexSets[*pNewBeginX]);
                    incrementJ=0;
                    break;
                }
                k++;
            }
            if(incrementJ) j++;
        }

        j = (*pBeginP);
        while(j<(*pBeginR))
        {
            int neighbor = vertexSets[j];
            int neighborLocation = j;

            //int numPotentialNeighbors = min(sizeOfP+1, numNeighbors[neighbor]);
            int numPotentialNeighbors = numNeighbors[neighbor];

            int k = 0;
            while(k<numPotentialNeighbors)
            {
                if(neighborsInP[neighbor][k] == vertex)
                {
                    swap(vertexLookup[neighbor], vertexLookup[vertexSets[*pNewBeginR]]);
                    swap(vertexSets[neighborLocation], vertexSets[*pNewBeginR]);
                    (*pNewBeginR)++;
                    break;
                }
                k++;
            }
            j++;
        }
}

void BKrcdMTP::moveFromRToX( int vertex, 
                          int* vertexSets, int* vertexLookup,               
                          int**  neighborsInP, int* numNeighbors,
                          int* pBeginX, int* pBeginP, int* pBeginR )
{
    int vertexLocation = vertexLookup[vertex];

    swap(vertexLookup[vertexSets[*pBeginP]], vertexLookup[vertex]);
    swap(vertexSets[vertexLocation], vertexSets[*pBeginP]);

    ++*pBeginP;
    ++*pBeginR;
}

void BKrcdMTP::BKReverseRecursive(long* cliqueCount,
                                list<int> &partialClique, 
                                int* vertexSets, int* vertexLookup,
                                int** neighborsInP, int* numNeighbors,
                                int beginX, int beginP, int beginR)
{
    // numberOfRecurrciveCall++;
    if(beginX >= beginP && beginP >= beginR)
    {
        (*cliqueCount)++;
        return;
    }

    if(beginP >= beginR)
        return;

    // add candiate vertices to the partial clique one at a time and 
    // search for maximal cliques
    list<int> pivots;
    while(true)
    {
        int pivot = 0;
        if (findMostNonNeighbors ( vertexSets, vertexLookup, neighborsInP, numNeighbors, beginX, beginP, beginR, &pivot))
        {
            if (!existCommonNeighborOfPinX(vertexSets, vertexLookup, neighborsInP, numNeighbors, beginX, beginP, beginR)) 
            {
                //R U P is a maximal clique now
                for (int j = beginP; j < beginR; ++j)
                    partialClique.push_back(vertexSets[j]);
                (*cliqueCount)++;

                for (int j = beginP; j < beginR; ++j)
                    partialClique.pop_back();
            }

            // swap vertices that were moved to X back into P, for higher recursive calls.
            for (auto p : pivots)
            {
                int vertexLocation = vertexLookup[p];
                swap(vertexSets[vertexLocation], vertexSets[--beginP]);
                swap(vertexLookup[p], vertexLookup[vertexSets[vertexLocation]]);
            }
            return ;
        }
        pivots.push_back(pivot);

        int newBeginX, newBeginP, newBeginR;

        // add vertex into partialClique, representing R.
        partialClique.push_back(pivot);

        // swaps pivot into R and update all data structures 
        moveToR( pivot, 
                 vertexSets, vertexLookup, 
                 neighborsInP, numNeighbors,
                 &beginX, &beginP, &beginR, 
                 &newBeginX, &newBeginP, &newBeginR);

        // recursively compute maximal cliques with new sets R, P and X
        BKReverseRecursive(cliqueCount, 
                           partialClique, 
                           vertexSets, vertexLookup,
                           neighborsInP, numNeighbors,
                           newBeginX, newBeginP, newBeginR);

        // remove vertex from partialClique
        partialClique.pop_back();

        moveFromRToX( pivot, 
                      vertexSets, vertexLookup,
                      neighborsInP, numNeighbors,
                      &beginX, &beginP, &beginR );
        int sizeOfP = beginR-beginP;
    }

}

void BKrcdMTP::moveToRDgcy( int vertex, 
                         int* vertexSets, int* vertexLookup, 
                         int** neighborsInP, int* numNeighbors,
                         int* pBeginX, int *pBeginP, int *pBeginR, 
                         int* pNewBeginX, int* pNewBeginP, int *pNewBeginR,
                         int& newPivot, int& newBoundSize, int boundSize, int* degreeInPArray)
{
    int vertexLocation = vertexLookup[vertex];

    (*pBeginR)--;
    vertexSets[vertexLocation] = vertexSets[*pBeginR];
    vertexLookup[vertexSets[*pBeginR]] = vertexLocation;
    vertexSets[*pBeginR] = vertex;
    vertexLookup[vertex] = *pBeginR;

    *pNewBeginX = *pBeginP;
    *pNewBeginP = *pBeginP;
    *pNewBeginR = *pBeginP;

    int sizeOfP = *pBeginR - *pBeginP;

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

    j = (*pNewBeginX);

    while(j < *pNewBeginR)
    {
        int thisVertex = vertexSets[j];

        int numPotentialNeighbors = min(boundSize, numNeighbors[thisVertex]);

        int numNeighborsInP = 0;

        int k = 0;
        while(k < numPotentialNeighbors)
        {
            int neighbor = neighborsInP[thisVertex][k];
            int neighborLocation = vertexLookup[neighbor];
            if(neighborLocation >= *pNewBeginP && neighborLocation < *pNewBeginR)
            {
                neighborsInP[thisVertex][k] = neighborsInP[thisVertex][numNeighborsInP];
                neighborsInP[thisVertex][numNeighborsInP] = neighbor;
                numNeighborsInP++;
            }
            k++;
        }
        degreeInPArray[thisVertex] = numNeighborsInP;
        if(numNeighborsInP>newBoundSize)
        {
            newBoundSize = numNeighborsInP;
            newPivot = thisVertex;
        }
        j++;
    }
}

void BKrcdMTP::findBestPivotNonNeighborsDgcy( int** pivotNonNeighbors, int* numNonNeighbors,
                                                int* vertexSets, int* vertexLookup,
                                                int** neighborsInP, int* numNeighbors,
                                                int beginX, int beginP, int beginR, int pivot, int boundSize, 
                                                int currentDepth, int* degreeInPArray)
{
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

        int numPivotNeighbors = min(beginR - beginP, numNeighbors[pivot]);

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

void BKrcdMTP::moveFromRToXDgcy( int vertex, 
                                    int* vertexSets, int* vertexLookup, 
                                    int* pBeginX, int* pBeginP, int* pBeginR )
{
    int vertexLocation = vertexLookup[vertex];

    //swap vertex into X and increment beginP and beginR
    vertexSets[vertexLocation] = vertexSets[*pBeginP];
    vertexLookup[vertexSets[*pBeginP]] = vertexLocation;
    vertexSets[*pBeginP] = vertex;
    vertexLookup[vertex] = *pBeginP;

    *pBeginP = *pBeginP + 1;
    *pBeginR = *pBeginR + 1;

}

void BKrcdMTP::listAllMaximalCliquesDgcyRecursive(long* cliqueCount,
                                               list<int> &partialClique, 
                                               int* vertexSets, int* vertexLookup,
                                               int** neighborsInP, int* numNeighbors,
                                               int beginX, int beginP, int beginR, 
                                               int pivot, int boundSize, int currentDepth, int* degreeInPArray)
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

    bool findPivot = true;
    if (numCandidatesToIterateThrough<=2)
    {
        findPivot = false;
    }

    if(findPivot)
    {
    findBestPivotNonNeighborsDgcy( &myCandidatesToIterateThrough,
                                         &numCandidatesToIterateThrough,
                                         vertexSets, vertexLookup,
                                         neighborsInP, numNeighbors,
                                         beginX, beginP, beginR, pivot, boundSize, currentDepth, degreeInPArray);
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
        int vertex = myCandidatesToIterateThrough[iterator];

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("%d ", vertex);
        #endif

        int newBeginX, newBeginP, newBeginR;

        partialClique.push_back(vertex);

        int newPivot = -1;
        int newBoundSize = -1;
        moveToRDgcy( vertex, 
                     vertexSets, vertexLookup, 
                     neighborsInP, numNeighbors,
                     &beginX, &beginP, &beginR, 
                     &newBeginX, &newBeginP, &newBeginR, 
                     newPivot, newBoundSize, boundSize, degreeInPArray);

        listAllMaximalCliquesDgcyRecursive(cliqueCount,
                                           partialClique, 
                                           vertexSets, vertexLookup,
                                           neighborsInP, numNeighbors,
                                           newBeginX, newBeginP, newBeginR,
                                           newPivot, newBoundSize,currentDepth+1, degreeInPArray);

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("b ");
        #endif

        partialClique.pop_back();

        moveFromRToXDgcy( vertex, 
                                vertexSets, vertexLookup,
                                &beginX, &beginP, &beginR );

        iterator++;
    }

    iterator = 0;

    while(iterator < numCandidatesToIterateThrough)
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
    }
    free(myCandidatesToIterateThrough);
}

long BKrcdMTP::run(){
    clock_t start, end;
    long cliqueCount = 0;
    list<int> partialClique;
    int* candidates;
    int numCandidates;

    start  = clock();
    int beginX = 0;
    int beginP = 0;
    int beginR = n;
    int currentDepth = 1;
    long scount = 0, kcount = 0, noise = 0, dive = 0;
    orderingArray = computeDegeneracyOrderArray();

    for(int i = 0; i<n ;i++)
    {
        int vertex = (int)orderingArray[i]->vertex;

        // add vertex to partial clique R
        int newBeginX, newBeginP, newBeginR;
        int newBeginX_, newBeginP_, newBeginR_;

        int pivot = -1;
        int boundSize = -1;
        fillInPandX( i, vertex, 
                        vertexSets, vertexLookup, 
                        orderingArray,
                        neighborsInP, numNeighbors,
                        &beginX, &beginP, &beginR, 
                        &newBeginX, &newBeginP, &newBeginR,
                        scount, kcount, noise, dive, pivot, boundSize, degreeInPArray);

        for (int j = newBeginR; j < n; ++j)
            partialClique.push_back(vertexSets[j]);

        long long sizeofP = scount+kcount;

        if (sizeofP == 0) {
            if (newBeginX == newBeginP) {
                cliqueCount++;
            }
            beginR = beginR + 1;
            for (int j = newBeginR; j < n; ++j)
                partialClique.pop_back();
            continue;
        }

        double benefit = 0.0;
        if (dive == 0) {
            benefit = scount+4.5-2.8*kcount;
        } else if (dive == 1) {
            benefit = scount+8.0-2.8*kcount;
        } else if (dive >= 2) {
            benefit = scount+11.0-2.8*kcount;
        }        

        if (benefit >= 0.0) {
            // BKreverseCount++;
            BKReverseRecursive(&cliqueCount, partialClique,
                                vertexSets, vertexLookup,
                                neighborsInP, numNeighbors,
                                newBeginX, newBeginP, newBeginR); 
        } else {
            listAllMaximalCliquesDgcyRecursive( &cliqueCount,
                                              partialClique, 
                                              vertexSets, vertexLookup,
                                              neighborsInP, numNeighbors,
                                              newBeginX, newBeginP, newBeginR,
                                              pivot, boundSize, currentDepth, degreeInPArray);
            // BKpivotCount++;
        }
        
        beginR = beginR + 1;
        for (int j = newBeginR; j < n; ++j)
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
    BKrcdMTP alg(filename);
    alg.run();
}