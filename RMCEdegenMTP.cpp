// local includes
#include "RMCEdegenMTP.h"
#include <queue>

using namespace std;

// static unsigned long numberOfRecurrciveCall=0;

RMCEdegenMTP::RMCEdegenMTP(string const &fileName){
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
    vector<vector<int>> adjList(n);
    // vector<list<int>> adjList(n);
    int u, v; // endvertices, to read edges.
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
        adjList[u].emplace_back(v);
        i++;
    }

    int e_idx = 0;
    inds = (int*)malloc(sizeof(int) * (n+1));
    vals = (int*)malloc(sizeof(int) * m);
    isEdgeDeleted = (int*)malloc(sizeof(int) * m);
    for (int i = 0; i < n; ++i) {
        inds[i] = e_idx;
        sort(adjList[i].begin(), adjList[i].end());
        for (auto u : adjList[i]) {
            isEdgeDeleted[e_idx] = -1;
            vals[e_idx++] = u;
        }
    }

    inds[n] = m;
    vertexArray = (int*)malloc(n * sizeof(int));
    vertexPos = (int*)malloc(n * sizeof(int));
    neighborsInP = (int**)malloc(n * sizeof(int*));
    // neighborsInP.resize(n);
    numNeighbors = (int*)malloc(n * sizeof(int));
    tmpNumNeighbors = (int*)malloc(n * sizeof(int));
    bfsParent = (bool*)malloc(sizeof(bool) * n);
    orders = (int*)malloc(sizeof(int) * (n));
    orderedVals = (int*)malloc(sizeof(int) * m);
    seps = (int*)malloc(sizeof(int) * (n));
    seps2 = (int*)malloc(sizeof(int) * (n));
    o2v = (int*)malloc(sizeof(int) * (n));
    // RParent = (bool*)malloc(sizeof(bool) * n);
    Xpruned = (int*)malloc(sizeof(int) * n);
    vertexDegree = (int*)malloc(sizeof(int) * n);
    
    i = 0;
    while(i<n)
    {
        bfsParent[i] = false;
        Xpruned[i] = n;
        // pruned[i] = 0;
        vertexPos[i] = i;
        vertexArray[i] = i;
        neighborsInP[i] = (int*)malloc(1 * sizeof(int));
        numNeighbors[i] = 1;
        vertexDegree[i] = 0;
        i++;
    }
}

RMCEdegenMTP::~RMCEdegenMTP() {
    free(inds);
    free(seps);
    free(seps2);
    free(vals);
    free(orderedVals);
    free(orders);
    free(o2v);
    free(vertexArray);
    free(vertexPos);
    free(Xpruned);
    free(bfsParent);
    free(vertexDegree);
    // free(pruned);
    for(int i=0; i<n; i++)
    {
        free(neighborsInP[i]);
    }
    free(numNeighbors);
}

void RMCEdegenMTP::computeDegeneracyOrder(){
    vector<list<int>> verticesByDegree(n);
    vector<list<int>::iterator> vertexLocator(n);
    vector<int> degree(n);

    for(int i=0; i<n; i++)
    {
        degree[i] = inds[i+1] - inds[i];
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

            orders[vertex] = numVerticesRemoved;
            o2v[numVerticesRemoved] = vertex;
            degree[vertex] = -1;

            int l = inds[vertex], r = inds[vertex+1]-1;

            for(int j=inds[vertex]; j < inds[vertex+1]; j++)
            {
                int neighbor = vals[j];
                if(degree[neighbor]!=-1)
                {
                    verticesByDegree[degree[neighbor]].erase(vertexLocator[neighbor]);
                    orderedVals[r--] = neighbor;

                    degree[neighbor]--;

                    if(degree[neighbor] != -1)
                    {
                        verticesByDegree[degree[neighbor]].push_front(neighbor);
                        vertexLocator[neighbor] = verticesByDegree[degree[neighbor]].begin();
                    }
                }
                else
                {
                    orderedVals[l++] = neighbor;
                }
            }
            seps[vertex] = l;
            numVerticesRemoved++;
            currentDegree = max(currentDegree-1, 0);
        }
        else
        {
            currentDegree++;
        }
    }
}

long RMCEdegenMTP::computeDegeneracyWithGlobalReduction(int* offset){
    long currentCliques = 0;
    vector<list<int>> verticesByDegree(n);
    vector<list<int>::iterator> vertexLocator(n);
    vector<int> degree(n);

    for(int i=0; i<n; i++)
    {
        degree[i] = inds[i+1] - inds[i];
        verticesByDegree[degree[i]].push_front(i);
        vertexLocator[i] = verticesByDegree[degree[i]].begin();
    }
    int numVerticesRemoved = 0;

    // non triangle edge removal
    for (int i = 0; i < n; i++)
    {
        if (degree[i] <= 0)
        {
            continue;
        }
        for (int j = inds[i]; j < inds[i+1]; j++)
        {
            int k = vals[j];
            if (i>k || degree[k] <= 0 || isEdgeDeleted[j] != -1)
            {
                continue;
            }
            int u_ind = inds[i], v_ind = inds[k];
            int reverse_j = -1;
            bool canDelete = true;
            while (u_ind < inds[i+1] && v_ind < inds[k+1]) {
                if (vals[v_ind] == i)
                {
                    reverse_j = v_ind;
                }
                if(degree[vals[u_ind]]<=0 || isEdgeDeleted[u_ind] == 1){
                    u_ind++;
                    continue;
                } 
                if(degree[vals[v_ind]]<=0 || isEdgeDeleted[v_ind] == 1){
                    v_ind++;
                    continue;
                } 
                if (vals[u_ind] == vals[v_ind]) {  
                    canDelete = false;  
                    if (isEdgeDeleted[u_ind] == -1 && isEdgeDeleted[v_ind] == -1)
                    {
                        int reverse1 = lower_bound(vals+inds[vals[u_ind]], vals+inds[vals[u_ind]+1], i) - vals;
                        int reverse2 = lower_bound(vals+reverse1, vals+inds[vals[u_ind]+1], k) - vals;
                        isEdgeDeleted[u_ind] = 0;
                        isEdgeDeleted[v_ind] = 0;
                        isEdgeDeleted[reverse1] = 0;
                        isEdgeDeleted[reverse2] = 0;
                    }
                    else if (isEdgeDeleted[u_ind] == -1){
                        isEdgeDeleted[u_ind] = 0;
                        isEdgeDeleted[lower_bound(vals+inds[vals[u_ind]], vals+inds[vals[u_ind]+1], k) - vals] = 0;
                    }
                    else if (isEdgeDeleted[v_ind] == -1){
                        isEdgeDeleted[v_ind] = 0;
                        isEdgeDeleted[lower_bound(vals+inds[vals[v_ind]], vals+inds[vals[v_ind]+1], i) - vals] = 0;
                    }
                    v_ind++;
                    u_ind++;
                    break;
                }
                else if (vals[u_ind] < vals[v_ind]) u_ind++;
                else v_ind++;
            }
            if (reverse_j==-1)
            {
                reverse_j = lower_bound(vals+v_ind, vals+inds[k+1], i) - vals; 
            }
            if (canDelete)
            {
                currentCliques++;
                isEdgeDeleted[j]=1;
                isEdgeDeleted[reverse_j]=1;
                verticesByDegree[degree[i]].erase(vertexLocator[i]);
                degree[i]--;
                verticesByDegree[degree[k]].erase(vertexLocator[k]);
                degree[k]--;
                verticesByDegree[degree[i]].push_front(i);
                vertexLocator[i] = verticesByDegree[degree[i]].begin();
                verticesByDegree[degree[k]].push_front(k);
                vertexLocator[k] = verticesByDegree[degree[k]].begin();
            }
            else{
                isEdgeDeleted[j] = 0;
                isEdgeDeleted[reverse_j] = 0;
            }
        }
    }

    // vertex reduction
    while (!verticesByDegree[1].empty() || !verticesByDegree[2].empty() || !verticesByDegree[0].empty()){

        while (!verticesByDegree[2].empty())
        {
            int const vertex = verticesByDegree[2].front();
            verticesByDegree[2].pop_front();
            int neighbor1=-1, neighbor2=-1;
            for(int j=inds[vertex]; j < inds[vertex+1]; j++)
            {
                int neighbor = vals[j];
                if(degree[neighbor] > 0 && isEdgeDeleted[j]!=1)
                {
                    if (neighbor1 == -1)
                    {
                        neighbor1 = neighbor;
                    }
                    else{
                        neighbor2 = neighbor;
                        break;
                    }
                }
            }
            int* p = lower_bound(vals+inds[neighbor1], vals+inds[neighbor1+1], neighbor2);
            bool found = (p != vals+inds[neighbor1+1] && !(neighbor2<*p));
            if (!found)
            {
                verticesByDegree[degree[neighbor1]].erase(vertexLocator[neighbor1]);
                degree[neighbor1]--;
                verticesByDegree[degree[neighbor2]].erase(vertexLocator[neighbor2]);
                degree[neighbor2]--;
                currentCliques = currentCliques + 2;
                verticesByDegree[degree[neighbor1]].push_front(neighbor1);
                vertexLocator[neighbor1] = verticesByDegree[degree[neighbor1]].begin();
                verticesByDegree[degree[neighbor2]].push_front(neighbor2);
                vertexLocator[neighbor2] = verticesByDegree[degree[neighbor2]].begin();
                verticesByDegree[2].erase(vertexLocator[vertex]);
                orders[vertex] = numVerticesRemoved;
                o2v[numVerticesRemoved] = -1;
                degree[vertex] = -2;
                numVerticesRemoved++;
                continue;
            }
            else{ // remove all deg 2 vertices
                verticesByDegree[degree[neighbor1]].erase(vertexLocator[neighbor1]);
                verticesByDegree[degree[neighbor2]].erase(vertexLocator[neighbor2]);
                int u_ind = inds[neighbor1], v_ind = inds[neighbor2];
                bool canDelete = true;
                while (u_ind < inds[neighbor1+1] && v_ind < inds[neighbor2+1]) {
                    if(degree[vals[u_ind]]<=0 || vals[u_ind]==vertex || isEdgeDeleted[u_ind] == 1){
                        u_ind++;
                        continue;
                    } 
                    if(degree[vals[v_ind]]<=0 || vals[v_ind]==vertex || isEdgeDeleted[v_ind] == 1){
                        v_ind++;
                        continue;
                    } 
                    if (vals[u_ind] == vals[v_ind]) {
                        canDelete = false;
                        break;
                    }
                    else if (vals[u_ind] < vals[v_ind]) u_ind++;
                    else v_ind++;
                }
                if(canDelete){
                    isEdgeDeleted[p-vals] = 1;
                    isEdgeDeleted[lower_bound(vals+inds[neighbor2], vals+inds[neighbor2+1], neighbor1) -vals] = 1;
                    degree[neighbor1]--;
                    degree[neighbor2]--;
                }
                degree[neighbor1]--;
                degree[neighbor2]--;
                currentCliques++;
                verticesByDegree[degree[neighbor1]].push_front(neighbor1);
                vertexLocator[neighbor1] = verticesByDegree[degree[neighbor1]].begin();
                verticesByDegree[degree[neighbor2]].push_front(neighbor2);
                vertexLocator[neighbor2] = verticesByDegree[degree[neighbor2]].begin();
                orders[vertex] = numVerticesRemoved;
                o2v[numVerticesRemoved] = -1;
                degree[vertex] = -2;
                numVerticesRemoved++;
                continue;
            }
        }

        while (!verticesByDegree[1].empty())
        {
            currentCliques++;
            int const vertex = verticesByDegree[1].front();
            verticesByDegree[1].pop_front();
            orders[vertex] = numVerticesRemoved;
            o2v[numVerticesRemoved] = -1;
            degree[vertex] = -2;
            for(int j=inds[vertex]; j < inds[vertex+1]; j++)
            {
                int neighbor = vals[j];
                if(degree[neighbor]!=-2 && isEdgeDeleted[j]!=1)
                {
                    verticesByDegree[degree[neighbor]].erase(vertexLocator[neighbor]);
                    degree[neighbor]--;
                    if(degree[neighbor] != -1)
                    {
                        verticesByDegree[degree[neighbor]].push_front(neighbor);
                        vertexLocator[neighbor] = verticesByDegree[degree[neighbor]].begin();
                    }
                }
            }
            numVerticesRemoved++;
        }

        // degree zero removal
        while (!verticesByDegree[0].empty())
        {
            int const vertex = verticesByDegree[0].front();
            verticesByDegree[0].pop_front();
            orders[vertex] = numVerticesRemoved;
            o2v[numVerticesRemoved] = -1;
            degree[vertex] = -2;
            numVerticesRemoved++;
        }
    }

    *offset = numVerticesRemoved;

    int currentDegree = 0;
    int type1 = 0;
    while(numVerticesRemoved < n)
    {
        if(!verticesByDegree[currentDegree].empty())
        {
            degeneracy = max(degeneracy,currentDegree);
            int const vertex = verticesByDegree[currentDegree].front();
            verticesByDegree[currentDegree].pop_front();
            orders[vertex] = numVerticesRemoved;
            o2v[numVerticesRemoved] = vertex;
            degree[vertex] = -1;
            int l = inds[vertex], r = inds[vertex+1]-1;
            for(int j=inds[vertex]; j < inds[vertex+1]; j++)
            {
                if(isEdgeDeleted[j] == 1){
                    // orderedVals[l++] = -1;
                    continue;
                }
                int neighbor = vals[j];
                if (degree[neighbor] == -2)
                {
                    continue; // mark deleted vertex as -1, which can not be accessed later
                }
                else if(degree[neighbor] ==-1)
                {
                    orderedVals[l++] = neighbor;
                }
                else
                {
                    type1++;
                    verticesByDegree[degree[neighbor]].erase(vertexLocator[neighbor]);
                    orderedVals[r--] = neighbor;
                    degree[neighbor]--;
                    if(degree[neighbor] != -1)
                    {
                        verticesByDegree[degree[neighbor]].push_front(neighbor);
                        vertexLocator[neighbor] = verticesByDegree[degree[neighbor]].begin();
                    }
                }
            }
            seps[vertex] = l;
            seps2[vertex] = r+1;
            numVerticesRemoved++;
            currentDegree = max(currentDegree-1, 0);
        }
        else
        {
            currentDegree++;
        }
    }
    return currentCliques;
}

void RMCEdegenMTP::setCandidates(int** candidates, int* numCandidates,
                                   int posX, int posP, int posR, int pivot, int boundSize, int* vertexDegree,
                                   long* cliqueCount, int currentDepth){
    // mtp
    if(currentDepth == 1)
    {
        *candidates = (int*)malloc((posR - posP) * sizeof(int));
        memcpy(*candidates, &vertexArray[posP], (posR-posP)*sizeof(int));
        *numCandidates = 0;

        std::queue<int> coverQueue;
        std::queue<int> candidateQueue;
        std::queue<int> pivotQueue;
        // int *vertexFlagSet = (int*)calloc((posR - posX), sizeof(int));
        vector<int> vertexFlagSet(posR - posX,0);

        int initialPivot = pivot;
        int initialpivotIndex = vertexPos[pivot];

        int unprocessedVertices = posR - posP;
        int threshold = 4;

        pivotQueue.push(initialPivot);

        vertexFlagSet[initialpivotIndex-posX] = 3; // 3 means pivot
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
                int neighborLocation = vertexPos[neighbor];

                if(neighborLocation >= posP && neighborLocation < posR)
                {
                    if(vertexFlagSet[neighborLocation-posX] == 0)
                    {
                        coverQueue.push(neighbor);
                        vertexFlagSet[neighborLocation-posX] = 1; // 1 means covered
                        unprocessedVertices--;
                        (*candidates)[neighborLocation-posP] = -1;
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
                        int candidateLocation = vertexPos[candidate];

                        if(candidateLocation >= posP && candidateLocation < posR)
                        {
                            int flag = vertexFlagSet[candidateLocation-posX];
                            if(flag == 0)
                            {
                                candidateQueue.push(candidate);
                                vertexFlagSet[candidateLocation-posX] = 2; // 2 means pillar vertex
                                vertexDegree[candidate]--;
                                unprocessedVertices--;
                            }
                            else if(flag == 2)
                            {
                                vertexDegree[candidate]--;
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

                int neighborCountLCS = vertexDegree[vertex];
                int neighborCountLPS = 0;

                int n = 0;
                while(n<numPotentialNeighbors)
                {
                    int neighbor = neighborsInP[vertex][n];
                    int neighborLocation = vertexPos[neighbor];

                    if(neighborLocation >= posP && neighborLocation < posR)
                    {
                        if(vertexFlagSet[neighborLocation-posX] == 0)
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
                        int candidateLocation = vertexPos[candidate];

                        if(candidateLocation >= posP && candidateLocation < posR)
                        {
                            int flag = vertexFlagSet[candidateLocation-posX];
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

            if(candidatePivotFound)
            {
                pivotQueue.push(candidatePivot);
                int loc = vertexPos[candidatePivot];
                if(loc>=posX && loc < posP)
                {
                    vertexFlagSet[loc-posX] = 3;
                }
                else if(loc>=posP && loc < posR)
                {
                    if (vertexFlagSet[loc-posX] == 0)
                    {
                        vertexFlagSet[loc-posX] = 3;
                        unprocessedVertices--;
                    }
                    else
                    {
                        vertexFlagSet[loc-posX] = 3;
                    }
                }
            }
            else
                break;

        }

        int tmpV = (*candidates)[0];
        int j = 1;
        while(j<posR-posP)
        {
            int vertex = (*candidates)[j];
            if(vertex < 0)
            {
                j++;
                continue;
            }
            (*candidates)[*numCandidates] = vertex;
            (*numCandidates)++;
            j++;
        }
        if(tmpV < 0)
        {
            return;
        }
        (*candidates)[*numCandidates] = tmpV;
        (*numCandidates)++;

        // for(int i = posP; i < posR; i++)
        // {
        //     int flag = vertexFlagSet[i-posX];
        //     if(flag == 3 ||flag == 2 || flag == 0)
        //     {
        //         (*candidates)[*numCandidates] = vertexArray[i];
        //         (*numCandidates)++;
        //     }
        // }
        // sum += (*numNonNeighbors);

    }
    else
    {
        // *subSize = 0;
        *candidates = (int*)malloc((posR-posP)*sizeof(int));
        memcpy(*candidates, &vertexArray[posP], (posR-posP)*sizeof(int));
        // *candidates = (int*)calloc(posR-posP,sizeof(int));
        *numCandidates = 0;
        
        int numPivotNeighbors = min(boundSize, numNeighbors[pivot]);
        int j = 0;
        while(j<numPivotNeighbors)
        {
            int neighbor = neighborsInP[pivot][j];
            int neighPos = vertexPos[neighbor];
            if(neighPos >= posP && neighPos < posR)
            {
                (*candidates)[neighPos-posP] = -1;
            }
            else
            {
                break;
            }
            j++;
        }

        int tmpV = (*candidates)[0];
        j = 1;
        while(j<posR-posP)
        {
            int vertex = (*candidates)[j];
            if(vertex < 0)
            {
                j++;
                continue;
            }
            (*candidates)[*numCandidates] = vertex;
            (*numCandidates)++;
            j++;
        }
        if(tmpV < 0)
        {
            return;
        }
        (*candidates)[*numCandidates] = tmpV;
        (*numCandidates)++;

        // for(int i = posP; i < posR; i++)
        // {
        //     if((*candidates)[i-posP] == -1)
        //     {
        //         continue;
        //     }
        //     (*candidates)[*numCandidates] = vertexArray[i];
        //     (*numCandidates)++;
        // }

    }

}

void RMCEdegenMTP::P2RandFindPivot(int vertex, 
                                     int* posX, int* posP, int* posR, 
                                     int* newPosX, int* newPosP, int* newPosR,
                                     int* pivot, int *pivotIntersectionSize, int boundSize, long* cliqueCount,
                                     int currentDepth){
    int vPos = vertexPos[vertex];

    (*posR)--;
    vertexArray[vPos] = vertexArray[*posR];
    vertexPos[vertexArray[*posR]] = vPos;
    vertexArray[*posR] = vertex;
    vertexPos[vertex] = *posR;

    *newPosX = *posP;
    *newPosP = *posP;
    int sizeOfP = *posR - *posP;
    int j = *posX;
    while(j<*newPosX)
    {
        int neighbor = vertexArray[j];
        int neiPos = j;
        int numPotentialNeighbors = min(boundSize, numNeighbors[neighbor]);
        bool find = false;

        int k = 0;
        while(k<numPotentialNeighbors)
        {
            if(neighborsInP[neighbor][k] == vertex)
            {
                (*newPosX)--;
                vertexArray[neiPos] = vertexArray[(*newPosX)];
                vertexPos[vertexArray[(*newPosX)]] = neiPos;
                vertexArray[(*newPosX)] = neighbor;
                vertexPos[neighbor] = (*newPosX);
                find=true;
                break;
            }
            k++;
        }
        if(!find) j++;
    }

    *newPosR = *posP;
    int numPotentialNeighbors = min(boundSize, numNeighbors[vertex]);
    j = 0;
    while(j <numPotentialNeighbors)
    {
        int neighbor = neighborsInP[vertex][j];
        int neighborLocation = vertexPos[neighbor];
        if(neighborLocation >= *posP && neighborLocation < *posR)
        {
            bfsParent[neighbor] = false;
            vertexArray[neighborLocation] = vertexArray[(*newPosR)];
            vertexPos[vertexArray[(*newPosR)]] = neighborLocation;
            vertexArray[(*newPosR)] = neighbor;
            vertexPos[neighbor] = (*newPosR);
            (*newPosR)++;
        }
        j++;
    }

    *pivot = -1;
    *pivotIntersectionSize = -1;
    for (j = (*newPosX); j < (*newPosP); j++)
    {
        int v = vertexArray[j];
        int numPotentialNeighbors = min(boundSize, numNeighbors[v]);
        tmpNumNeighbors[v] = 0;
        int k = 0;
        while(k < numPotentialNeighbors)
        {
            int neighbor = neighborsInP[v][k];
            int neiPos = vertexPos[neighbor];
            if(neiPos >= *newPosP && neiPos < *newPosR)
            {
                bfsParent[neighbor] = true;
                neighborsInP[v][k] = neighborsInP[v][tmpNumNeighbors[v]];
                neighborsInP[v][tmpNumNeighbors[v]] = neighbor;
                tmpNumNeighbors[v]++;
            }
            k++;
        }
        if (tmpNumNeighbors[v] > *pivotIntersectionSize)
        {
            *pivotIntersectionSize = tmpNumNeighbors[v];
            *pivot = v;
            if (tmpNumNeighbors[v] == (*newPosR - *newPosP))
            {
                *pivot = -2;
                return;
            }
        }
    }

    int numEdges=0; // earlystopping
    j = (*newPosP);
    int degZeroUpdated = 0;
    while (j < *newPosR)
    {
        int v = vertexArray[j];
        int numPotentialNeighbors = min(boundSize, numNeighbors[v]);
        tmpNumNeighbors[v] = 0;
        int k = 0;
        while(k < numPotentialNeighbors)
        {
            int neighbor = neighborsInP[v][k];
            int neiPos = vertexPos[neighbor];
            if(neiPos >= *newPosP && neiPos < *newPosR)
            {
                neighborsInP[v][k] = neighborsInP[v][tmpNumNeighbors[v]];
                neighborsInP[v][tmpNumNeighbors[v]] = neighbor;
                tmpNumNeighbors[v]++;
                numEdges++;
            }
            k++;
        }

        if (tmpNumNeighbors[v] == 0)
        {
            degZeroUpdated++;
            if (!bfsParent[v])
            {
                (*cliqueCount)++;
            }
            (*newPosR)--;
            vertexArray[j] = vertexArray[*newPosR];
            vertexPos[vertexArray[*newPosR]] = j;
            vertexArray[*newPosR] = v;
            vertexPos[v] = *newPosR;
            continue;
        }
        j++;
    }
    int maxEcount = (*newPosR - *newPosP) * (*newPosR - *newPosP - 1);
    if (numEdges == maxEcount && degZeroUpdated == 0)
    {
        *pivot = -3; // already maximal
        return;
    }
    if (numEdges == 0 && degZeroUpdated > 0){
        *pivot = -2; 
        return;
    }

    j = (*newPosP);
    while (j < *newPosR){
        int v = vertexArray[j];
        if(tmpNumNeighbors[v] == 1)
        {
            if (!bfsParent[v]){ //can return a clique directly
                degZeroUpdated++;
                (*cliqueCount)++;
                (*newPosR)--;
                vertexArray[j] = vertexArray[(*newPosR)];
                vertexPos[vertexArray[(*newPosR)]] = j;
                vertexArray[(*newPosR)] = v;
                vertexPos[v] = (*newPosR);
                int v2 = neighborsInP[v][0];
                if (tmpNumNeighbors[v2] == 1)
                {
                    int j2 = vertexPos[v2];
                    (*newPosR)--;
                    vertexArray[j2] = vertexArray[(*newPosR)];
                    vertexPos[vertexArray[(*newPosR)]] = j2;
                    vertexArray[(*newPosR)] = v2;
                    vertexPos[v2] = (*newPosR);
                    if (j2 < j && j < (*newPosR))
                    {
                        j--;
                        int tmp = vertexArray[j2];
                        vertexArray[j2] = vertexArray[j];
                        vertexPos[vertexArray[j]] = j2;
                        vertexArray[j] = tmp;
                        vertexPos[tmp] = j;
                    }
                }
                else {
                    for (int tmpIter = 0; tmpIter < tmpNumNeighbors[v2]; tmpIter++)
                    {
                        if (neighborsInP[v2][tmpIter] == v)
                        {
                            tmpNumNeighbors[v2]--;
                            neighborsInP[v2][tmpIter] = neighborsInP[v2][tmpNumNeighbors[v2]];
                            neighborsInP[v2][tmpNumNeighbors[v2]] = v;
                            break;
                        }
                    }  
                }
                continue;
            }
        }
        j++;
    }
    if (*newPosR-*newPosP <= 0){
        *pivot = -2;
        return;
    }
    
    j = (*newPosP);
    int newNewPosR = *newPosR;
    while(j < newNewPosR){
        int v = vertexArray[j];
        if (tmpNumNeighbors[v] == (*newPosR - *newPosP - 1))
        {
            newNewPosR--;
            vertexArray[j] = vertexArray[newNewPosR];
            vertexPos[vertexArray[newNewPosR]] = j;
            vertexArray[newNewPosR] = v;
            vertexPos[v] = newNewPosR;
            continue;
        }
        if (tmpNumNeighbors[v] > *pivotIntersectionSize)
        {
            *pivotIntersectionSize = tmpNumNeighbors[v];
            *pivot = v;
        }
        j++;
    }
    if (newNewPosR != *newPosR || degZeroUpdated)
    {
        j = *newPosX;
        *pivot = -1;
        *pivotIntersectionSize = -1;
        while(j < *newPosP){
            int v = vertexArray[j];
            int newNeiCount = 0;
            int rCount = 0;
            if (tmpNumNeighbors[v] < (*newPosR -newNewPosR))
            {
                vertexArray[j] = vertexArray[*newPosX];
                vertexPos[vertexArray[*newPosX]] = j;
                vertexArray[*newPosX] = v;
                vertexPos[v] = *newPosX;
                (*newPosX)++;
                j++;
                continue;
            }
            for (int k = 0; k < tmpNumNeighbors[v]; k++)
            {
                int neighbor = neighborsInP[v][k];
                int neiPos = vertexPos[neighbor];
                if (neiPos >= *newPosP && neiPos < newNewPosR)
                {
                    neighborsInP[v][k] = neighborsInP[v][newNeiCount];
                    neighborsInP[v][newNeiCount] = neighbor;
                    newNeiCount++;
                }
                else if (neiPos >= newNewPosR && neiPos < *newPosR){
                    rCount++;
                }
            }
            if (rCount != (*newPosR - newNewPosR))
            {
                vertexArray[j] = vertexArray[*newPosX];
                vertexPos[vertexArray[*newPosX]] = j;
                vertexArray[*newPosX] = v;
                vertexPos[v] = *newPosX;
                (*newPosX)++;
                j++;
                continue;
            }
            if (newNeiCount > *pivotIntersectionSize)
            {
                *pivot = v;
                *pivotIntersectionSize = newNeiCount;
            }
            if (newNeiCount == (newNewPosR - *newPosP))
            {
                *pivot = -2; // no maximal clique
                return;
            }
            j++;
        }
        j = *newPosP;
        while (j < newNewPosR)
        {
            int v = vertexArray[j];
            int newNeiCount = 0;
            for (int k = 0; k < tmpNumNeighbors[v]; k++)
            {
                int neighbor = neighborsInP[v][k];
                int neiPos = vertexPos[neighbor];
                if (neiPos >= *newPosP && neiPos < newNewPosR)
                {
                    neighborsInP[v][k] = neighborsInP[v][newNeiCount];
                    neighborsInP[v][newNeiCount] = neighbor;
                    newNeiCount++;
                }
                else if (neiPos < *newPosP || neiPos >= *newPosR)
                {
                    break;
                }
            }

            // if(currentDepth == 1)
            // {
            //     vertexDegree[v] = newNeiCount;
            // }

            if (newNeiCount > *pivotIntersectionSize)
            {
                *pivot = v;
                *pivotIntersectionSize = newNeiCount;
            }
            j++;
        }
        *newPosR = newNewPosR;
    }
}

void RMCEdegenMTP::R2X(int vertex, 
                         int* posX, int* posP, int* posR){
    int vPos = vertexPos[vertex];
    vertexArray[vPos] = vertexArray[*posP];
    vertexPos[vertexArray[*posP]] = vPos;
    vertexArray[*posP] = vertex;
    vertexPos[vertex] = *posP;
    *posP = *posP + 1;
    *posR = *posR + 1;
}

long RMCEdegenMTP::run(){
    
    clock_t start, mid, end;
    long cliqueCount = 0;
    list<int> partialClique;
    // int* candidates;
    // int numCandidates;
    int offset;

    start  = clock();

    int posX = 0;
    int posP = 0;
    int posR = n;

    // computeDegeneracyOrder()
    
    cliqueCount = cliqueCount + computeDegeneracyWithGlobalReduction(&offset);
    
    for(int i=offset; i<n ;i++)
    {
        int vertex = o2v[i];
        if (vertex == -1) 
        {
            continue;
        }

        int newPosX, newPosP, newPosR;

        int vPos = vertexPos[vertex];
        posR--;
        vertexArray[vPos] = vertexArray[posR];
        vertexPos[vertexArray[posR]] = vPos;
        vertexArray[posR] = vertex;
        vertexPos[vertex] = posR;

        newPosR = posR;
        newPosP = posR;
        for (int j = inds[vertex+1]-1; j >= seps2[vertex]; j--)
        {
            int neighbor = orderedVals[j];
            int neiPos = vertexPos[neighbor];
            bfsParent[neighbor] = false; //degree zero fliter
            newPosP--;
            vertexArray[neiPos] = vertexArray[newPosP];
            vertexPos[vertexArray[newPosP]] = neiPos;
            vertexArray[newPosP] = neighbor;
            vertexPos[neighbor] = newPosP;

            numNeighbors[neighbor] = 0;
            // neighborsInP[neighbor].resize(min(inds[vertex+1]-seps2[vertex], inds[neighbor+1]+seps[neighbor] -seps2[neighbor] - inds[neighbor])); 
            free(neighborsInP[neighbor]);
            neighborsInP[neighbor]=(int*)malloc(min(inds[vertex+1]-seps2[vertex], inds[neighbor+1]+seps[neighbor] -seps2[neighbor] - inds[neighbor]) * sizeof(int));
        }

        int pivot = -1;
        int maxIntersectionSize = -1;
        newPosX = newPosP;
        for (int j = inds[vertex]; j < seps[vertex]; j++)
        {
            int neighbor = orderedVals[j];
            if (i > Xpruned[neighbor])
            {
                continue;
            }
            int neiPos = vertexPos[neighbor];

            newPosX--;
            vertexArray[neiPos] = vertexArray[newPosX];
            vertexPos[vertexArray[newPosX]] = neiPos;
            vertexArray[newPosX] = neighbor;
            vertexPos[neighbor] = newPosX;

            free(neighborsInP[neighbor]);
            neighborsInP[neighbor] = (int*)malloc(min(newPosR-newPosP,inds[neighbor+1]-seps2[neighbor]) * sizeof(int));
            // neighborsInP[neighbor].resize(min(newPosR-newPosP,inds[neighbor+1]-seps2[neighbor])); 
            numNeighbors[neighbor] = 0;

            for (int k = seps2[neighbor]; k < inds[neighbor+1]; k++)
            {
                int laterNeighbor = orderedVals[k];
                int laterNeiPos = vertexPos[laterNeighbor];
                if(laterNeiPos >= newPosP && laterNeiPos < newPosR)
                {
                    // if (i == 558441) cout << "x:" << neighbor << " " << laterNeighbor << endl;
                    bfsParent[laterNeighbor] = true;
                    neighborsInP[neighbor][numNeighbors[neighbor]] = laterNeighbor;
                    numNeighbors[neighbor]++;
                }
            }
            if (numNeighbors[neighbor] > maxIntersectionSize)
            {
                pivot = neighbor;
                maxIntersectionSize = numNeighbors[neighbor];
            }
        }
        
        int numEdges = 0;
        int j = newPosP;
        while(j<newPosR)
        {
            int v = vertexArray[j];
            int isPruned = 1;
            for (int k = seps2[v]; k < inds[v+1]; k++)
            {
                int laterNeighbor = orderedVals[k];
                int laterNeiPos = vertexPos[laterNeighbor];
                if(laterNeiPos >= newPosP && laterNeiPos < newPosR)
                {
                    numEdges = numEdges + 2;
                    neighborsInP[v][numNeighbors[v]] = laterNeighbor;
                    numNeighbors[v]++;
                    neighborsInP[laterNeighbor][numNeighbors[laterNeighbor]] = v;
                    numNeighbors[laterNeighbor]++;
                }
                else{
                    isPruned = 0;
                }
            }
            if (isPruned && o2v[orders[v]] > 0)
            {
                o2v[orders[v]] = -1;
                Xpruned[v] = i;
            } 
            j++;
        }
        if (numEdges == (newPosR-newPosP)*(newPosR-newPosP-1))
        {
            cliqueCount++;
        }
        else{
            // degree zero filter
            int numDegZeroVertices = 0;
            bool degOneUpdated = false;
            for (int j = newPosP; j < newPosR;)
            {
                int v = vertexArray[j];
                if (numNeighbors[v] == 0)
                {
                    numDegZeroVertices++;
                    if (!bfsParent[v])
                    {
                        cliqueCount++;
                    }
                    newPosR--;
                    vertexArray[j] = vertexArray[newPosR];
                    vertexPos[vertexArray[newPosR]] = j;
                    vertexArray[newPosR] = v;
                    vertexPos[v] = newPosR;
                    continue;
                }
                else if (numNeighbors[v] == 1)
                {
                    if (!bfsParent[v]){ //can return a clique directly
                        degOneUpdated = true;
                        numDegZeroVertices++;
                        cliqueCount++;
                        newPosR--;
                        vertexArray[j] = vertexArray[newPosR];
                        vertexPos[vertexArray[newPosR]] = j;
                        vertexArray[newPosR] = v;
                        vertexPos[v] = newPosR;
                        int v2 = neighborsInP[v][0];
                        if (numNeighbors[v2] == 1)
                        {
                            int j2 = vertexPos[v2];
                            newPosR--;
                            vertexArray[j2] = vertexArray[newPosR];
                            vertexPos[vertexArray[newPosR]] = j2;
                            vertexArray[newPosR] = v2;
                            vertexPos[v2] = newPosR;
                            if (j2 < j && j < newPosR)
                            {
                                j--;
                                int tmp = vertexArray[j2];
                                vertexArray[j2] = vertexArray[j];
                                vertexPos[vertexArray[j]] = j2;
                                vertexArray[j] = tmp;
                                vertexPos[tmp] = j;
                            }
                        }
                        else {
                            for (int tmpIter = 0; tmpIter < numNeighbors[v2]; tmpIter++)
                            {
                                if (neighborsInP[v2][tmpIter] == v)
                                {
                                    numNeighbors[v2]--;
                                    neighborsInP[v2][tmpIter] = neighborsInP[v2][numNeighbors[v2]];
                                    break;
                                }
                            }  
                        }
                        continue;
                    }
                }
                j++;
            }
            if (newPosR == newPosP)
            {
                posR++;
                continue;
            }
            numEdges = 0;
            j = newPosP;
            int newNewPosR = newPosR;
            while(j<newNewPosR){
                int v = vertexArray[j];
                if (numNeighbors[v] == (newPosR-newPosP-1))
                {
                    if(!degOneUpdated && o2v[orders[v]] != -1) Xpruned[vertex] = min(orders[v], Xpruned[vertex]);
                    newNewPosR--;
                    vertexArray[j] = vertexArray[newNewPosR];
                    vertexPos[vertexArray[newNewPosR]] = j;
                    vertexArray[newNewPosR] = v;
                    vertexPos[v] = newNewPosR;
                    continue;
                }
                else 
                {
                    // P中的degree第一次初始化写在了这里
                    vertexDegree[v] = numNeighbors[v];
                    if (numNeighbors[v] > maxIntersectionSize)
                    {
                        pivot = v;
                        maxIntersectionSize = numNeighbors[v];
                        numEdges = numEdges + numNeighbors[v];
                    }

                }
                j++;
            }
            if (newNewPosR != newPosR || numDegZeroVertices > 0)
            {
                j = newPosX;
                pivot = -1;
                maxIntersectionSize = -1;
                while(j<newPosP){
                    int v = vertexArray[j];
                    int rCount = 0;
                    int zeroCount = 0;
                    for (int k = 0; k < numNeighbors[v];)
                    {
                        int neighbor = neighborsInP[v][k];
                        int neiPos = vertexPos[neighbor];
                        if (neiPos >= newPosP && neiPos < newNewPosR)
                        {
                            k++;
                            continue;
                        }
                        else if (neiPos >= newNewPosR && neiPos < newPosR){
                            neighborsInP[v][k] = neighborsInP[v][numNeighbors[v]-1];
                            numNeighbors[v]--;
                            rCount++;
                        }
                        else{
                            neighborsInP[v][k] = neighborsInP[v][numNeighbors[v]-1];
                            numNeighbors[v]--;
                            zeroCount++;
                        }
                        if (zeroCount == numDegZeroVertices && rCount == (newPosR - newNewPosR))
                        {
                            break;
                        }
                    }
                    // if (i == 558441) cout << v << " " << rCount << " " << newNewPosR-newPosP << endl;
                    if (rCount != (newPosR - newNewPosR))
                    {
                        vertexArray[j] = vertexArray[newPosX];
                        vertexPos[vertexArray[newPosX]] = j;
                        vertexArray[newPosX] = v;
                        vertexPos[v] = newPosX;
                        newPosX++;
                        j++;
                        continue;
                    }
                    if (numNeighbors[v] == newNewPosR-newPosP)
                    {
                        pivot = -2;
                        break;
                    }
                    if (numNeighbors[v] > maxIntersectionSize)
                    {
                        pivot = v;
                        maxIntersectionSize = numNeighbors[v];
                    }
                    j++;
                }
                if (pivot == -2)
                {
                    continue;
                }
                // 因为有点被移入R了，所以这里要重新计算一下Pneighbors
                j = newPosP;
                while (j < newNewPosR)
                {
                    int v = vertexArray[j];
                    for (int k = 0; k < numNeighbors[v];)
                    {
                        int neighbor = neighborsInP[v][k];
                        int neiPos = vertexPos[neighbor];
                        if (neiPos >= newPosP && neiPos < newNewPosR)
                        {
                            k++;
                            continue;
                        }
                        else{
                            neighborsInP[v][k] = neighborsInP[v][numNeighbors[v]-1];
                            numNeighbors[v]--;
                            vertexDegree[v] = numNeighbors[v]; //在这里更新
                        }
                    }
                    if (numNeighbors[v] > maxIntersectionSize)
                    {
                        pivot = v;
                        maxIntersectionSize = numNeighbors[v];
                    }
                    j++;
                }
                newPosR = newNewPosR;
            }
            int currentDepth = 1;
            int boundSize = maxIntersectionSize;
            partialClique.push_back(vertex);
            bpRecur(&cliqueCount,
                    partialClique, 
                    newPosX, newPosP, newPosR, pivot, boundSize, vertexDegree, currentDepth);
            partialClique.pop_back();
        }
        posR++;
        
    }

    partialClique.clear();

    end = clock();
    cout << "******************total time count: " << (double)(end-start)/CLOCKS_PER_SEC << "s" << endl;
    cout << "******************Number of maximal cliques: " << cliqueCount << endl;
    // cout << "******************Number of recursive calls: " << numberOfRecurrciveCall << endl;

    return cliqueCount;
}

void RMCEdegenMTP::bpRecur(long* cliqueCount,
                             list<int> &partialClique, 
                             int posX, int posP, int posR, int pivot, int boundSize, int* vertexDegree,int currentDepth){
    // numberOfRecurrciveCall++;
    
    if(posX >= posP && posP >= posR)
    {
        (*cliqueCount)++;
        return;
    }

    if(posP >= posR)
        return;

    int* candidates;
    int numCandidates = posR - posP;
    // vector<int> candidates(numCandidates,0);

    bool findPivot = true;
    if (numCandidates<=2)
    {
        findPivot = false;
    }

    if(findPivot)
    {
    setCandidates(&candidates,
                  &numCandidates,
                  posX, posP, posR, pivot, boundSize, vertexDegree,
                  cliqueCount, currentDepth);
    }
    else
    {
        candidates = (int*)malloc((posR-posP)*sizeof(int));
        memcpy(candidates, &vertexArray[posP], (posR-posP)*sizeof(int));
    }

    if(numCandidates != 0)
    {
        int iterator = 0;
        while(iterator < numCandidates)
        {
            int vertex = candidates[iterator];
            
            int newPosX, newPosP, newPosR;

            int newPivot = -1;
            int newBoundSize = -1;
            P2RandFindPivot(vertex, 
                            &posX, &posP, &posR, 
                            &newPosX, &newPosP, &newPosR, 
                            &newPivot, &newBoundSize, boundSize, cliqueCount, currentDepth);
            if (newPivot <= -3)
            {
               (*cliqueCount) = (*cliqueCount) + (-2 - newPivot);
            }
            else if (newPivot >= -1)
            {
                partialClique.push_back(vertex);
                bpRecur(cliqueCount,
                        partialClique,
                        newPosX, newPosP, newPosR, newPivot, newBoundSize, vertexDegree,currentDepth+1);
                partialClique.pop_back();
            }

            if(iterator == numCandidates-1)
            {
                posR = posR + 1;
                break;
            }

            R2X(vertex, 
                &posX, &posP, &posR);
            iterator++;
        }
        iterator = 0;

        while(iterator < numCandidates-1)
        {
            int vertex = candidates[iterator];
            int vPos = vertexPos[vertex];

            posP--;
            vertexArray[vPos] = vertexArray[posP];
            vertexArray[posP] = vertex;
            vertexPos[vertex] = posP;
            vertexPos[vertexArray[vPos]] = vPos;
            iterator++;
        }
    }
    free(candidates);
}

int main(int argc, char** argv) {
    string const filename = argv[1];
    RMCEdegenMTP alg(filename);
    alg.run();
}
