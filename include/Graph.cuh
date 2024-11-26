#pragma once
#include <bits/stdc++.h>
#include <cmath>
#include <ctime>
#include <iostream>
#include <queue>
#include <vector>
#include <thrust/host_vector.h>

namespace cuDFNsys
{
    class Graph
    {
    public:
        // The number of fractures
        int V;

        // Adjacent list
        std::list<int> *Adj;

    public:
        // constructor
        Graph(const int &NumOfFractures_1,
              const thrust::host_vector<int2> &IntersectionFracturePairs);
        // implement DFS
        void UseDFS(thrust::host_vector<thrust::host_vector<int>> &S);
        // destructor
        ~Graph()
        {
            delete[] this->Adj;
            this->Adj = NULL;
        };

    private:
        void DFS(thrust::host_vector<thrust::host_vector<int>> &ListOfClusters);
        void DFSUtil(int s, thrust::host_vector<bool> &visited, thrust::host_vector<int> &onecluster);
        void addEdge(int v, int w);
    };

    Graph::Graph(const int &NumOfFractures_1,
                 const thrust::host_vector<int2> &IntersectionFracturePairs)
    {
        Adj = new std::list<int>[NumOfFractures_1];
        V = NumOfFractures_1;
        // std::map<pair<int, int>, pair<cuDFNsys::Vector3<T>,
        //                                     cuDFNsys::Vector3<T>>>::iterator i;
        for (auto i : IntersectionFracturePairs)
        {
            addEdge(i.x, i.y);
            addEdge(i.y, i.x);
        }
    }; // Graph

    void Graph::DFS(thrust::host_vector<thrust::host_vector<int>> &ListOfClusters)
    {
        // Mark all the vertices as not visited
        thrust::host_vector<bool> visited(V, false);
        ListOfClusters.reserve(V);

        for (int i = 0; i < V; i++)
        {
            thrust::host_vector<int> onecluster;
            onecluster.reserve(V);

            if (!visited[i])
            {
                DFSUtil(i, visited, onecluster);
            }

            onecluster.shrink_to_fit();

            if (onecluster.size() > 0)
            {
                ListOfClusters.push_back(onecluster);
            }
        }
        ListOfClusters.shrink_to_fit();
    };

    void Graph::DFSUtil(int s, thrust::host_vector<bool> &visited, thrust::host_vector<int> &onecluster)
    {
        // Create a stack for DFS
        std::stack<int> stack;

        // Push the current source node.
        stack.push(s);

        while (!stack.empty())
        {
            // Pop a vertex from stack and print it
            s = stack.top();
            stack.pop();

            // Stack may contain same vertex twice. So
            // we need to print the popped item only
            // if it is not visited.
            if (!visited[s])
            {
                // cout << s << " ";
                onecluster.push_back((int)s);
                visited[s] = true;
            }

            // Get all adjacent vertices of the popped vertex s
            // If a adjacent has not been visited, then push it
            // to the stack.
            for (auto i = Adj[s].begin(); i != Adj[s].end(); ++i)
            {
                if (!visited[*i])
                {
                    stack.push(*i);
                };
            }
        }
    }; // DFSUtil

    inline void Graph::UseDFS(thrust::host_vector<thrust::host_vector<int>> &S)
    {
        S.clear();
        DFS(S);
        S.shrink_to_fit();
    };

    inline void Graph::addEdge(int v, int w)
    {
        // Add w to v’s list.
        Adj[v].push_back(w);
    };
};