#include<iostream>
#include<vector>
#include<unordered_map>
#include<stack>
#include<queue>
#include<list>
#include<algorithm>
#include<limits.h>
#include<set>
//#ifndef ONLINE_JUDGE
#include <cstdio>
using namespace std;
void InputadjMatrix(vector<int>adjmat[],int V)
{
    for(int i=1; i<=V; i++)
    {
        for(int j=1; j<=V; j++)
        {
//            cin>>n;
			cin>>adjmat[i][j];
        }
    }
}
void BSFtravundir(vector<int>adjmat[],int V,int src,unordered_map<int,bool> &visited)
{
    queue<int>q;
    q.push(src);
    visited[src]=true;
   vector<int> distance;
   unordered_map<int,int>parent;
   
   int treeEdges=0;
   int crossEdge=0;
   parent[src] = -1;
    while(!q.empty())
    {
        int count=0;
        int front=q.front();
        q.pop();
       // cout<<front<<" ";
        for(int nbr =1; nbr<=V; nbr++)
        {
            
            if( visited[nbr]==false && adjmat[front][nbr]!=0)
            {
                q.push(nbr);
                visited[nbr]=true;
                //parent[nbr]=true;
                parent[adjmat[front][nbr]]=front;
                count++;
                treeEdges++;
            }
        }
        if(count!=0)
        {
            distance.push_back(count);
        }
        //cout<<count<<endl;
        //count=0;

    }
    for(auto i: distance)
    {
        //treeEdges=treeEdges+i;
        cout<<i<<" ";
    }
    cout<<0<<" ";
    int totaledges=0;
    for(int i=1; i<=V; i++)
    {
        for(int j=1; j<=V; j++)
        {
            if(adjmat[i][j]==1)
            {
                totaledges++;
            }
        }
    }
    crossEdge=((totaledges/2)-treeEdges);
    cout<<treeEdges<<" "<<crossEdge;
    cout<<endl;
    
}
void BFSDir(vector<int> adjmat[], int V, int src, unordered_map<int, bool> &visited) {
    queue<int> q;
    q.push(src);
    visited[src] = true;
    int x=0;
    vector<int> distance;
    vector<int> parent(V+1, -1);
    //vector<int> edgeType(V, -1); // -1 represents unclassified edges
    //int crossedge=0;
    int baed=0;
    int fred=0;
    int cred=0;
    while (!q.empty()) {
        int count=0;
        int front1 = q.front();
        q.pop();
        for (int nbr = 1; nbr <=V; nbr++) {
            if (!visited[nbr] && adjmat[front1][nbr] != 0) {
                q.push(nbr);
                visited[nbr] = true;
                count++;
                x++;

                parent[nbr] = front1;
                
            }
            else if (nbr != parent[front1]&& adjmat[front1][nbr] != 0) {
                    cred++;
                    //baed++;
                }
                //If the neighbor is visited and is the parent, it's a back edge
                else if (nbr == parent[front1] &&adjmat[front1][nbr] != 0) {
                    baed++;
                }

        }

        if(count!=0)
        {
            distance.push_back(count);
        }

    }
    

    for(auto i: distance)
    {
        cout<<i<<" ";
    }
    cout<<0<<" "<<x<<" "<<baed<<" "<<fred<<" "<<cred;


    // cout<<treeEdges<<" "<<crossEdge;
    cout<<endl;
}

void dfsdir(vector<int>adj[],int V,int src, unordered_map<int,bool>&visited,int &time,int &treeEdge,int &baed, int &fred, int &cred,int &time1,unordered_map<int,int>&start_time,unordered_map<int,int>&end_time)
{
    // unordered_map<int,int>start_time;
    // unordered_map<int,int>end_time;
    start_time[src] = time1;
    time1++;
    //int starttime=1;
   // cout<<src<<" ";
    visited[src]=true;
    for(int i=1; i<=V; i++)
    {
        if(adj[src][i]==1 && !visited[i])
        {
            time++;
            //cout<<time<<" ";
            treeEdge++;
            //time1++;
            dfsdir(adj,V,i,visited,time,treeEdge,baed,fred,cred,time1,start_time,end_time);
            end_time[src] = time1;
            time++;
            time1++;
            
        }
        else if(adj[src][i]==1 && visited[i])
        {
            
            if(start_time[src] < start_time[i] && end_time[src] > end_time[i])
            {
                // forward edge
                // cout << src << " -> " << i << " is a forward edge" << endl;
                
                fred++;
            }
            else if(start_time[src] > start_time[i] && end_time[src] < end_time[i])
            {
                // backward edge
                // cout << src << " -> " << i << " is a backward edge" << endl;
                baed++;
            }
            else if(start_time[src] < start_time[i] && end_time[src] < end_time[i])
            {
                // cross edge
                // cout << src << " -> " << i << " is a cross edge" << endl;
                //cout<<"Hi"<<endl;
                cred++;
            }
        }
    }
}

void dfsundir(vector<int>adj[],int V,int src, unordered_map<int,bool> &visited,int &time,int &treeEdge)
{
    // unordered_map<int,int>start_time;
    // unordered_map<int,int>end_time;
    //int starttime=1;
   // cout<<src<<" ";
   //start_time[src] = time1;
    //time1++;
    visited[src]=true;
    for(int i=1; i<=V; i++)
    {
        if(adj[src][i]==1 && !visited[i])
        {
            time++;
            //cout<<time<<" ";
            treeEdge++;
            dfsundir(adj,V,i,visited,time,treeEdge);
            time++;
            //cout<<"Hi";


        }
        // else if(adj[src][i]==1 && visited[i])
        // {
        //     //cout<<"Hi"<<endl;
        //     if(start_time[src] > start_time[i] && end_time[src] < end_time[i])
        //     {
        //         // backward edge
        //         // cout << src << " -> " << i << " is a backward edge" << endl;
        //         //cout<<"Hi"<<endl;
        //         baed++;
        //     }
        // }
    }
   // cout<<treeEdge<<" ";
    
}

void addEdge(unordered_map<int,list<int>> &adjList,int u, int v,int dir)
{
    if(dir==1)
    {
        adjList[u].push_back(v);
    }
    else
    {
        adjList[u].push_back(v);
        adjList[v].push_back(u);
    }
}
void addEdgefortopo(unordered_map<int,list<int>> &adjList,int u, int v)
{
    adjList[u].push_back(v);
}
bool checkCyclicDirectedGraphUsingDfs(unordered_map<int,list<int>> &adjList,int src, unordered_map<int,bool>& visited,unordered_map<int,bool> dfsVisited) 
{

		visited[src] = true;
		dfsVisited[src] = true;

		for(auto nbr: adjList[src]) {
			if(!visited[nbr]) {
				bool aageKaAnswer = checkCyclicDirectedGraphUsingDfs(adjList,nbr, visited, dfsVisited);
				if(aageKaAnswer == true)
					return true;
			}
			if(visited[nbr] == true && dfsVisited[nbr] == true) {
				return true;
			}
		}
		//yaha hi galti hoti h 
		dfsVisited[src] = false;
		return false;
	}
void toposortbfs(unordered_map<int,list<int> >adjList,int n,vector<int>&v)
{
    unordered_map<int,int>indegree;
    //queue<int>q;
    priority_queue<int,vector<int>,greater<int>>pq;
    //now calculate indegree
    for(auto i: adjList)
    {
        int src=i.first;
        for(auto nbr : i.second)
        {
            indegree[nbr]++;
        }
    }
    //now put node inside queue which has indegree zero
    for(int i=1; i<=n; i++)
    {
        if(indegree[i]==0)
        {
            pq.push(i);
        }
    }

    //now apply bfs logic to find the topological sort
    while(!pq.empty())
    {
        int fnode=pq.top();
        pq.pop();
       // cout<<fnode<<" ";
        v.push_back(fnode);
        for(auto nbr : adjList[fnode])
        {
            indegree[nbr]--;
            if(indegree[nbr]==0)
            {
                pq.push(nbr);
            }
        }

    }

}


//Dijkstra Algorithm
void addEdgeforDij(unordered_map<int,list<pair<int,int>>> &adjList ,int u, int v,int weight,int direction)
{
    //adjlist[u].push_back({v,weight});
    //let direction =1 for directed graph
    // if(direction==1)
    // {
    //     adjList[u].push_back({v,weight});
    // }
    adjList[u].push_back({v,weight});
    // if(direction!=1)
    // {
    //     //adjList[u].push_back({v,weight});
    //     adjList[v].push_back({u,weight});
    // }
}
void printadjList(unordered_map<int,list<pair<int,int>>> &adjList)
{
    for(auto i: adjList)
    {
        cout<<i.first<<"-> ";
        for(auto j: i.second)
        {
            cout<<"("<<j.first<<","<<j.second<<"),";
        }
        cout<<endl;
    }
}
//*******Implementing Dikstra Algorithm
void sortestDistDijkstra(unordered_map<int,list<pair<int,int>>> adjlist,int n,int src)
{
    for(int i=1; i<=n; i++)
    {
        for(auto nbr:adjlist[i])
        {
            if(nbr.second < 0)
            {
                cout<<-1<<endl;
                return;
            }
        }
    }
    vector<int>dist(n+1,INT_MAX);
    set<pair<int ,int>>st;
    //initial step 
    dist[src]=0;
    st.insert(make_pair(0,src));
    while(!st.empty())
    {
        //fetch the first or say smallest element from set
        auto topElement= *(st.begin());//note that set.begin() return iterator so we neeed to derefrence it so topElement will be our pair of integer
        int nodeDist=topElement.first;
        int nodeval=topElement.second;
        //because we take element from the set so we need to erase it
        st.erase(st.begin());//in erase function we send iteratro not value
        for(auto nbr: adjlist[nodeval])
        {
            if((nodeDist + nbr.second) < dist[nbr.first])
            {
                //now need to update distance and 
                //Please Node that
                //but if nbr already in set the need to update on set as well as distance array

                //find entry in set
                auto result = st.find(make_pair(dist[nbr.first],nbr.first));
                
                //if found then erase that iterater data
                if(result!=st.end())
                {
                    //need to erase
                    st.erase(result);

                }

                //uptdation in distance array and set
                dist[nbr.first]=(nodeDist + nbr.second);
                //Now need to insert in both case if erase the also and not erase mean it not in set so need to insert
                st.insert(make_pair(dist[nbr.first],nbr.first));

            }
        }
    }

    //cout<<"Now printing the value of distance array for ansewer: "<<endl;
    // for(auto i:dist)
    // {
    //     cout<<i<<" ";
    // }
    for(int i=1; i<=n; i++)
    {
        if(dist[i]==INT_MAX)
        {
            cout<<999999<<" ";
        }
        else
        {
            cout<<dist[i]<<" ";
        }
    }
    cout<<endl;




}
//Bellman Ford
void bellmanFord(unordered_map<int, list<pair<int, int>>> &adjList, int N, int source) {
    vector<int> distance(N + 1, 999999); // Initialize distances to a large value
    distance[source] = 0;
    int relax_operation=0;
    int modification=0;
    // Relax all edges N-1 times
    for (int i = 1; i <= N; ++i) {
        for (auto &it : adjList) {
            int u = it.first;
            //count++;
            for (auto &edge : it.second) {
                int v = edge.first;
                int weight = edge.second;

                relax_operation++;
                // if (distance[u] != 999999 && distance[u] + weight== distance[v]) {
                //     relax_operation++;
                // }
                if (distance[u] != 999999 && distance[u] + weight < distance[v]) {
                    distance[v] = distance[u] + weight;
                    modification++;
                }
            }
        }
    }
    
    // Check for negative weight cycles
    for (auto &it : adjList) {
        int u = it.first;

        for (auto &edge : it.second) {
            int v = edge.first;
            int weight = edge.second;

            if (distance[u] != 999999 && distance[u] + weight < distance[v]) {
                cout << -1 << endl;
                return;
            }
        }
    }

    // Print the shortest distances
    for (int i = 1; i <= N; ++i) {
        if (distance[i] == 999999) {
            cout << 999999 << " ";
        } else {
            cout << distance[i] << " ";
        }
    }
    cout<<relax_operation<<" "<<modification;
    cout << endl;
}

int main()
{
    //#ifndef 0NLINE_JUDGE
    // freopen("input.txt","r",stdin);
    // freopen("output1.txt","w",stdout);
    
    int T;
   // cout<<"Enter the number of test cases: ";
    cin>>T;
    //cout<<endl;
    while(T--)
    {

        int q;
       // cout<<"Enter query type: ";
        cin>>q;
        //cout<<endl;
        if(q==1)
        {
            int N,D,src;
            //cout<<"Enter Number of Nodes,direction,source: "<<endl;
            cin>>N>>D>>src;
            vector<int> adjmat[N+1];
            for(int i=1; i<=N ; i++)
            {
                for(int j=1; j<=N; j++)
                {
                    adjmat[i].push_back(0);
                }
            }
            //cout<<"Enter the matrix:"<<endl;
            InputadjMatrix(adjmat,N);
            unordered_map<int,bool>visited;
            if(D==1)
            {
                BFSDir(adjmat,N,src,visited);
            }
            else{
                BSFtravundir(adjmat,N,src,visited);
            }
        }
        else if(q==2)
        {
            int N,D,src;
            //cout<<"Enter Number of Nodes,direction,source: "<<endl;
            cin>>N>>D>>src;
            vector<int> adjmat[N+1];
            for(int i=1; i<=N ; i++)
            {
                for(int j=1; j<=N; j++)
                {
                    adjmat[i].push_back(0);
                }
            }
            //cout<<"Enter the matrix:"<<endl;
            InputadjMatrix(adjmat,N);
            unordered_map<int,bool>visited;
            int time=1;
            int treeEdge=0;
            int baed=0;
            int fred=0;
            int cred=0;
            int time1=1;
            if(D==1)
            {
                unordered_map<int,int>start_time;
                unordered_map<int,int>end_time;
                ///* dfsdirected done 
                dfsdir(adjmat,N,src,visited,time,treeEdge,baed,fred,cred,time1,start_time,end_time);
                cout<<(++time)<<" "<<treeEdge<<" "<<baed<<" "<<fred<<" "<<cred<<endl;
            }
            else
            {
                dfsundir(adjmat,N,src,visited,time,treeEdge);
                int totaledges=0;
                for(int i=1; i<=N; i++)
                {
                    for(int j=1; j<=N; j++)
                    {
                        if(adjmat[i][j]==1)
                        {
                            totaledges++;
                        }
                    }
                }
                baed=((totaledges/2)-treeEdge);
                cout<<(++time)<<" "<<treeEdge<<" "<<baed<<endl;
            }


        }
        else if(q==3)
        {
            int N;
            //cout<<"Enter the number of nodes: ";
            cin>>N;
            //cout<<endl;
            vector<int> adjmat[N+1];
            for(int i=1; i<=N ; i++)
            {
                for(int j=1; j<=N; j++)
                {
                    adjmat[i].push_back(0);
                }
            }
            //cout<<"Enter the matrix:"<<endl;
            InputadjMatrix(adjmat,N);
            unordered_map<int,bool>visited;
            unordered_map<int,list<int>>adjList;
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                {
                    if(adjmat[i][j]==1)
                    {
                        addEdgefortopo(adjList,i,j);
                    }
                }
            }
            //cout<<endl;
            // for(auto node: adjList)
            // {
            //     cout<<node.first<<"->";
            //     for(auto nbr:node.second)
            //     {
            //         cout<<nbr<<",";
            //     }
            //     cout<<endl;
            // }
            //cout<<endl;
            vector<int>v;
            toposortbfs(adjList,N,v);
            if(v.size()==N)
            {
                for(int i=0; i<v.size(); i++)
                {
                cout<<v[i]<<" ";
                }
                cout<<endl;
            }
            else
            {
                cout<<-1<<endl;
            }

        }
        else if(q==4)
        {
            int N,D,src;
            //cout<<"Enter Number of Nodes,direction,source: "<<endl;
            cin>>N>>D>>src;
            vector<int> adjmat[N+1];
            for(int i=1; i<=N ; i++)
            {
                for(int j=1; j<=N; j++)
                {
                    adjmat[i].push_back(0);
                }
            }
            //cout<<"Enter the matrix:"<<endl;
            InputadjMatrix(adjmat,N);
            unordered_map<int,bool>visited;
            unordered_map<int,list<pair<int,int>>>adjListDij;
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                {
                    if(adjmat[i][j]!=999999 && adjmat[i][j]!=0)
                    {
                        addEdgeforDij(adjListDij,i,j,adjmat[i][j],D);
                    }
                }
            }
           // printadjList(adjListDij);
            //cout<<endl;
            sortestDistDijkstra(adjListDij,N,src);
        }
        else if(q==5)
        {
            int N,D,src;
            //cout<<"Enter Number of Nodes,direction,source: "<<endl;
            cin>>N>>D>>src;
            vector<int> adjmat[N+1];
            for(int i=1; i<=N ; i++)
            {
                for(int j=1; j<=N; j++)
                {
                    adjmat[i].push_back(0);
                }
            }
            //cout<<"Enter the matrix:"<<endl;
            InputadjMatrix(adjmat,N);
            unordered_map<int,bool>visited;
            unordered_map<int,list<pair<int,int>>>adjListDij;
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                {
                    if(adjmat[i][j]!=999999 && adjmat[i][j]!=0)
                    {
                        addEdgeforDij(adjListDij,i,j,adjmat[i][j],D);
                    }
                }
            }
            //printadjList(adjListDij);
            bellmanFord(adjListDij,N,src);
            //cout<<endl;
        }

    }








}

