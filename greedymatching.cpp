#include <vector>
//#include <unordered_map>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <cstdlib>
#include "IO.h"
#include "HASH.h"

using namespace std; 
//using std::unordered_map;

/* Author: Charalampos E. Tsourakakis, babis@seas.harvard.edu*/


const int MAXV = 10000000;
const int MAXE = 60000000;
const int QUERYBUFFER = 6000000;

int  MAXDEG=0; 
// we shall assume that all vertices are numbered  from 1 to n 
// For convenience, we also name the edges e_1 through e_m

int E, V;
int eList[MAXE+1][2]; // we load all edges on memoer
int degrees[MAXV+1];  
int NumTriangles=0;
// remove multiple edges, self-loops
void SimplifyGraph() 
{
  int E1 = 1;
  int multipleEdges=0;
  HASH::init();
  for(int i = 1; i <= E; ++i) {
    //printf("Testing edge (%d,%d)\n", eList[i][0], eList[i][1]);
    if(1 <= eList[i][0] && eList[i][0] <= V &&
         1 <= eList[i][1] && eList[i][1] <= V &&
         eList[i][0] != eList[i][1] &&
         HASH::insert(eList[i][0], eList[i][1])) {
      eList[E1][0] = eList[i][0];
      eList[E1][1] = eList[i][1];
	  degrees[eList[i][0]]++;
	  degrees[eList[i][1]]++; 
      E1++;
    }
	else
		multipleEdges++; //printf("has appeared before\n");
  }
  E = E1-1;
  printf("number of edges in simple graph: %d\n", E);
  printf("%d\n",multipleEdges);
}



void GraphIn() {
  int u, v;
  IO::read(V);
  IO::read(E);
  printf("Number of vertices and edges (%d,%d)\n",V,E);
  for(int i = 1; i <= E; ++i) 
  {
    IO::read(u);
    IO::read(v);
	if(v>u)
	{
		eList[i][0] = u;
		eList[i][1] = v;
	}
    if(u>v)	
	{
		eList[i][0] = v;
		eList[i][1] = u;
	}
  }
}

void PrintEdgeList()
{
   for(int i=1; i<=E; i++) 
   {
        printf("(%d,%d)\n",eList[i][0],eList[i][1]); 
   }
}

void PrintDegrees()
{ 
   printf("***************************\n");
   printf("Vertex id \t Degree\n"); 
   for(int i = 1; i <= V; i++)
      printf("%d \t %d\n",i, degrees[i]); 
   printf("***************************\n");
}



void MaximumDegree()
{  
   for(int i = 1; i <= V; i++)
    if(MAXDEG < degrees[i] )
    	MAXDEG=degrees[i];
    
}


void graphinStdio() {
  scanf("%d%d", &V, &E);
  for(int i = 0; i < E; ++i) {
    scanf("%d%d", &eList[i][0], &eList[i][1]);
  }
}

void ElapsedTime(double startTime)
{ 
   printf("Elapsed time to read input: %f\n", (clock()-startTime)/CLOCKS_PER_SEC );
}

void ElapsedTime(double startTime, char* s)
{ 
   printf("Elapsed time form method %s: %f\n", s, (clock()-startTime)/CLOCKS_PER_SEC );
}
///////////////adjacency list
int eStart[MAXV];
int eNext[MAXE];

const int NOEDGE = -1;

void BuildGraph() {
  for(int i = 1; i <= V; ++i) {
    eStart[i] = NOEDGE;
  }
  for(int i = 1; i <= E; ++i) {
    eNext[i] = eStart[eList[i][0]];
    eStart[eList[i][0]] = i;
  }
}



 
vector<int> AdjList[MAXV];
 
void BuildAdjList()
{ 
  for(int i = 1; i<=E; i++)
  {
       AdjList[ eList[i][0] ].push_back( eList[i][1] ); 
	   AdjList[ eList[i][1] ].push_back( eList[i][0] );   
  }
}

// This function deletes a vertex from the graph that has the minimum degree 

int permutation[MAXV+1]; //permutation[1] through permutation[V] contains the desired perm 

double  EDGEDENSITY = -1; 
int CharikarSize=-1; 
double CharikarFe = 0.0; 


void SampleWedge()
{
    double r = ((double) rand() / (RAND_MAX));
}


int Greedy = 0 ;
int RGreedy = 0 ;

/* This function takes as input the graph and greedily adds 2-paths
 to the matching in an arbitrary way. */
void GreedyMatching()
{
    int ind1,ind2;
    // go through the vertices in order and if the degree is greater than 2
    // then choose two random neighbors. Increase the size of the matching by 1
    // and then destroy the three matched vertices
    for(int u = 1; u<=V; u++)
    {
        if(degrees[u]>=2)
        {
            do
            {
                 ind1 = rand()%degrees[u];
                 ind2 = rand()%degrees[u];
            }
            while(ind1 == ind2);
            Greedy++;
            
            int u1 = AdjList[u][ind1];
            int u2 = AdjList[u][ind2];
            //cout<<"2-path ("<<u1<<","<<u<<","<<u2<<") chosen."<<endl;
			//cout<<u1<<" "<<u<<" "<<u2<<""<<endl;
            degrees[u1]=0; // these three vertices will never be considered again
            degrees[u2]=0;
            degrees[u]=0;
            
            
            // now let's erase these three vertices and
            // update the graph and the degrees accordingly
            int w,i=0;
            while(i < AdjList[u].size())
            {
                if( i !=ind1 && i !=ind2 )
                {
                    w = AdjList[u][i];
                    for(int j = 0; j < AdjList[w].size(); j++)
                    {
                        if( AdjList[w][j] == u   )
                        {
                            degrees[w]--;
                            AdjList[w].erase(AdjList[w].begin()+j);
                            break;
                        }
                    }
                }
                i++;
            }
            
            i = 0;
            while( i < AdjList[u1].size() )
            {
                w = AdjList[u1][i];
                if(w!=u && w!=u2)
                {
                    for(int j = 0; j < AdjList[w].size(); j++)
                    {
                        if(   AdjList[w][j] == u1   )
                        {
                            degrees[w]--;
                            AdjList[w].erase(AdjList[w].begin()+j);
                            break;
                        }
                    }
                }
                i++;
            }
            
            i = 0;
            while( i < AdjList[u2].size() )
            {
                w = AdjList[u2][i];
                if(w!=u && w!=u1)
                {
                    for(int j = 0; j < AdjList[w].size(); j++)
                    {
                        if(  AdjList[w][j] == u2 )
                        {
                            degrees[w]--;
                            AdjList[w].erase(AdjList[w].begin()+j);
                            break;
                        }
                    }
                }
                i++;
            } // end while
            
        } // end if
    }
    cout<<"Greedy matching results in a matching of size "<<Greedy<<endl;
    
}




void PrintAdjList()
{
	for(int j = 1; j<=V; j++)
	{
	    cout<< "Neighbors of "<<j<<" [";
		for (int i=0; i<AdjList[j].size();i++)
		{
		cout << AdjList[j][i]<< " ";
		}
        cout<< "]" <<endl;
		
	}
}




int main(int argc, char **argv)
{
  double startTime = clock();
  GraphIn();
  ElapsedTime(startTime);
  //PrintEdgeList();
  SimplifyGraph(); 
 // PrintDegrees();
 // MaximumDegree();
  BuildGraph();
  BuildAdjList();
  //PrintAdjList();
    GreedyMatching();
 
 
  return 0;
}
 
