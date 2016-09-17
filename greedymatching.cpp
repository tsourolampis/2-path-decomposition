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



int deg1deg2=0;
void Degree1Degree2() 
{
  for(int i = 1; i <=E ; i++)
  {
     if( (degrees[eList[i][0] ]==1 && degrees[eList[i][1]]==2) || (degrees[eList[i][0] ]==2 && degrees[eList[i][1]]==1) )
			deg1deg2++;
  }
  printf("There exist %d pairs of deg1-deg2 edges\n",deg1deg2);
  return;
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
int KSGreedy = 0;
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

/* This function finds a maximal collection of disjoint paths. When a vertex
 v with degree 1 is connected to a vertex with degree 2, then this path is added
 to the collection. Otherwise we do as GreedyMatching()*/
void KSGreedyMatching()
{
    // we shall keep the degree 1 vertices with degree 2 neighbors in a special bucket. every time we
    // remove vertices then we can
    vector<int> Degree1Bucket;
	vector<int> Bucket;
    int u,w1,w2,w;
    int c12=0;
    for(int i=1; i<=V; i++)
    {
        if(degrees[i] == 1 && degrees[AdjList[i][0]] == 2)
        {
            Degree1Bucket.push_back(i);
            c12++;
        }
		else 
			Bucket.push_back(i);
		 
    }
    
    
    cout<<"There exist "<<Degree1Bucket.size()<<" vertices of degree 1, connected to vertices of degree 2 in G\n";
    Degree1Degree2();
    if( c12!= deg1deg2 )
    {
        cerr<<"Debug now!\n";
        exit(-1);
    }
    
    
    int active = V; // initially all vertices are active
    for(int mm=1; mm<=V; mm++)
        if(degrees[mm]==0)
            active--; //degree 0 vertices are considered inactive	
	//cout<<" Initially there are "<< active << " vertices."<<endl; 
	
    while(active>2)
    {
		//cout<<" Currently there exist "<< active << " vertices."<<endl; 
		
        while(Degree1Bucket.size())
        {
            // get vertex u with degree 1
            u = Degree1Bucket.back();
            Degree1Bucket.pop_back();
			w1 = AdjList[u][0]; 
            while( degrees[u]!=1 || degrees[w1]!=2 ) // during updates in later stages, a vertex u may have become inactive, but still be inside here
            {
				u = Degree1Bucket.back();
				Degree1Bucket.pop_back();
				w1 = AdjList[u][0]; // even if we make the degree of inactive  vertices of degree 1 zero, we do not modify their adjacency list 
			}
                
                
			if(degrees[u]==1 && degrees[w1]==2)
			{
					 degrees[u]  = 0;
					 active--; 
					 degrees[w1]  = 0;
					 active--; 
					if( AdjList[w1][0] == u)
						w2 = AdjList[w1][1];
					else
						w2 = AdjList[w1][0];
					degrees[w2] = 0; 
					active--; // now we have de-activated the vertices in the path u,w1,w2 
					cout<<"2-path ("<<u<<","<<w1<<","<<w2<<") chosen."<<endl;
					KSGreedy++; //increase the size of the matching
					
					// Update(w2) 
					int i; 
					for(i=0; i< AdjList[w2].size(); i++)
					{
						int Neighbor = AdjList[w2][i]; // the i-the neighbor of w2 
						for(int j = 0; j < AdjList[Neighbor].size(); j++ )
						{
							if(AdjList[Neighbor][j]==w2)
								AdjList[Neighbor].erase(AdjList[Neighbor].begin()+j);
						}
						
						if ( degrees[Neighbor] > 0 ) 
						{
							degrees[Neighbor]--;
							if( degrees[Neighbor] == 0 )
								active--; 
						}
					}
					
					for(i=0; i< AdjList[w2].size(); i++)
					{
						int x  = AdjList[w2][i];
						if( degrees[x] == 1 && degrees[ AdjList[x][0] ] == 2 )
						{
							Degree1Bucket.push_back(x);
						}
						else if(degrees[x]==2)
						{
							if(degrees[ AdjList[x][0] ] == 1 )
								Degree1Bucket.push_back( AdjList[x][0] );
							if(degrees[ AdjList[x][1] ] == 1 )
								Degree1Bucket.push_back( AdjList[x][1] );
						}
					}
			}
		} // end while (Degree1Bucket 
		
		/*at this point we are guaranteed that Degree1bucket is empty
		For this reason we continue popping things out of Bucket, until
		 * something fills in again Degree1bucket. This notification is 
		 * given by flag 
		 */
		 bool flag = true; 
		 
		 while(  Bucket.size() && flag )
		 {
			u =  Bucket.back();
		    Bucket.pop_back();
			//cout<<"Popped out u="<<u<<" with current degree " << degrees[u] << endl;
            while( degrees[u]==0 && Bucket.size()>0 )  
            {
				u =  Bucket.back();
				Bucket.pop_back();
			}
			//cout<< " Bucket size now is "<< Bucket.size()<< endl;
			if( degrees[u] > 0 )
			{
				// this is a single isolated edge
				if( degrees[u]==1 && degrees[AdjList[u][0]] == 1 )
				{
					degrees[u] = 0; 
					degrees[ AdjList[u][0] ] = 0; 
					active -= 2; 
				}
				else if(degrees[u]==1) 
				{
					w1 = AdjList[u][0]; 
					if( AdjList[w1][0] == u)
						w2 = AdjList[w1][1];
					else
						w2 = AdjList[w1][0];
					if(degrees[w1]>0 && degrees[w2]>0)
					{	
						degrees[u]   = 0;
						degrees[w1]  = 0;
						degrees[w2]  = 0;
						active-=3; 
						//cout<<"2-path ("<<u<<","<<w1<<","<<w2<<") chosen."<<endl;
						KSGreedy++; //increase the size of the matching
						
						
						// Update( ) 
						int i; 
						for(i=0; i< AdjList[w2].size(); i++)
						{
							int Neighbor = AdjList[w2][i]; // the i-the neighbor of w2 
							for(int j = 0; j < AdjList[Neighbor].size(); j++ )
							{
								if(AdjList[Neighbor][j]==w2)
									AdjList[Neighbor].erase(AdjList[Neighbor].begin()+j);
							}
							
							if(degrees[Neighbor]>0)
							{
								degrees[Neighbor]--;
								if( degrees[Neighbor] == 0 )
									active--; 
							}
						}
						
						
						for(i=0; i< AdjList[w2].size(); i++)
						{
							int x  = AdjList[w2][i];
							if( degrees[x] == 1 && degrees[ AdjList[x][0] ] == 2 )
							{
								Degree1Bucket.push_back(x);
								flag=false; 
							}
							else if(degrees[x]==2)
							{
								if(degrees[ AdjList[x][0] ] == 1 )
								{
									Degree1Bucket.push_back( AdjList[x][0] );
									flag=false;
								}
								
								if(degrees[ AdjList[x][1] ] == 1 )
								{
									Degree1Bucket.push_back( AdjList[x][1] );
									flag=false;
								}
							}
						}
					}
				}
				else if( degrees[u]>=2 )
				{
					int ind1,ind2; 
					do
					{
						ind1 = rand()%degrees[u];
						ind2 = rand()%degrees[u];
					}
					while(ind1 == ind2);
					KSGreedy++;
					int u1 = AdjList[u][ind1];
					int u2 = AdjList[u][ind2];
					//cout<<"2-path ("<<u1<<","<<u<<","<<u2<<") chosen."<<endl;
					degrees[u1]=0; // these three vertices will never be considered again
					degrees[u2]=0;
					degrees[u]=0;
					active-=3; 
					// update neighborhood of u 
					// Update(u) 
					int i; 
					for(i=0; i< AdjList[u].size(); i++)
					{
						int Neighbor = AdjList[u][i]; // the i-the neighbor of w2 
						for(int j = 0; j < AdjList[Neighbor].size(); j++ )
						{
							if(AdjList[Neighbor][j]==u)
								AdjList[Neighbor].erase(AdjList[Neighbor].begin()+j);
						}
						
						 if(degrees[Neighbor]>0)
						{
							degrees[Neighbor]--;
							if( degrees[Neighbor] == 0 )
								active--; 
						}
					}
					
					for(i=0; i< AdjList[u1].size(); i++)
					{
						int Neighbor = AdjList[u1][i]; // the i-the neighbor of w2 
						for(int j = 0; j < AdjList[Neighbor].size(); j++ )
						{
							if(AdjList[Neighbor][j]==u)
								AdjList[Neighbor].erase(AdjList[Neighbor].begin()+j);
						}
						
						if(degrees[Neighbor]>0)
						{
							degrees[Neighbor]--;
							if( degrees[Neighbor] == 0 )
								active--; 
						}
					}
					
					for(i=0; i< AdjList[u2].size(); i++)
					{
						int Neighbor = AdjList[u2][i]; // the i-the neighbor of w2 
						for(int j = 0; j < AdjList[Neighbor].size(); j++ )
						{
							if(AdjList[Neighbor][j]==u)
								AdjList[Neighbor].erase(AdjList[Neighbor].begin()+j);
						}
						
						if(degrees[Neighbor]>0)
						{
							degrees[Neighbor]--;
							if( degrees[Neighbor] == 0 )
								active--; 
						}
					}
					
					
					for(i=0; i< AdjList[u].size(); i++)
					{
						int x  = AdjList[u][i];
						if( degrees[x] == 1 && degrees[ AdjList[x][0] ] == 2 )
						{
							Degree1Bucket.push_back(x);
							flag = false; 
						}
						else if(degrees[x]==2)
						{
							if(degrees[ AdjList[x][0] ] == 1 )
							{
								Degree1Bucket.push_back( AdjList[x][0] );
								flag=false;
							}
							if(degrees[ AdjList[x][1] ] == 1 )
							{
								Degree1Bucket.push_back( AdjList[x][1] );
								flag=false; 
							}
						}
					}
					
					
					for(i=0; i< AdjList[u1].size(); i++)
					{
						int x  = AdjList[u1][i];
						if( degrees[x] == 1 && degrees[ AdjList[x][0] ] == 2 )
						{
							Degree1Bucket.push_back(x);
							flag = false; 
						}
						else if(degrees[x]==2)
						{
							if(degrees[ AdjList[x][0] ] == 1 )
							{
								Degree1Bucket.push_back( AdjList[x][0] );
								flag=false;
							}
							if(degrees[ AdjList[x][1] ] == 1 )
							{
								Degree1Bucket.push_back( AdjList[x][1] );
								flag=false; 
							}
						}
					}
					
					for(i=0; i< AdjList[u2].size(); i++)
					{
						int x  = AdjList[u2][i];
						if( degrees[x] == 1 && degrees[ AdjList[x][0] ] == 2 )
						{
							Degree1Bucket.push_back(x);
							flag = false; 
						}
						else if(degrees[x]==2)
						{
							if(degrees[ AdjList[x][0] ] == 1 )
							{
								Degree1Bucket.push_back( AdjList[x][0] );
								flag=false;
							}
							if(degrees[ AdjList[x][1] ] == 1 )
							{
								Degree1Bucket.push_back( AdjList[x][1] );
								flag=false; 
							}
						}
					}
					
				}
			}
		 }
		
		
	} // end while (active>0)
	cout<<"KSGreedy "<<KSGreedy<<endl;
} // end function 
					 
            

/* This function takes as input the graph and greedily adds 2-paths
 to the matching by sampling u-v-w uar from the set of all such 2-paths. */
void RandomGreedyMatching()
{
    
    
    
    
    
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
    /*CountTrianglesNaively();
  startTime = clock();
  PeelTriangles();
  cout<<"************** Triangle  PEELING ***************"<<endl;
  ElapsedTime(startTime,"Triangle Peeling"); 
  cout<<"Triangle Peeling's Results"<<endl; 
  cout<<"Avg. Triangle density:"  <<TRIANGLEDENSITY<<endl;
  cout<<"Size:"  <<TrianglePeelSize<<endl;
  cout<<"fe:"  <<TrianglePeelFe<<endl;
  cout<<"ft:"  <<TrianglePeelTe<<endl;
  cout<<"*******************************************"<<endl;
  ElapsedTime(startTime,"Triangle peeling"); 
*/
 
  return 0;
}
 
