
// Flávio Keidi Miyazawa
#include <mc859-strip-packing.h>
#include <stdlib.h>

int main(int argc, char *argv[]) 
{
  int seed;
  int box_width,box_height;
  Ranking_Entrada C;
  double random_factor;
  DNode V[MAXNODES];

  // double cutoff;   // used to prune non promissing branches (of the B&B tree)
  //cout << "argc = " << argc << endl;
  if (argc!=5) {
    cout<<"Usage: "<< argv[0]<<
      " <seed> <number_of_rectangles> <strip_width> <lower_item_width> <upper_item_width>"<<endl<<endl<<
      " Example:" << endl <<
      argv[0]<< " 1 10  100 0 100"<< endl <<
      argv[0]<< " 1 100 1000 251 750"<< endl <<
      argv[0]<< " 1 3 20 0.5"<< endl <<
      endl<<
      "Obs.: Se <random_factor> for igual a 0, os pesos sao mais influenciados pela distancia."<<endl;
    exit(0);}
  seed = atoi(argv[1]);
  srand48(seed);
  C.k = atoi(argv[2]);
  C.nnodes = atoi(argv[3]);
  if (C.nnodes>MAXNODES) {
    cout << "ERROR: maximum number of nodes is " << MAXNODES << endl; exit(1); }
  C.narcs = C.nnodes*(C.nnodes-1)/2;
  random_factor = atof(argv[4]); // a factor between 0 and 1.
  box_width  = MAXDIST; // todos vertices tem coordenadas x,y em {0,...,MAXDIST}
  box_height = MAXDIST;
  
  for (int i=0;i<C.nnodes;i++) {
    DNode v=C.g.addNode();
    V[i] = v;
    C.vname[v] = IntToString(i);
    C.px[v] = (int) (1.0+drand48()*box_width);
    C.py[v] = (int) (1.0+drand48()*box_height);
  }
  for (int i=0;i<C.nnodes-1;i++) {
    for (int j=i+1;j<C.nnodes;j++) {
      Arc a=C.g.addArc(V[i],V[j]);
      int dist = (int) (drand48()*sqrt(pow(C.px[V[i]]-C.px[V[j]],2)+pow(C.py[V[i]]-C.py[V[j]],2)))+ (drand48()*sqrt(pow(C.px[V[i]],2)+pow(C.py[V[i]],2)));
	// Coloca alguma variacao nas distancias
      //cout << dist <<" ";
      C.antes[a]  = (int) (dist  +  drand48()*random_factor*dist);
      C.depois[a] = (int) (drand48()*random_factor*(2*MAXDIST - C.antes[a]));
      C.igual[a]  = (int) (drand48()*random_factor*(C.antes[a] + C.depois[a])/2.0);
      //cout << C.antes[a]<<"__";
      //cout << C.igual[a]<<"__";
      //cout << C.depois[a]<<endl;
    }
  }
  C.ImprimeFormatoArquivo();
  return 0;
}

