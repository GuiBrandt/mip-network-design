
// Project and Analysis of Algorithms
// Flávio Keidi Miyazawa
// Problems with connectivity: Generate Triangulated Digraph
#include <mc859-compras.h>

int main(int argc, char *argv[]) 
{
  int seed;
  int box_width,box_height;
  Compras_Entrada C;
  double edge_factor;
  DNode V[MAXNODES];

  // double cutoff;   // used to prune non promissing branches (of the B&B tree)
  //cout << "argc = " << argc << endl;
  if (argc!=5) {
    cout<<"Usage: "<< argv[0]<<
      " <seed> <number_of_nodes_in_graph> <number_of_itens> <edge_factor>"<<endl<<endl<<
      " Example:" << endl <<
      argv[0]<< " 1 3 4 0.2"<< endl <<
      argv[0]<< " 1 100 100 0"<< endl <<endl<<
      "Obs.: Se <edge_factor> for igual a 0, o peso nao e' contado."<<endl;
    exit(0);}
  seed = atoi(argv[1]);
  srand48(seed);
  C.nnodes = atoi(argv[2]);
  C.nitens = atoi(argv[3]);
  if (C.nnodes>MAXNODES) {
    cout << "ERROR: maximum number of nodes is " << MAXNODES << endl;
    exit(1);
  }
  C.narcs = C.nnodes*(C.nnodes-1);
  if (C.nitens>MAXITENS) {
    cout << "ERROR: maximum number of itens is " << MAXITENS << endl;
    exit(1);
  }
  edge_factor = atof(argv[4]); // a factor between 0 and 1.
  box_width = 10000;// todos vertices tem coordenadas x,y em {0,...,10000}
  box_height = 10000;
  C.peso.resize(C.nitens);
  for (int i=0;i<C.nitens;i++) C.peso[i] = (int) (1.0 + drand48()*MAXPESO);
  
  for (int i=0;i<C.nnodes;i++) {
    DNode v=C.g.addNode();
    V[i] = v;
    C.vname[v] = IntToString(i);
    C.preco[v].resize(C.nitens);
    for (int j=0;j<C.nitens;j++)
      C.preco[v][j] = (int) (1.0 + drand48()*MAXPRECO);
    C.px[v] = (int) (1.0+drand48()*box_width);
    C.py[v] = (int) (1.0+drand48()*box_height);
  }
  C.Deposito = V[0];
  for (int i=0;i<C.nnodes;i++) {
    for (int j=0;j<C.nnodes;j++) {
      if (i!=j) {
	Arc a=C.g.addArc(V[i],V[j]);
	C.dist[a] = (int) (MAXDIST*sqrt(pow(C.px[V[i]]-C.px[V[j]],2)+pow(C.py[V[i]]-C.py[V[j]],2)));
	// Coloca alguma variacao nas distancias
	C.dist[a] = (int) (1.0 + 0.9*C.dist[a]  +  drand48()*0.1*C.dist[a]);
	C.fatorpeso[a] = (int) (drand48() * edge_factor * C.dist[a]);
      }
    }
  }
  for (int i=0;i<C.nitens;i++) C.preco[C.Deposito][i] = INT_MAX;
  //for (DNodeIt v(C.g); v!=INVALID; ++v) {
    //cout <<C.vname[v] << " " <<C.px[v] << " " <<C.py[v] << " ";
    //for (int i=0;i<C.nitens;i++) cout << C.preco[v][i] << " ";
    //cout << endl;
    //  }
  C.ImprimeFormatoArquivo();
  return 0;
}

