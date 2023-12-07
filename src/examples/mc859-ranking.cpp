
// Project and Analysis of Algorithms
// Flávio Keidi Miyazawa
#include <mc859-ranking.h>
  
// Uma heuristica ruim so' para exemplificar como construir solucao.
void Heuristica(Ranking_Entrada &C)
{
  DNode V[C.nnodes];
  // constroi uma solucao aleatoria
  int i=0,nit,melhorvalor,novovalor;
  for (DNodeIt v(C.g); v!=INVALID; ++v) {
    V[i] = v;
    i++;
    C.inserenaparte(v,(int) (drand48()*C.k));}

  // Tenta melhor a solucao atual com uma busca local probabilistica ingenua
  nit = 10000;
  melhorvalor = C.custosolucao();
  // repete nit vezes: Sorteia um vertice e coloca na melhor parte.
  for (int i=0;i<nit;i++) {
    DNode No;
    int p1,p2;
    No = V[(int) drand48()*C.nnodes]; // sorteia um no'
    p1 = C.solucao[No]; // guarda parte anterior
    p2 = (int) drand48()*C.k; // sorteia uma parte;
    if (p1==p2) continue;
    C.inserenaparte(No,p2);
    novovalor = C.custosolucao(); // pega o custo da nova solucao
    if (novovalor < melhorvalor) {
      melhorvalor = novovalor;
    }else{
      C.inserenaparte(No,p1);} // volta a parte anterior
  }
  //C.imprimesolucao();
  if (!C.verificasolucao()){
    cout << "Erro: Nao construiu solucao valida na fase 2."<<endl;exit(1);}
}


int main(int argc, char *argv[]) 
{
  //cout << "argc = " << argc << endl;
  string filename;
  if (argc!=2) {
    cout<<"Usage: "<< argv[0]<<" <filename>"<< endl<<endl<<
      "Example:" << endl <<
      "\t" << argv[0] << " "<< getpath(argv[0])+
      "../instances/mc859-ranking_1_10_100_0.5" << endl;
    exit(0);}
  filename = argv[1];
  Ranking_Entrada C(filename);
  // C.ImprimeFormatoArquivo();
  
  Heuristica(C);
  
  if (!C.verificasolucao()){
    cout<<"Erro: Ocorreu erro na producao de uma solucao."<<endl;}
  else{
    cout << "Custo da solucao obtida eh " << C.custosolucao() << endl;
  }
  C.imprimesolucao();
  C.visualizasolucao();
  return 0;
}
