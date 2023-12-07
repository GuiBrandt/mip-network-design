
// Project and Analysis of Algorithms
// Flávio Keidi Miyazawa
#include <mc859-compras.h>

void TrocaDNode(DNode &a,DNode &b){DNode aux; aux=a; a=b; b=aux;}
  
// Uma heuristica ruim so' para exemplificar como construir solucao.
void Heuristica(Compras_Entrada &C)
{
  vector <int> LojaEscolhida(MAXITENS);
  int menorpreco,nit;
  long long int valorsolucaoatual,valorsolucaonova;
  DNodeBoolMap TemItemComprado(C.g);
  DNode No;
  for (DNodeIt v(C.g); v!=INVALID; ++v) TemItemComprado[v] = false;
  for (int i=0;i<C.nitens;i++){
    menorpreco = INT_MAX;
    No = INVALID;
    for (DNodeIt v(C.g); v!=INVALID; ++v) {
      //cout << "Preco do item "<<i<<" na loja "<<C.vname[v]<<": "<<C.preco[v][i]<<endl;
      if (C.preco[v][i] < menorpreco){
	No = v;
	menorpreco = C.preco[v][i];
      }
    }
    if (No==INVALID) {
      cout << "Erro: Nao eh possivel comprar item " << i << "." << endl;
      exit(3);
    }
    C.lojaitem[i] = No;
    TemItemComprado[No] = true;
    //cout << "Item "<< i <<" comprado na loja "<< C.vname[C.lojaitem[i]] <<
    //  " por " << C.preco[C.lojaitem[i]][i] << "." << endl;
  }
  C.esvaziasolucao();
  // Primeiro noh da solucao eh o deposito.
  C.inserelojasolucao(C.Deposito);
  // Insere na solucao todas lojas onde comprou algo
  for (DNodeIt v(C.g); v!=INVALID; ++v) 
    if (TemItemComprado[v]) C.inserelojasolucao(v);
  // Ultimo noh da solucao eh o deposito.
  C.inserelojasolucao(C.Deposito);

  // Tenta melhor a solucao atual com uma busca local probabilistica ingenua
  nit = 10000;
  // repete nit vezes: Troca dois no's quaisquer. Se troca nao melhorou, volta atras.
  valorsolucaoatual = C.custosolucao();
  for (int i=0;i<nit;i++) {
    if (i%(nit/10)==0) cout << "It. "<<i<<" de "<<nit<<endl;
    int ia,ib;
    ia = 1+(int) (drand48()*(C.nnodesequence-2));
    ib = 1+(int) (drand48()*(C.nnodesequence-2));
    //if ((ia==ib)||(ia<1)||(ia>C.nnodesequence-2)||(ib<1)||(ib>C.nnodesequence-2)) continue;
    TrocaDNode(C.nodesequence[ia],C.nodesequence[ib]);
    valorsolucaonova=C.custosolucao();
    if (valorsolucaoatual < valorsolucaonova){
      TrocaDNode(C.nodesequence[ia],C.nodesequence[ib]);
    }else{
      valorsolucaoatual=valorsolucaonova;
    }
  }
  //C.imprimesolucao();
  if (!C.verificasolucao()){
    cout << "Erro: Nao construiu solucao valida na fase 2."<<endl;exit(1);}
}


int main(int argc, char *argv[]) 
{
  //cout << "argc = " << argc << endl;
  if (argc!=2) {
    cout<<"Usage: "<< argv[0]<<" <filename>"<< endl<<endl<<
      "Example:" << endl <<
      "\t" << argv[0] << " "<< getpath(argv[0])+
      "../instances/mc859-compras_1_100_100_0" << endl;
    exit(0);}

  Compras_Entrada C(argv[1]);
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
