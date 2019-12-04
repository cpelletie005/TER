#include <iostream>
#include <fstream> // Inclusion de "fstream"
using namespace std;
#include "Eigen"
using namespace Eigen;

double deplacement_front(VectorXd rho, VectorXd rho_1,VectorXd rho_2, VectorXd maille, VectorXd maille_1, VectorXd maille_2, int nb_maille)
{
  double dx=1./nb_maille;
  // renvoie le deplacement du front de pyrolise entre les deux denieres iterations
  double rho_max(rho_2(nb_maille-1));
  double rho_min(rho_2(0));
  double rho_1max(rho_1(nb_maille-1));
  double rho_1min(rho_1(0));
  // cout << "rho min " << rho_max << endl;

  int i(0);
  while(rho_2(i)<((rho_max+rho_min)/2))
    i++;
  // cout << rho_2(0) << endl;
  // cout << rho_max << endl;
//   cout << i << endl;
// cout << "dx:" << dx << endl;
  double y2=i*dx;
  // cout <<"y2 =" << y2 << endl;
  i=0;
  while(rho_1(i)<((rho_1max+rho_1min)/2))
    i++;
  double y1=i*dx;
  // cout <<"y1 =" << y1;

  return(abs(y2-y1));
}


VectorXd Adaptation_maillage(VectorXd rho, VectorXd rho_1, VectorXd rho_2, VectorXd maille, VectorXd maille_1, VectorXd maille_2,int iteration, double nb_maille, double l, double dx)
{
  // maillage initial
  double alpha=l/dx-l/8;

  if ( iteration <=5 )
    {

    for (int i=1 ; i<int(nb_maille/4) ; i++)
     {
       maille(i)=dx/2*i;
     }
    for (int i=int(nb_maille/4) ; i<nb_maille ; i++)
      {
        maille(i)=alpha*dx*i;
      }
      // cout << maille << endl;
    }
  // maillage adapte chaque iteration
  else
  {
    double distance_deplacement=deplacement_front(rho,rho_1,rho_2 , maille, maille_1, maille_2, nb_maille);
    cout << "distance_deplacement :" << endl;
    cout << distance_deplacement << endl;
      // on adapte le maillage
      for(int i=0; i<nb_maille; i++)
        maille(i)+=distance_deplacement;
  }
  // cout << "maille : " << endl;
  // cout << maille << endl;
  return maille;
}

int main()
{
  int nb_maille(100);
  VectorXd maille(nb_maille), rho(nb_maille);
  VectorXd maille_1(nb_maille), rho_1(nb_maille);
  VectorXd maille_2(nb_maille), rho_2(nb_maille);
  double i, rho_i, temp_i;
  double l=2;  // longueur du segment etudie
  double dx=l/nb_maille;


  for(int iteration=1; iteration<80; iteration++)
  {
    /////////////////////////////////////////////////////////////////
    //////         recuperation du fichier de valeur        /////////
    /////////////////////////////////////////////////////////////////

    //nom du fichier
    char num[10];
    sprintf(num, "%d", iteration+10);
    string name("solution_"), extension(".txt");
    name+=num+extension;
    // cout << name << endl;
    // lecture du fichier
    ifstream mon_flux(name);
    for (int ligne =0; ligne <nb_maille; ligne ++)
    {
      mon_flux >> i >> rho_i >> temp_i;
      // maille(ligne)=i;
      rho(ligne)=rho_i;
    }
    // 3ème façon de lire : caractère par caractère
    mon_flux.close();


    /////////////////////////////////////////////////////////////////
    //////              adaptation de maillage              /////////
    /////////////////////////////////////////////////////////////////


    maille = Adaptation_maillage(rho, rho_1, rho_2, maille, maille_1, maille_2, iteration, nb_maille, l, dx);
    // cout << maille << endl;
    /////////////////////////////////////////////////////////////////
    //////               decalage des indices               /////////
    /////////////////////////////////////////////////////////////////

    maille_2=maille_1;
    maille_1 =maille;

    rho_2=rho_1;
    rho_1=rho;

    // if (iteration ==10)
    // {
    //   cout << "1 -> " << endl << rho_2 << endl;
    //   cout << "2 -> " << endl << rho_1 << endl;
    //   cout << "3 -> " << endl << rho << endl;
    //
    // }
    // cout << "iteration : " << endl;
    // cout << iteration << endl;
    // cout << maille << endl;


  }
}
