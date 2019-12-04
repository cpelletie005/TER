#include <iostream>
#include <fstream> // Inclusion de "fstream"
using namespace std;
// #include <vtkstream>
#include "Eigen"
using namespace Eigen;

double deplacement_front(VectorXd rho, VectorXd rho_1,VectorXd rho_2, VectorXd maille, VectorXd maille_1, VectorXd maille_2, int nb_maille)
{
  //coucou
  double dx=1./nb_maille;
  // renvoie le deplacement du front de pyrolise entre les deux denieres iterations
  double rho_max(rho_2(nb_maille-1));
  double rho_min(rho_2(0));
  double rho_1max(rho_1(nb_maille-1));
  double rho_1min(rho_1(0));

  int i(0);
  while(rho_2(i)<((rho_max+rho_min)/2))
    i++;
  // cout << rho_2(0) << endl;
  // cout << rho_max << endl;
//   cout << i << endl;
// cout << "dx:" << dx << endl;
  double y2=2*i*dx;
  // cout <<"y2 =" << y2 << endl;
  i=0;
  while(rho_1(i)<((rho_1max+rho_1min)/2))
    i++;
  double y1=2*i*dx;

  // cout <<"y1 =" << y1;

  return(abs(y2-y1));
}


VectorXd Adaptation_maillage(VectorXd rho, VectorXd rho_1, VectorXd rho_2, VectorXd maille, VectorXd maille_1, VectorXd maille_2,int iteration, double nb_maille, double l, double dx)
{
  // maillage initial
  double a=1.5;
  double alpha=(1-1/(4*a))*4/3;

  if ( iteration <=5 )
    {

    for (int i=0 ; i<int(nb_maille/4) ; i++)
     {
       maille(i)=dx/a*i;
     }
    for (int i=int(nb_maille/4); i<nb_maille ; i++)
      {
        maille(i)=maille(i-1)+alpha*dx;
      }
      // cout << "coucou "<<endl;
      // cout << alpha*dx << endl;
      // cout << dx/a;
    }
  // maillage adapte chaque iteration
  else
  {
    double distance_deplacement=deplacement_front(rho,rho_1,rho_2 , maille, maille_1, maille_2, nb_maille);
    // cout << "distance_deplacement :" << endl;
    // cout << distance_deplacement << endl;
      // on adapte le maillage
      for(int i=0; i<nb_maille; i++)
        {maille(i)+=distance_deplacement;
        if (maille(i)>1)
          maille(i)-=1;
        }
      cout << maille(1) << endl;
  }
  return maille;
}

void save_maille(VectorXd maille,int iteration, int nb_maille)
{



  //
  ofstream mon_flux; // Contruit un objet "ofstream"
  string name_file("maille_" + to_string(iteration) + ".vtk"); // Le nom de mon fichier
  mon_flux.open(name_file, ios::out); // Ouvre un fichier appelé name_file


  if(mon_flux) // Vérifie que le fichier est bien ouvert
  {
    mon_flux << "# vtk DataFile Version 3.0" << endl; // Remplit le fichier
    mon_flux << "cell" << endl; // Remplit le fichier
    mon_flux << "ASCII" << endl; // Remplit le fichier
    mon_flux << "DATASET STRUCTURED_POINTS" << endl; // Remplit le fichier
    mon_flux << "DIMENSIONS" << " " << nb_maille << " 1" << " 1"<< endl; // Remplit le fichier
    mon_flux << "ORIGIN 0 0 0" << endl; // Remplit le fichier
    mon_flux << "SPACING 1 1 1" << endl; // Remplit le fichier
    mon_flux << "POINT_DATA" << " " << nb_maille << endl; // Remplit le fichier
    mon_flux << "SCALARS cell float" << endl; // Remplit le fichier
    mon_flux << "LOOKUP_TABLE default" << endl; // Remplit le fichier
    mon_flux.open("maille_" + to_string(iteration) + ".txt");
    for (int i=0; i<nb_maille;i++)
      {
        mon_flux << maille(i) <<  " 1" << endl ;
      }

    // for (int i=0; i<nb_maille;i++)
    //   {
    //     mon_flux << maille(i) <<  " 0" ;
    //     mon_flux << '\n';
    //   }

  }


}


int main()
{
  int nb_maille(500);
  VectorXd maille(nb_maille), rho(nb_maille);
  VectorXd maille_1(nb_maille), rho_1(nb_maille);
  VectorXd maille_2(nb_maille), rho_2(nb_maille);
  double i, rho_i, temp_i;
  double l=2;  // longueur du segment etudie
  double dx=1./nb_maille;


  for(int iteration=1; iteration<100; iteration++)
  {
    /////////////////////////////////////////////////////////////////
    //////         recuperation du fichier de valeur        /////////
    /////////////////////////////////////////////////////////////////

    //nom du fichier
    char num[10];
    sprintf(num, "%d", iteration);
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
    //////              adaptation de maillage              ///////// U= a completer
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


    save_maille(l*maille, iteration, nb_maille);

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
