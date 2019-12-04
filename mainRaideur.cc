#include <iostream>
#include <math.h>
#include "Eigen"

using namespace std;

Eigen::MatriXd adv_elem1D(const double nu, const double c,
  int it,Eigen::VectorXd p,Eigen::VectorXd t)
{
  Eigen::MatriXd Ke(2,2); Ke.setZero();
  double i1=t(it,1), i2=t(it,2);
  double x1=p(i1,1), x2=p(i2,1);
  double l=abs(x1-x2);

  Ke(0,0)=nu/l - c/2; Ke(0,1)=-nu/l + c/2;
  Ke(1,0)=-nu/l - c/2; Ke(1,1)=nu/l + c/2;

  return Ke;
}
Eigen::MatriXd mdv_elem1D(const double nu, const double c,
  int it,Eigen::VectorXd p,Eigen::VectorXd t)
{
  Eigen::MatriXd Ke(2,2); Ke.setZero();
  double i1=t(it,1), i2=t(it,2);
  double x1=p(i1,1), x2=p(i2,1);
  double l=abs(x1-x2);

  Me(0,0)=l/3; Me(0,1)=l/6;
  Me(1,0)=l/6; Ke(1,1)=l/3;

  return Me;
}
Eigen::MatriXd assemblage_raideur(const int N,const double nu, const double c,
  Eigen::VectorXd p,Eigen::VectorXd t)
{
  Eigen::MatriXd K(N+1,N+1); K.setZero();
  for (int elem=0; elem<N; elem++)
  {
    Eigen::MatriXd Kelem=adv_elem1D(nu, c, elem, p, t);
    for (int i=0; i<2; i++)
    {
      for (j=0; j<2; j++)
      {
        K(t(elem,i),t(elem,j))+=Kelem(i,j);
      }
    }
  }
  K(0,0)=1e12; K(N,N)=1e12;

  return K;
}
Eigen::MatriXd assemblage_masse(const int N,const double nu, const double c,
  Eigen::VectorXd p,Eigen::VectorXd t)
{
  Eigen::MatriXd M(N+1,N+1); M.setZero();
  for (int elem=0; elem<N; elem++)
  {
    Eigen::MatriXd Kelem=mdv_elem1D(nu, c, elem, p, t);
    for (int i=0; i<2; i++)
    {
      for (j=0; j<2; j++)
      {
        M(t(elem,i),t(elem,j))+=Melem(i,j);
      }
    }
  }
  return M;
}

Eigen::VectorXd compute_derivative(const int N,
  Eigen::VectorXd U,Eigen::VectorXd p,Eigen::VectorXd t)
{
  Eigen::VectorXd D(N+1); D.setZero();

  for (int i=1; i<N; i++)
    D(i)=(U(i+1)-U(i-1))/(p(i+1)-p(i-1));
  D(0)=(U(1)-U(0))/(p(1)-p(0));
  D(N)=(U(N)-U(N-1))/(p(N)-p(N-1));

  return D;
}
Eigen::VectorXd compute_so_derivative(const int N,
  Eigen::VectorXd U,Eigen::VectorXd p,Eigen::VectorXd t)
{
  Eigen::VectorXd D2(N+1); D.setZero();
  Eigen::VectorXd D(N+1);D=compute_derivative(N,U,p,t);
  double l1(0),uppi1(0),l2(0),uppi2(0);

  for (int i=1; i<N; i++)
  {
    l1=p(i+1)-p(i); upp1=(U(i+1)-U(i)-D(i)*l1)*2/pow(l1,2);
    l2=p(i-1)-p(i); upp2=(U(i-1)-U(i)-D(i)*l2)*2/pow(l2,2);
    D2=(uppi1*l1 - uppi2*l2)/(l1-l2);
  }
  D2(0)=(U(1)-U(0)-D(0)*(p(1)-p(0)))*2/pow(p(1)-p(0),2);
  D2(N)=(U(N-1)-U(N)-D(N)*(p(N-1)-p(N)))*2/pow(p(N)-p(N-1),2);

  return D2;
}

Eigen::VectorXd compute_metric(const int N,
  Eigen::VectorXd U,Eigen::VectorXd p,Eigen::VectorXd t)
{
  Eigen::VectorXd M(N+1); M.setZero();
  Eigen::VectorXd D2;

  D2=compute_so_derivative(N,U,p,t);
  for (int i=0; i<N+1; i++)
  {
    M(i)=max(abs(D2(i)),1);
    M(i)=sqrt(M(i));
  }

  return M;
}
Eigen::MatriXd met_elem1D(Eigen::VectorXd M,
  int it,Eigen::VectorXd p,Eigen::VectorXd t)
{
  Eigen::MatriXd Me(2,2);
  double ki=(M(it-1)+M(it))/2;
  Me(0,0)=ki; Me(0,1)=-ki; Me(1,0)=-ki; Me(1,1)=ki;
  return Me;
}
Eigen::MatriXd assemblage_metrique(const int N,Eigen::VectorXd M,
  Eigen::VectorXd p,Eigen::VectorXd t)
{
  Eigen::MatriXd Met(N+1,N+1); Met.setZero();
  Eigen::MatriXd Melem;

  for(int elem=0; elem<N; elem++)
  {
    Melem=met_elem1D(M, elem, p, t);
    for (int i=0; i<2; i++)
    {
      for (int j=0; j<2; j++)
      {
        Met(t(elem,i),t(elem,j))+=Melem(i,j);
      }
    }
  }
  Met(0,0)=1e12; Met(N,N)=1e12;

  return Met;
}

Eigen::VectorXd solve_system(const int N,const double nu, const double c,
   const double L, Eigen::VectorXd p,Eigen::VectorXd t)
{
  Eigen::VectorXd U(N+1);
  Eigen::MatriXd K, M;

  K=assemblage_raideur(N,nu,c,p,t);
  M=assemblage_masse(N,nu,c,p,t);

  Eigen::VectorXd B(N+1),F(N+1); F.setOnes();
  B=M*F; B(1)=0; B(N)=0;

  // TODO Résolution de K*U=B
  U(0)=0; U(N)=0;

  return U;
}

int main()
{
  // Définition des constantes
  const double nu(0.01), c(0.5), L(1);
  const int N(50);

  // Construction vecteur de correpondances élément -> indices sommets
  Eigen::MatriXd corres(N,2); corres.setZero();
  for (int i=0, i<N, i++)
  {
    corres(i,0)=i+1;
    corres(i,1)=i+2;
  }
  // Construction table de coordonnées
  Eigen::VectorXd coord(N+1);
  for (int i=0; i<N+1; i++)
    coord(i)=L*i/N;

  Eigen::VectorXd U(N+1); U.setZero();
  Eigen::VectorXd X(N+1), new_coord(N+1), D(N+1);
  Eigen::VectorXd Metrique;
  Eigen::MatriXd Mat_metrique;

  int test(1), nb_iter(0);
  while((test==1)&&(nb_iter<50))
  {
    U=solve_system(N, nu, c, L, coord, corres);
    Metrique=compute_metric(N,U,coord,corres);
    Mat_metrique=assemblage_metrique(N,Metrique,coord,corres);

    X.setZero();
    X(N)=1e12;
    // TODO Résoudre Mat_metrique*new_coord=X
    new_coord(0)=0; new_coord(N)=1;

    D=new_coord-coord;
    if(D.norm()<1e-2)
    {
      test=0;
    }
    else
    {
      coord=new_coord;
      nb_iter+=1;
    }
  }

  return 0;
}
