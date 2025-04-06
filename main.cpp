#include <iostream>
#include "Eigen/Eigen"
#include <iomanip>

using namespace std;
using namespace Eigen;

/*
	Ho scelto di utilizzare Palu e non full perchè...
*/
bool MatrixIsSingular(const Matrix2d& A)
{
    JacobiSVD<Matrix2d> svd(A);
    Vector2d singularValuesA = svd.singularValues();

    if( singularValuesA.minCoeff() < 1e-16)
        return false;
    else
        return true;
}


void Err_rel_PALU(const Matrix2d& A, const Vector2d& b, const Vector2d& sol){
	cout<<"Calcolo l'errore relativo della soluzione del sistema lineare Ax = b, usando decomposizione PALU e QR"<<endl;
	cout<<"Matrice A:"<<endl;
	cout <<setprecision(16) <<scientific << A << endl;
	cout<<"Vettore b:"<<endl;
	cout<<b<<endl;
	
	// NB: occorre verificare se la matrice è invertibile peerchè...
	
	if (MatrixIsSingular(A)){
		Vector2d x=A.lu().solve(b);
		
		double er_rel = (sol-x).norm()/sol.norm();
		cout<<"Erorre relativo usando PALU: "<<er_rel<<endl;
	}else{
		cout<<"La matrice relativa a questo sistema lineare non è invertibile"<<endl;
	}
}

// return true if the matrix A is invertible, else it returns false

int main()
{
	
	
	Vector2d sol = -Vector2d::Ones();
	Matrix2d A1 = Matrix2d::Zero();
	A1 << 5.547001962252291e-01, -3.770900990025203e-02, 8.320502943378437e-01, -9.992887623566787e-01;	
	Vector2d b1 = Vector2d::Zero();
	b1 << -5.169911863249772e-01, 1.672384680188350e-01;
	// Nb meglio usare una reference?
	Err_rel_PALU(A1,b1,sol);
	
	Matrix2d A2 = Matrix2d::Zero();
	A2 <<5.547001962252291e-01, -5.540607316466765e-01, 8.320502943378437e-01, -8.324762492991313e-01;
	Vector2d b2 = Vector2d::Zero();
	b2 << -6.394645785530173e-04, 4.259549612877223e-04;
	// Nb meglio usare una reference?
	Err_rel_PALU(A2,b2,sol);
	
	Matrix2d A3 = Matrix2d::Zero();
	A3 << 5.547001962252291e-01, -5.547001955851905e-01, 8.320502943378437e-01,-8.320502947645361e-01;
	Vector2d b3 = Vector2d::Zero();
	b3 << -6.400391328043042e-10, 4.266924591433963e-10;
	// Nb meglio usare una reference?
	Err_rel_PALU(A3,b3,sol);
	
	
	
	
	
    return 0;
}
