#include <iostream>
#include "Eigen/Eigen"
#include <iomanip>

using namespace std;
using namespace Eigen;




/*
	MatrixIsSingular(A) è la funzione implementata ad esercitazione per capire
	se la matrice A è invertibile o no.
	Ritorna vero se A è invertibile, altrimenti ritorna falso
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

/*
	Err_rel(A,b,sol) calcola l'errore relativo della soluziona del sistema lineare Ax=b
	rispetto alla soluzione di riferimento sol. In particolare:
	A matrice 2x2 di double
	b vettore 2x1 di double
	sol vettore 2x1 di double
	L'errore relativo viene calcolato nel caso in cui la matrice sia invertibile, per più motivi:
	 - E' garantita l'esistenza e l'unicità della soluzione esatta (dunque di riferimento) rispetto alla quale calcolare l'errore
	 - Da documentazione PartialPivLU non effettua questa verifica, perciò la impelmento prima di usarlo
*/
void Err_rel(const Matrix2d& A, const Vector2d& b, const Vector2d& sol){
	cout<<"Calcolo l'errore relativo della soluzione del sistema lineare Ax = b, usando decomposizione PALU e QR"<<endl;
	cout<<"Matrice A:"<<endl;
	cout <<setprecision(16) <<scientific << A << endl;
	cout<<"Vettore b:"<<endl;
	cout<<b<<endl;
	
	if (MatrixIsSingular(A)){
		Vector2d x_PALU=A.lu().solve(b);
		Vector2d x_QR=A.householderQr().solve(b);
		
		double er_rel = (sol-x_PALU).norm()/sol.norm();
		cout<<"Erorre relativo usando PALU: "<<er_rel<<endl;
		
		er_rel = (sol-x_QR).norm()/sol.norm();
		cout<<"Erorre relativo usando QR: "<<er_rel<<endl;
		
	}else{
		cout<<"La matrice relativa a questo sistema lineare non è invertibile"<<endl;
	}
}
/*
	Ho scelto di utilizzare la classe Eigen::PartialPivLU perchè nella consegna 
	dell'esercizio si fa riferimento alla fattorizzazione PALU, che è effettivamente quella 
	implementata da questa classe. La decomposizione effettuata dalla classe più generale 
	Eigen::FullPivLU è invece A=(P^−1)LU(Q^−1) e vale per ogni matrice. (P e Q di permutazione)
	
	Allo stesso modo la classe Eigen::HouseholderQR va ad effettuare proprio la fattorizzazione 
	A=QR richiesta. Le fattorizzazioni implementate da Eigen::FullPivHouseholderQR o 
	Eigen::ColPivHouseholderQR sono invece più artivolate.
*/


int main()
{
	Vector2d sol = -Vector2d::Ones();
	
	Matrix2d A1 = Matrix2d::Zero();
	A1 << 5.547001962252291e-01, -3.770900990025203e-02, 8.320502943378437e-01, -9.992887623566787e-01;	
	Vector2d b1 = Vector2d::Zero();
	b1 << -5.169911863249772e-01, 1.672384680188350e-01;
	
	Err_rel(A1,b1,sol);
	
	Matrix2d A2 = Matrix2d::Zero();
	A2 <<5.547001962252291e-01, -5.540607316466765e-01, 8.320502943378437e-01, -8.324762492991313e-01;
	Vector2d b2 = Vector2d::Zero();
	b2 << -6.394645785530173e-04, 4.259549612877223e-04;

	Err_rel(A2,b2,sol);
	
	Matrix2d A3 = Matrix2d::Zero();
	A3 << 5.547001962252291e-01, -5.547001955851905e-01, 8.320502943378437e-01,-8.320502947645361e-01;
	Vector2d b3 = Vector2d::Zero();
	b3 << -6.400391328043042e-10, 4.266924591433963e-10;

	Err_rel(A3,b3,sol);
	
    return 0;
}
