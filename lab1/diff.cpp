#include<fstream>
#include<iostream>
#include<array>
#include<vector>
#include<cmath>
// using namespace std;

template<typename RealType, unsigned int N>
struct DerivativeCoef
{
  RealType centralCoef;
  std::array<RealType, N> otherCoeffs;
};

template<typename RealType, unsigned int N>
DerivativeCoef<RealType, N>
calcDerivativeCoef(const std::array<RealType, N>& points) noexcept
{
    std::vector<std::vector<RealType>> A(N + 1, std::vector<RealType>(N + 1, 0));
    for (unsigned int i = 0; i <= N; i++){
        A[0][i] = 1;
    }
    for (unsigned int i = 1; i <= N; i++){
        A[i][0] = 0;
        for (unsigned int j = 1; j <= N; j++){
        A[i][j] = A[i - 1][j] * points[j - 1] / i;
        }
    }

    std::vector<RealType> b(N + 1, 0);
    b[1] = 1;

    std::vector<RealType> x(N + 1, 0);
    for (unsigned int i = 0; i <= N; i++){
        for (unsigned int j = i + 1; j <= N; j++){
            RealType factor = A[j][i] / A[i][i];
            for (unsigned int k = i; k <= N; k++){
                A[j][k] -= factor * A[i][k];
            }
            b[j] -= factor * b[i];
        }
    }     
    for (int i = N; i >= 0; i--){
        RealType sum = b[i];
        for (unsigned int j = i + 1; j <= N; j++){
            sum -= A[i][j] * x[j];
        }
        x[i] = sum / A[i][i];
    }

    RealType centralCoef = x[0];

    std::array<RealType, N> coefs;
    for (unsigned int i = 0; i < N; i++){
        coefs[i] = x[i + 1];
    }

    return DerivativeCoef<RealType, N>{centralCoef, coefs};
}

template<typename RealType, unsigned int N>
RealType dif(const RealType x0, const RealType h, const std::array<RealType, N>& points)
{
    DerivativeCoef<double, N> coefs = calcDerivativeCoef<RealType, N>(points);

    RealType dif_x0 = coefs.centralCoef * exp(x0) / h;
    for (unsigned int i = 0; i < N; i++){
    dif_x0 += coefs.otherCoeffs[i] * exp(x0 + points[i] * h) / h;
    }

    return dif_x0;
}


int main(){
	const double x0 = 1;
	const std::array< double, 16> step = { 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15 };

	std::ofstream data2("data2points.txt");
	data2.precision(15);

	for (int i = 0; i < 16; i++){
		data2 << step[i] << "	" << dif<double, 2>(x0, step[i], { -1, 1 }) << std::endl;
	}

	data2.close();


	std::ofstream data3("data3points.txt");
	data3.precision(15);

	for (int i = 0; i < 16; i++)
	{ 
		data3 << step[i] << "	" << log(std::abs(dif<double, 3>(x0, step[i], { -1, 1, 2 }) - exp(x0)))<< std::endl;
	}

	data3.close();
    
	std::ofstream data4("data4points.txt");
	data4.precision(15);

	for (int i = 0; i < 16; i++)
	{ 
		data4 << step[i] << "	" << log(std::abs(dif<double, 4>(x0, step[i], { -2, -1, 1, 2 })-exp(x0))) << std::endl;
	}

	data4.close();

	std::ofstream data5("data5points.txt");
	data5.precision(15);

	for (int i = 0; i < 16; i++)
	{ 
		data5 << step[i] << "	" << log(std::abs(dif<double, 5>(x0, step[i], { -2, -1, 1, 2, 3 })- exp(x0))) << std::endl;
	}

	data5.close();

	return 0;
}