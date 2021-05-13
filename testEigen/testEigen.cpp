#include <Eigen/Dense>
#include  <iostream>

int main()
{
	Eigen::Matrix<double, 10, 10> A;
	A.setZero();
	//A(9, 0) = 1.234;
	for (size_t i = 0; i < 10; ++i)
	{
		for (size_t j = 0; j < 10; ++j) {
			A(i, j) = i * j;
		}
	}
	Eigen::Matrix<double, 10, 10> B;
	B.setZero();
	//B(9, 0) = 1.234;
	for (size_t i = 0; i < 10; ++i)
	{
		for (size_t j = 0; j < 10; ++j) {
			B(i, j) = i * j;
		}
	}
	
	std::cout << " Matrix A \t: " << std::endl;
	std::cout << A << std::endl;
	std::cout << " Matrix B \t: " << std::endl;
	std::cout << B << std::endl;

	std::cout << " result multiplication \t: " << std::endl;
	Eigen::Matrix<double, 10, 10> C;
	C.setZero();
	C = A * B;
	std::cout << "\n Hasil Kali Matrix a * b \t: " << C << std:: endl;;


	return 0;
}