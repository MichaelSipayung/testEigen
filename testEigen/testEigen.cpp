#include <Eigen/Dense>
#include  <iostream>

int main()
{
	 /*
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

	*/ 
	Eigen::Matrix<int, 2, 3> MatrixA;
	Eigen::Matrix<int, 3, 4> MatrixB;
	//begin 
	Eigen::MatrixXd  m = Eigen::MatrixXd::Random(3, 3); // matrixxd= matrix with type double and it reperesent arbitary size
	std::cout << "show matrix \t: ["<<std::endl;
	std::cout << m;
	std::cout << "]" << std::endl;
	std::cout << "modify the matrix \t: [" << std::endl;
	m = (m + Eigen::MatrixXd::Constant(3, 3, 1.2)) * 50;
	std::cout << "show matrix \t: [" << std::endl;
	std::cout << m;
	std::cout << "]" << std::endl;

	//vector
	Eigen::VectorXd v(3); //vector 1 dimension 
	//Eigen::VectorXd v1(1, 2); //2 dimension vector 

	v << 1, 2, 3;
	std::cout << " vector v value \t: " << std::endl << v << std::endl;
	std::cout << " Multiply matrix a * vector v \t: [" << std::endl;
	std::cout << m * v;
	std::cout << "]" << std::endl;

	std::cout << " Shwow matrix constant \t: [" << std::endl;
	std::cout << Eigen::MatrixXd::Constant(3, 3, 1.2) << std::endl;
	std::cout << "]" << std::endl;

	std::cout << " adding matrix constant \t: [" << std::endl<<std::endl;
	std::cout << " First Matrix \t: [" << std::endl;
	std::cout << Eigen::MatrixXd::Constant(3, 3, 1.0) << std::endl;
	std::cout << " Second Matrix \t: [" << std::endl;
	std::cout << Eigen::MatrixXd::Constant(3, 3, 2.0) << std::endl;
	std::cout << " Result Matrix \t: [" << std::endl;
	std::cout << Eigen::MatrixXd::Constant(3, 3, 1.0)+ Eigen::MatrixXd::Constant(3, 3, 2.0) << std::endl;
	std::cout << "]" << std::endl;
	
	//set at compile time
	Eigen::Matrix3d compileTime = Eigen::Matrix3d::Random();
	compileTime = (compileTime + Eigen::Matrix3d::Constant(1.2))*50;
	std::cout << " initialization at compile time \t: [" << std::endl;
	std::cout << compileTime;
	std::cout << std::endl;

	Eigen::Matrix<int, 2, 3> matrixAddL;
	Eigen::Matrix<int, 3, 4> matrixAddR;
	matrixAddL(0, 0) = 1;
	matrixAddL(0, 1) = 2;
	matrixAddL(0, 2) = 4;
	matrixAddL(1, 0) = 2;
	matrixAddL(1, 1) = 6;
	matrixAddL(1, 2) = 0;

	matrixAddR(0, 0) = 4, matrixAddR(0, 1) = 1, matrixAddR(0, 2) = 4, matrixAddR(0, 3) = 3;
	matrixAddR(1, 0) = 0, matrixAddR(1, 1) = -1, matrixAddR(1, 2) = 3, matrixAddR(1, 3) = 1;
	matrixAddR(2, 0) = 2, matrixAddR(2, 1) = 7, matrixAddR(2, 2) = 5, matrixAddR(2, 3) = 2;

	std::cout << " show first matrix \t: [" << std::endl;
	std::cout << matrixAddL << std::endl;
	std::cout << "]" << std::endl;

	std::cout << " show second matrix \t: [" << std::endl;
	std::cout << matrixAddR << std::endl;
	std::cout << "]" << std::endl;

	Eigen::Matrix<int, 2, 4> multiplyMat;
	multiplyMat = matrixAddL * matrixAddR;
	std::cout << " show after multiply  matrix \t: [" << std::endl;
	std::cout << multiplyMat << std::endl;
	std::cout << "]" << std::endl;
	std::cout<<" And it's transpose \t: ["<<std::endl<<multiplyMat.transpose()<<std::endl;







	
	return 0;
}