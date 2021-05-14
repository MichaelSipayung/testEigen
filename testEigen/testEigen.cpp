#include <Eigen/Dense>
#include  <iostream>
#include <fstream>
#include <cmath>

typedef Eigen::Matrix<float, 100, 100> typedefForS;
int main()
{

	std::ofstream make("solveLinear.txt", std::ios::app);
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

	/*
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
	*/
	//vector
	Eigen::VectorXd v(3); //vector 1 dimension 
	//Eigen::VectorXd v1(1, 2); //2 dimension vector 
	/*
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

	Eigen::Matrix<double, 2, 3> matrixAddL;
	Eigen::Matrix<double, 3, 4> matrixAddR;
	matrixAddL(0, 0) = 1.0;
	matrixAddL(0, 1) = 2.0;
	matrixAddL(0, 2) = 4.0;
	matrixAddL(1, 0) = 2.0;
	matrixAddL(1, 1) = 6.0;
	matrixAddL(1, 2) = 0.0;

	matrixAddR(0, 0) = 4.0, matrixAddR(0, 1) = 1.0, matrixAddR(0, 2) = 4.0, matrixAddR(0, 3) = 3.0;
	matrixAddR(1, 0) = 0.0, matrixAddR(1, 1) = -1.0, matrixAddR(1, 2) = 3.0, matrixAddR(1, 3) = 1.0;
	matrixAddR(2, 0) = 2.0, matrixAddR(2, 1) = 7.0, matrixAddR(2, 2) = 5.0, matrixAddR(2, 3) = 2.0;

	std::cout << " show first matrix \t: [" << std::endl;
	std::cout << matrixAddL << std::endl;
	std::cout << "]" << std::endl;

	std::cout << " show second matrix \t: [" << std::endl;
	std::cout << matrixAddR << std::endl;
	std::cout << "]" << std::endl;

	Eigen::Matrix<double, 2, 4> multiplyMat;
	multiplyMat = matrixAddL * matrixAddR;
	std::cout << " show after multiply  matrix \t: [" << std::endl;
	std::cout << multiplyMat << std::endl;
	std::cout << "]" << std::endl;
	std::cout<<" And it's transpose \t: ["<<std::endl<<multiplyMat.transpose()<<std::endl;
	make << "matrix a \t:" << std::endl;
	make << matrixAddL<<std::endl;
	make << "matrix b \t:" << std::endl;
	make << matrixAddR << std::endl;
	make << "hasil" << std::endl;
	make << multiplyMat << std::endl;
	make << "then transpose \t: \n" << multiplyMat.transpose() << std::endl;
	Eigen::Matrix<double, 3, 3> invers;
	invers(0, 0) = 1; invers(0, 1) = 2; invers(0, 2) = 3;
	invers(1, 0) = 2; invers(1, 1) = 5; invers(1, 2) = 3;
	invers(2, 0) = 1; invers(2, 1) = 0; invers(2, 2) = 8;
	make << "\nMatriks Awal \t: \n" << invers << std::endl << "\nMatriks Invers \t: \n" << invers.inverse() << std::endl;

	

	make.close();

	make.open("solveEquation.txt", std::ios::app);
	make << "\nx\t+2y\t+3z\t=5\n";
	make << "2x\t+5y\t+3z\t=3\n";
	make << "x\t+0y\t+8z\t=17\n";

	//solve
	Eigen::Matrix<double, 3, 3> aMatriks;
	aMatriks(0, 0) = 1, aMatriks(0,1) = 2, aMatriks(0, 2) = 3;
	aMatriks(1, 0) = 2, aMatriks(1, 1) = 5, aMatriks(1, 2) = 3;
	aMatriks(2, 0) = 1, aMatriks(2, 1) = 0, aMatriks(2, 2) = 8;
	Eigen::Matrix<double,3,1> bMatriks;
	bMatriks(0, 0) = 5, bMatriks(1, 0) = 3, bMatriks(2, 0) = 17;
	make << "Hasil \t:\n";
	make << aMatriks.inverse() * bMatriks;



	make.close();
	*/
	
	std::cout << "Generate matrix 2  dimension\t:" << std::endl;
	Eigen::Matrix<double, 2, 2>  name2dimen;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dynamicSize;
	Eigen::MatrixXi m(10, 10);//generate arbitary mattrix
	m(0, 1) = 1;//set matrix value
	//std::cout << "Using Typedef \t: \n" << m << std::endl;
	Eigen::VectorXd vec(10);//generate 10x1 matrix
	vec[0] = 12;
	vec[1] = 13;
	//std::cout << "Test Vector \t: \n" << vec << std::endl;
	Eigen::Matrix3f ThreeTime;
	ThreeTime << 1, 2, 3, 1, 2, 3, 1, 2, 3;
	std::cout << " Initialization Matrix 2 dimension \t: \n" << ThreeTime << std::endl;
	Eigen::MatrixXd ten(3, 3);
	ten << 2, 3, 1, 11, 21, 32, 13, 41, 13;
	std::cout << "Show Matrix \t: \n" << ten << std::endl;
	Eigen::Matrix<double, 3, 3> findIn;
	findIn << 1, 2, 3,2,5,3,1,0,8;
	std::cout << "Show Matrix Element \t: [\n"<<findIn<<std::endl ;
	std::cout << "And  it's invers \t: [\n" << findIn.inverse() << std::endl;
	std::cout << " Find Determinant \t:[" << std::endl;
	Eigen::Matrix<double, 4, 4> findDet;
	findDet << 3, 5, -2, 6, 1, 2, -1, 1, 2, 4, 1, 5, 3, 7, 5, 3;
	std::cout << "Show Matrix \t:[" << std::endl << findDet << std::endl;
	std::cout << "Show Determinant \t:[\n" << findDet.determinant() << std::endl;
	std::cout << "Show norm matrix\t: " << findDet.norm() << std::endl;
	Eigen::Matrix<float, 3, 1> columnVe;
	columnVe << -3,2,1;
	std::cout << "Original Vector \t:[\n" << columnVe << std::endl;
	std::cout << "Find Norm Vector \t:[\n" << columnVe.norm()<<std::endl;
	std::cout << "Square  root 14 is  equal to norm vector \t:  [" << sqrt(14)<<"]" << std::endl;
	std::cout << "Dynamic matrix size  \t:[\n";
	Eigen::Matrix<float, 4, 1 > dynamicSi;
	dynamicSi << 2, -1, 3, -5;//column vector
	std::cout << "Test dynamic size , column vector \t: " << dynamicSi << std::endl <<
		"and it's  norm \t[\n" << dynamicSi.norm() << "]\n";
	std::cout << "test diagonal \t:[\n" << dynamicSi.diagonal() << "]" << std::endl;
	std::cout << "Matrix Persegi \t:[\n";
	Eigen::Matrix<float, 3, 3> diag;
	diag << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	std::cout << diag.diagonal() << "]" << std::endl;
	std::cout << "Test fixed size column or row " << std::endl;
	Eigen::Vector2d aFix(1, 2);//row vector 
	std::cout << "Test fixed size \t:[\n" << aFix << "]" << std::endl;
	//Eigen::Matrix<int, 5, 1> fixedSizeF{ 1,2,3,4,5 };//for c++-0x . if c++-11 is enabled;
	/*Eigen::MatrixX2d aDoub{
		{1,2},
		{3,4}
	};*/ //for c++11

	std::cout << "Assign Matrix \t:[\n";
	Eigen::Matrix2d doubleTwo(2, 2);
	doubleTwo(0, 0) = 3, doubleTwo(0, 1) = 0, doubleTwo(1, 0) = 8, doubleTwo(1,1) = -1;
	std::cout << "2 x 2 matrix assign \t:[\n"<<doubleTwo<<"]"<<std::endl;
	std::cout << "Find Eigen value \t: [\n" << doubleTwo.eigenvalues() << "]" << std::endl;
	std::cout << "Eigen value for 3x3 matrices \t:[" << std::endl;
	Eigen::Matrix<double, 3, 3>eigenVal;
	eigenVal << 0, 1, 0, 0, 0, 1, 4, -17, 8;
	std::cout << "Matrix element\t:[\n";
	std::cout << eigenVal << std::endl;
	std::cout << "Ant it's eigen value \t:[\n " << eigenVal.eigenvalues()<<"]" << std::endl;
	
	std::cout << "Resize Matrix \t:[\n";
	Eigen::MatrixXd matrixRes(3,2);
	std::cout << "Before Resizing Matrix, \n[Size Rows \t: " << matrixRes.rows() << "]" << "\n[Size Cols\t:" << matrixRes.cols() << "]" << std::endl;
	matrixRes.resize(4, 1);
	std::cout << "After Resizing Matrix, \n[Size Rows \t: " << matrixRes.rows() << "]" << "\n[Size Cols\t:" << matrixRes.cols() << "]" << std::endl;
	std::cout << "\nResize vector" << std::endl;
	Eigen::VectorXd vecRe(3, 1);
	std::cout << "Before Resizing vector, \n[Size Rows \t: " << vecRe.rows() << "]" << "\n[Size Cols\t:" << vecRe.cols() << "]" << std::endl;
	vecRe.resize(9, 1);
	std::cout << "After Resizing vector, \n[Size Rows \t: " << vecRe.rows() << "]" << "\n[Size Cols\t:" << vecRe.cols() << "]" << std::endl;
	





	make.close();

	return 0;
}

