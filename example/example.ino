#include <Matrix.h>

void setup() {
  Serial.begin(115200);
  Serial.println(":)");

  Matrix::Matrix A(10,10);
  for(int i=0; i<10; i++)
  {
    A.at(i,i,1);
  }
  A.add(A).print("A:");
  A.add(A).invert().print("Invert:");
  A.one();
  A.print("A:");
  

  Matrix::Matrix B(10);
  B.print("B:");
  for(int i=0; i<10; i++)
  {
    B.at(i,i,1);
  }
  B.print("B:");

  Matrix::one(3).print("one:");

  A.Block(2,3,Matrix::zero(3));
  A.print("A:");

  Matrix::Matrix C = Matrix::one(3);
  C.print("C:");

  (Matrix::one(3)+Matrix::zero(3)).print("test:");

  (Matrix::one(3)*Matrix::zero(3)).print("test:");

  (Matrix::one(3)-Matrix::eye(3)).print("test:");

  Matrix::Matrix D = Matrix::eye(3);
  D.swapRow(0,2);
  D.print("swapRow");
}

void loop() {
  // put your main code here, to run repeatedly:

}
