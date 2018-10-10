/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author:
 * @Date  : 03.03.2018
 */

#include <cmath>
#include <limits>
#include <functional>
#include <algorithm>

#include "Exception.hpp"
#include "Matrix.hpp"

#ifndef ANPI_LU_DOOLITTLE_HPP
#define ANPI_LU_DOOLITTLE_HPP

namespace anpi {


  /**
   * Auxiliary method used to debug LU decomposition.
   *
   * It separates a packed LU matrix into the lower triangular matrix
   * L and the upper triangular matrix U, such that the diagonal of L
   * is composed by 1's.
   */
  template<typename T>
  void unpackDoolittle(const Matrix<T>& LU,Matrix<T>& L,Matrix<T>& U) {
    //Ajustando los valores de L y U.
    L.allocate(LU.rows(),LU.cols());
    L.fill(T(0));
    U.allocate(LU.rows(),LU.cols());
    U.fill(T(0));
    //Fin del ajuste
    if (LU.cols()==LU.rows()){ //square Matrix (LU, L and U).
      for(unsigned int i=0;i<LU.cols();++i){
        for(unsigned int j=0;j<LU.cols();++j){
          if(i<j){ //Save on Upper matrix and store a zero on Lower matrix
            U[i][j] = LU[i][j];
            L[i][j] = T(0);
          }else{
            if(i>j){ //Save on Lower matrix and store a zero on Upper matrix
              U[i][j] = T(0);
              L[i][j] = LU[i][j];
            }else{//Diagonal. Store a '1' on Lower and the value on Upper
              U[i][j] = LU[i][j];
              L[i][j] = T(1);
            }
          }
        }
      }
    }
  }

  /**
   * Auxiliary method used to permute the matrix.
   *
   * Takes the place (column) and check under the diagonal
   * the mayor element to make the pivot.
   * Interchange the rows and copy the change in p.
   */
  template<typename T>
  void permutation(unsigned int col,Matrix<T>& LU,std::vector<size_t>& permut) {
    T mayorValue     = LU[col][col]; //Initial value (on diagonal).
    unsigned int row = col;
    for(unsigned int tempRow=col;tempRow<LU.rows();++tempRow){
      if (LU[tempRow][col] > mayorValue){
        mayorValue = LU[tempRow][col];
        row = tempRow;
      }
    }
    if (row != col){
      T temp = permut[row];
      permut[row] = permut[col];
      permut[col] = temp;
      for(unsigned int j=0;j<LU.cols();++j){
        temp       = LU[row][j];
        LU[row][j] = LU[col][j];
        LU[col][j] = temp;
      }
    }
  }

  /**
   * Decompose the matrix A into a lower triangular matrix L and an
   * upper triangular matrix U.  The matrices L and U are packed into
   * a single matrix LU.
   *
   * The L matrix will have in the Doolittle's LU decomposition a
   * diagonal of 1's
   *
   * @param[in] A a square matrix
   * @param[out] LU matrix encoding the L and U matrices
   * @param[out] permut permutation vector, holding the indices of the
   *             original matrix falling into the corresponding element.
   *             For example if permut[5]==3 holds, then the fifth row
   *             of the LU decomposition in fact is dealing with the third
   *             row of the original matrix.
   *
   * @throws anpi::Exception if matrix cannot be decomposed, or input
   *         matrix is not square.
   */
  template<typename T>
  void luDoolittle(const Matrix<T>& A, Matrix<T>& LU,
                   std::vector<size_t>& permut) {
    //Ajustando los valores de permut y LU.
    permut.resize(A.rows());
    for(unsigned int i=0;i<permut.size();++i){
      permut[i] = i;
    }
    LU.allocate(A.rows(),A.cols());
    LU.fill(T(0));
    //Fin del ajuste
    if ((A.rows()==permut.size()) && (A.cols()==A.rows())){
      LU = A;        //Se realiza la copia de la matriz A en LU
      for(unsigned int k=0;k<LU.cols()-1;++k){//Recorrido de columnas
        permutation(k,LU,permut);
        for(unsigned int i=k+1;i<LU.cols();++i){//Bajo la diagonal
          const T factor = LU[i][k]/LU[k][k];
          for(unsigned int j=k;j<LU.cols();++j){
            LU[i][j]-= factor*LU[k][j];
          }
          LU[i][k] = factor;
        }
      }
    }else{
      throw anpi::Exception("Dimensiones incompatibles\n");
    }
  }

  /*
  @brief : Metodo que calcula la matriz inversa a partir de la factorizacion LU.
  @param inv: Matriz donde se almacenara el resultado.
  @param LU : Matriz factorizada y obtenida de A.
  */
  template<typename T>
  void luinv(Matrix<T>& inv, Matrix<T>& LU){
    std::vector<T> b;
    std::vector<T> x;
    b.resize(LU.rows());
    x.resize(LU.rows());
    inv.allocate(LU.rows(),LU.cols());
    inv.fill(T(0));

    for(unsigned int i=0; i<LU.rows();++i){
      // Vector de términos independientes, columna i de la identidad
      for(unsigned int j=0;j<LU.rows();++j){
        b[j] = (j == i) ? 1 : 0;
      }
      //Sustitucion hacia adelante
      for(unsigned int j=1;j<LU.rows();++j){    //Recorre columna
        for(unsigned int k=0;k<j;++k){          //Recorre fila
          b[j]-=LU[j][k]*b[k];
        }
      }

      //Sustitucion hacia atras.
      x[LU.rows()-1] = b[LU.rows()-1]/LU[LU.rows()-1][LU.rows()-1];
      for(int j=LU.rows()-2;j>=0;--j){
        x[j] = b[j];
        for(int k=LU.rows()-1;k>j;--k){
          x[j]-=LU[j][k]*x[k];
        }
        x[j] = x[j]/LU[j][j];
      }
      //Almacenamiento
      for(unsigned int j=0;j<LU.rows();++j){
        inv[j][i] = x[j];
      }
    }
  }

}

#endif
