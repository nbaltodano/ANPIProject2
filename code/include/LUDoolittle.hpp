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
#include <iostream>       //cout
#include "Exception.hpp"
#include "Matrix.hpp"

// #include "Intrinsics.hpp"
#include <immintrin.h>
#include <type_traits>
#include <limits>               //Max value of unsigned int 32 bits (Mask on permut)

#ifndef ANPI_LU_DOOLITTLE_HPP
#define ANPI_LU_DOOLITTLE_HPP

namespace anpi{

  /*
@brief Metodo utilizado para imprimir una matriz
@param c: Matriz que se desea imprimir.
*/
template<class T>
void imprimirA(anpi::Matrix<T> c){
  std::cout << "\nImpresion sin dcols\n" << std::endl;
  for (unsigned int i=0;i<c.rows();i++){
    for (unsigned int j=0;j<c.cols();j++){
     std::cout << c[i][j] << "      ";
    }
    std::cout << "\n";
  }
}

template<class T>
void imprimirB(anpi::Matrix<T> c){
  std::cout << "\nImpresion con dcols\n" << std::endl;
  for (unsigned int i=0;i<c.rows();i++){
    for (unsigned int j=0;j<c.dcols();j++){
     std::cout << c[i][j] << "      ";
    }
    std::cout << "\n";
  }
}

/*
@brief Metodo utilizado para imprimir un vector transpuesto.
@param res Vector que se desea imprimir
*/
void imprimir(std::vector<size_t> res){
  for(unsigned int i=0;i<=res.size()-1;i++){
      std::cout << res[i] << "  ";
    }
    std::cout << "\n";
}



  namespace fallback {

    template<typename T>
    void solveLU(const anpi::Matrix<T>& A,std::vector<size_t>& x, const std::vector<size_t>& b);
    template<typename T>
    void unpackDoolittle(const Matrix<T>& LU,Matrix<T>& L,Matrix<T>& U);
    template<typename T>
    void permutation(unsigned int col,Matrix<T>& LU,std::vector<size_t>& permut);
    template<typename T>
    void luDoolittle(const Matrix<T>& A, Matrix<T>& LU,
                   std::vector<size_t>& permut);
    template<typename T>
    void luinv(Matrix<T>& inv, Matrix<T>& LU);


     /**
   * Solve a sistem of equations with QR decomposition.
   * @param[in]  A :matrix to decompose
   * @param[out] x :vector of solutions
   * @param[out] b :vector of results
   *
   */
  template<typename T>
  void solveLU(const anpi::Matrix<T>& A,std::vector<size_t>& x, const std::vector<size_t>& b){
    anpi::Matrix<T> LU;
    std::vector<size_t> permut;
    anpi::fallback::luDoolittle(A,LU,permut);
    std::vector<size_t> temp;
    temp.resize(b.size());
    x.resize(b.size());

    //Sustitucion hacia adelante
    for(unsigned int j=1;j<LU.rows();++j){    //Recorre columna
      for(unsigned int k=0;k<j;++k){          //Recorre fila
        temp[j]-=LU[j][k]*b[k];
      }
    }

    //Sustitucion hacia atras.
    x[LU.rows()-1] = temp[LU.rows()-1]/LU[LU.rows()-1][LU.rows()-1];
    for(int j=LU.rows()-2;j>=0;--j){
      x[j] = temp[j];
      for(int k=LU.rows()-1;k>j;--k){
        x[j]-=LU[j][k]*x[k];
      }
      x[j] = x[j]/LU[j][j];
    }
  }

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

  } // namespace fallback

  namespace simd {
  
  template<typename T>
  void solveLU(const anpi::Matrix<T>& A,std::vector<size_t>& x, const std::vector<size_t>& b);
  template<typename T>
  void unpackDoolittle(const Matrix<T>& LU,Matrix<T>& L,Matrix<T>& U);
  template<typename T>
  void permutation(unsigned int col,Matrix<T>& LU,std::vector<size_t>& permut);
  template<typename T>
  void luDoolittle(const Matrix<T>& A, Matrix<T>& LU,
                   std::vector<size_t>& permut);
  template<typename T>
  void mm_permut(T *a,T* b,unsigned int* j);
  template<typename T>
  int tamanoDePaso();
  template<typename T>
  void opOverRows(const T factor,T *a,T *b,unsigned int* i);

  #ifdef __AVX512F__
    template<>
    void mm_permut<double>(double* a, double* b,unsigned int* j) {
      __m128i atmp = _mm_load_epi64(a);                       //Load "a" into a __128i register
      __m128i btmp = _mm_load_epi64(b);                       //Load "b" into a __128i register
      _mm_store_epi64(a,btmp);
      _mm_store_epi64(b,atmp);
      *j+=2;
    }
    template<>
    void mm_permut<float>(float* a, float* b,unsigned int* j) {
      __m128i atmp = _mm_load_epi32(a);                       //Load "a" into a __128i register
      __m128i btmp = _mm_load_epi32(b);                       //Load "b" into a __128i register
      _mm_store_epi32(a,btmp);
      _mm_store_epi32(b,atmp);
      *j+=4; 
    }
    template<>
    void opOverRows<double>(const double factor,double *a,double *b,unsigned int* i){
      std::allocator<double> alloc;
      double* factors = alloc.allocate(2);
      factors[0] = factor; factors[1] = factor;
      uint8_t omask = std::numeric_limits<uint8_t>::max();      //Get a mask to multiply
      __mmask8 otherMask = _load_mask8 (&omask);                //Mask to multiply
      __m128d tmpA = _mm_maskz_loadu_pd (omask, a);             //Load "a"       into a __128i register
      __m128d tmpB = _mm_maskz_loadu_pd (omask, b);             //Load "b"       into a __128i register
      __m128d tmpF = _mm_maskz_loadu_pd (omask, factors);       //Load "factors" into a __128i register

      __m128d FB  = _mm_maskz_mul_pd (otherMask, tmpF,tmpB);    //Multiply Factors with "b".
      __m128d res = _mm_maskz_sub_pd (otherMask,tmpA,FB);       //Subtract A - Factors * b
     _mm_store_epi64(a,res);                                    //Store the value in memory

      alloc.deallocate(factors,2);                              //Free memory
      *i+=2;
    }
    template<>
    void opOverRows<float>(const float factor,float *a,float *b,unsigned int* i){
      std::allocator<float> alloc;
      float* factors = alloc.allocate(4);
      factors[0] = factor; factors[1] = factor; factors[2] = factor; factors[3] = factor;
      __m128i tmpA = _mm_load_epi32(a);                           //Load "a"       into a __128i register
      __m128i tmpB = _mm_load_epi32(b);                           //Load "b"       into a __128i register
      __m128i tmpF = _mm_load_epi32(factors);                     //Load "factors" into a __128i register
      
      __m128d FB  = _mm_maskz_madd_epi16 (otherMask, tmpF,tmpB);  //Multiply Factors with "b".
      __m128d res = _mm_maskz_sub_epi32 (otherMask,tmpA,FB);      //Subtract A - Factors * b
      _mm_store_epi32(a,res);                                     //Store the value in memory      

      alloc.deallocate(factors,4);                                //Free memory
      *i+=4;
    }
  /////////////////////////////////////////////////////////////////////
  #elif defined __AVX__
    template<>
    void mm_permut<double>(double* a, double* b,unsigned int* j) {
      __m256 tmpA = _mm256_load_pd (a);
      __m256 tmpB = _mm256_load_pd (b);

      _mm256_store_pd (a,tmpB);
      _mm256_store_pd (b,tmpA);
      
      *j+=4;
    }
    template<>
    void mm_permut<float>(float* a, float* b,unsigned int* j) {
      __m256 tmpA = _mm256_load_ps (a);
      __m256 tmpB = _mm256_load_ps (b);

      _mm256_store_ps (a,tmpB);
      _mm256_store_ps (b,tmpA);
      *j+=8;
    }
    template<>
    void opOverRows<double>(const double factor,double *a,double *b,unsigned int* i){
      std::allocator<double> alloc;
      double* factors = alloc.allocate(4);
      factors[0] = factor; factors[1] = factor; factors[2] = factor; factors[3] = factor;
      __m256 tmpA = _mm256_load_pd (a);                                    //Load the value of "a"       in a register of 256 bits
      __m256 tmpB = _mm256_load_pd (b);                                    //Load the value of "b"       in a register of 256 bits
      __m256 tmpF = _mm256_load_pd (factors);                              //Load the value of "factors" in a register of 256 bits     

      __m256 FB   = _mm256_mul_pd (tmpF,tmpB);                             //Multiply factors * b
      __m256 res  = _mm256_sub_pd (tmpA,FB);                               //Sub a-factors*b
      _mm256_store_ps (a,res);                                             //Store the result in "a"

      alloc.deallocate(factors,4);                                         //Free memory
      *i+=2;
    }
    template<>
    void opOverRows<float>(const float factor,float *a,float *b,unsigned int* i){
      std::allocator<float> alloc;
      float* factors = alloc.allocate(8);
      factors[0] = factor; factors[1] = factor; factors[2] = factor; factors[3] = factor;
      factors[4] = factor; factors[5] = factor; factors[6] = factor; factors[7] = factor;
      uint8_t omask = std::numeric_limits<uint8_t>::max();                 //Get a mask to multiply
      __m256 tmpA = _mm256_load_ps (a);                                    //Load the value of "a"       in a register of 256 bits
      __m256 tmpB = _mm256_load_ps (b);                                    //Load the value of "b"       in a register of 256 bits
      __m256 tmpF = _mm256_load_ps (factors);                              //Load the value of "factors" in a register of 256 bits     

      __m256 FB   = _mm256_mul_ps (tmpF,tmpB);                             //Multiply factors * b
      __m256 res  = _mm256_sub_ps (tmpA,FB);                               //Sub a-factors*b
      _mm256_store_ps (a,res);                                             //Store the result in "a"

      alloc.deallocate(factors,4);                                         //Free memory
      *i+=4;
    }

  /////////////////////////////////////////////////////////////////////    
  #elif  defined __SSE2__
    template<>
    void mm_permut<double>(double* a, double* b,unsigned int* j) {
      __m128d atmp = _mm_load_pd(a);                            //Load "a" into a __128d register
      __m128d btmp = _mm_load_pd(b);                            //Load "b" into a __128d register
      _mm_store_pd (a,btmp);                                    //Change the positions on memory (two doubles of 64 bits)
      _mm_store_pd (b,atmp);                                    //Change the positions on memory (two doubles of 64 bits)
      *j+=2;
    }
    template<>
    void mm_permut<float>(float* a, float* b,unsigned int* j) {
      __m128d atmp  = _mm_load_pd(reinterpret_cast<double*>(a));//Cast pointers to double and load into __m128d
      __m128d btmp  = _mm_load_pd(reinterpret_cast<double*>(b));//Cast pointers to double and load into __m128d
      __m128i atmp1 = _mm_castpd_si128(atmp);                   //Cast pointers to __128i from __m128d
      __m128i btmp1 = _mm_castpd_si128(btmp);                   //Cast pointers to __128i from __m128d
      __m128d atmp2 = _mm_castsi128_pd(atmp1);                  //Cast pointers to __128d from __m128i
      __m128d btmp2 = _mm_castsi128_pd(btmp1);                  //Cast pointers to __128d from __m128i
      _mm_store_pd (reinterpret_cast<double*>(a),btmp2);        //Change the positions on memory (4 floats of 64 bits)
      _mm_store_pd (reinterpret_cast<double*>(b),atmp2);        //Change the positions on memory (4 floats of 64 bits)
      *j+=4;
    }
    template<>
    void opOverRows<double>(const double factor,double *a,double *b,unsigned int* i){
      std::allocator<double> alloc;
      double* factors = alloc.allocate(2);
      factors[0] = factor; factors[1] = factor;
      __m128d tmpA = _mm_load_pd(a);                            //Load "a"       into a __128d register
      __m128d tmpB = _mm_load_pd(b);                            //Load "b"       into a __128d register
      __m128d tmpF = _mm_load_pd(factors);                      //Load "factors" into a __128d register

      __m128d FB  = _mm_mul_pd (tmpF,tmpB);                     //Multiply Factors with "b".
      __m128d res = _mm_sub_pd (tmpA,FB);                       //Subtract A - Factors * b
      _mm_store_pd (a,res);                                     //Store res in memory (two doubles of 64 bits)

      alloc.deallocate(factors,2);                              //Free memory
      *i+=2;
    }
    template<>
    void opOverRows<float>(const float factor,float *a,float *b,unsigned int* i){
      std::allocator<float> alloc;
      float* factors = alloc.allocate(4);
      factors[0] = factor; factors[1] = factor; factors[2] = factor; factors[3] = factor;
      __m128 tmpA = _mm_load_ps(a);                            //Load "a"       into a __128d register
      __m128 tmpB = _mm_load_ps(b);                            //Load "b"       into a __128d register
      __m128 tmpF = _mm_load_ps(factors);                      //Load "factors" into a __128d register

      __m128 FB  = _mm_mul_ps (tmpF,tmpB);                     //Multiply Factors with "b".
      __m128 res = _mm_sub_ps (tmpA,FB);                       //Subtract A - Factors * b
      _mm_store_ps (a,res);                                     //Store res in memory (two doubles of 64 bits)      

      alloc.deallocate(factors,4);                                        //Free memory
      *i+=4;
    }

  #else
    template<>
    void mm_permut<double>(double* a, double* b,unsigned int* j) {
      //First double
      double tmp = *a;
      *a = *b;
      *b = tmp;
      ++a;
      ++b;
      //Second double
      tmp = *a;
      *a = *b;
      *b = tmp;
      ++a;
      ++b;
      *j+=2;
    }
    template<>
    void mm_permut<float>(float* a, float* b,unsigned int* j) {
      //First double
      double tmp = *a; *a = *b; *b = tmp; ++a; ++b;
      //Second double
      tmp = *a; *a = *b; *b = tmp; ++a; ++b;
      //Third double
      tmp = *a; *a = *b; *b = tmp; ++a; ++b;
      //Fourth double
      tmp = *a; *a = *b; *b = tmp; ++a; ++b;

      *j+=4;
    }
    template<>
    void opOverRows<double>(const double factor,double *a,double *b,unsigned int* i){
      // First operation
      a-= factor*b; ++a; b++;
      // Second operation
      a-= factor*b; ++a; b++;

      *i+=2;
    }
    template<>
    void opOverRows<float>(const float factor,float *a,float *b,unsigned int* i){
      // First operation
      a-= factor*b; ++a; b++;
      // Second operation
      a-= factor*b; ++a; b++;
      // Third operation
      a-= factor*b; ++a; b++;
      // Fourth operation
      a-= factor*b; ++a; b++;

      *i+=4;
    }
  #endif

    #ifdef __AVX__
      /*
      Return the size step to operate with AVX512, AVX or SS2 
      according to the type that we are using
      */
      template<>
      int tamanoDePaso<double>() {
        return 4;
      }

      /*
      Return the size step to operate with AVX512, AVX or SS2 
      according to the type that we are using
      */
      template<>
      int tamanoDePaso<float>() {
        return 8;
      }
    #else
      /*
      Return the size step to operate with AVX512, AVX or SS2 
      according to the type that we are using
      */
      template<>
      int tamanoDePaso<double>() {
        return 2;
      }

      /*
      Return the size step to operate with AVX512, AVX or SS2 
      according to the type that we are using
      */
      template<>
      int tamanoDePaso<float>() {
        return 4;
      }
    #endif

    /**
     * Auxiliary method used to permute the matrix.
     *
     * Takes the place (column) and check under the diagonal
     * the mayor element to make the pivot.
     * Interchange the rows and copy the change in p.
     */
    template<typename T>
    void permutation(unsigned int col,Matrix<T>& LU,std::vector<size_t>& permut) {
      T mayorValue     = LU[col][col];          //Initial value (on diagonal).
      unsigned int row = col;                   //row = row to interchange with the row accord to col.
      for(unsigned int tempRow=col;tempRow<LU.rows();++tempRow){
        if (LU[tempRow][col] > mayorValue){     //Find a mayor value on the column
          mayorValue = LU[tempRow][col];        //Store the mayor value to check in other iteration.
          row = tempRow;                        //Store the number of the row.
        }
      }
      if (row != col){                          //In case that we found something
        T temp = permut[row];                   
        permut[row] = permut[col];              //Interchange the positions on the permut vector.
        permut[col] = temp;                     //Interchange the positions on the permut vector.
        T *init = &LU[row][0];
        T *end  = &LU[row][LU.dcols()-1];
        T *temporary = &LU[col][0];
        unsigned int j=0;
        for(;init<=end;){
          #ifdef __AVX512F__                    //Support to AVX512
            mm_permut(init,temporary,&j);    
          #elif  __AVX__                        //Support to AVX
            mm_permut(init,temporary,&j);
          #elif  __SSE2__                       //Support to SSE2
            mm_permut(init,temporary,&j);
          #endif
          init      = &LU[row][j];
          temporary = &LU[col][j];
        }
      }
    }

     /**
   * Solve a sistem of equations with QR decomposition.
   * @param[in]  A :matrix to decompose
   * @param[out] x :vector of solutions
   * @param[out] b :vector of results
   *
   */
  template<typename T>
  void solveLU(const anpi::Matrix<T>& A,std::vector<size_t>& x, const std::vector<size_t>& b){
    anpi::fallback::solveLU(A,x,b);   //Does not implement Inttrinsics 
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
      //Adjusting the values of permut
      permut.resize(A.rows());
      for(unsigned int i=0;i<permut.size();++i){
        permut[i] = i;
      }
      //End of the adjust
      if ((A.rows()==permut.size()) && (A.cols()==A.rows())){ //Square matriz and rows of A and size of vector are equals.
        LU = A;                                               //Copy of "A" in "LU".
        for(unsigned int k=0;k<LU.cols()-1;++k){              //Running on the cols
          permutation(k,LU,permut);                           //Permut the matrix to get a better pivot.
          for(unsigned int i=k+1;i<LU.cols();++i){            //Under the diagonal.
            const T factor = LU[i][k]/LU[k][k];               //Store the factor (pivot).
            for(unsigned int j=k;j<LU.cols();){
              if (j%anpi::simd::tamanoDePaso<T>()==0){        //Operation with AVX512 || AVX || SSE2
                anpi::simd::opOverRows<T>(factor,&LU[i][j],&LU[k][j],&j);
              }else{                                          //Normal operation over the matrix
                LU[i][j]-= factor*LU[k][j];
                ++j;
              }
            }
            LU[i][k] = factor;
          }
        }
      }else{
        throw anpi::Exception("Dimensiones incompatibles\n");
      }
    }

    /**
     * Auxiliary method used to debug LU decomposition.
     *
     * It separates a packed LU matrix into the lower triangular matrix
     * L and the upper triangular matrix U, such that the diagonal of L
     * is composed by 1's.
     */
    template<typename T>
    void unpackDoolittle(const Matrix<T>& LU,Matrix<T>& L,Matrix<T>& U) {
      anpi::fallback::unpackDoolittle(LU,L,U);
    }
  

  } // namespace simd

  // The arithmetic implementation (aimpl) namespace
  // dispatches to the corresponding methods
#ifdef ANPI_ENABLE_SIMD
  namespace aimpl=simd;
#else
  namespace aimpl=fallback;
#endif
} // namespace anpi
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
