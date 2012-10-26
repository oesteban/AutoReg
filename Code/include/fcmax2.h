#ifndef FCMAX2__H
#define FCMAX2__H

//Versi�: 3/10/00
//fcmax.h q calcula el factor de comptes entre l'activada i la basal 
//ajustant una par�bola al perfil de quocients
//Optimitzaci� de cerca del m�xim (no falsos pics a l'extrem de l'histograma)

#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <stdlib.h>
#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkMultiplyByConstantImageFilter.h>
#include <itkCastImageFilter.h>
#include <vector>

//Define the vector types
typedef std::vector <double> VectorDouble;

//Defining image types
typedef double InternalPixel;

const unsigned int DIMENSION = 3;	
typedef itk::Image<InternalPixel, DIMENSION> InternalImageType;


double ComputeNormalizationFactor(InternalImageType::Pointer im_i, InternalImageType::Pointer im_ii, InternalImageType::Pointer mask);

void determinant(double a1,double a2,double a3,double a4,double a5,double a6,double a7,double a8,double a9,double *det);

#endif /* FCMAX2__H */
