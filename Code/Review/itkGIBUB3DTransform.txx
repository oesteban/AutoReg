/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGIBUB3DTransform.txx,v $
  Language:  C++
  Date:      $Date: 2010-03-30 15:20:02 $
  Version:   $Revision: 1.16 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGIBUB3DTransform_txx
#define __itkGIBUB3DTransform_txx

#include "itkGIBUB3DTransform.h"


namespace itk
{

// Constructor with default arguments
template<class TScalarType>
GIBUB3DTransform<TScalarType>::
GIBUB3DTransform() :
  Superclass(OutputSpaceDimension, ParametersDimension)
{
  m_isComputedZYX = false;
  m_AngleFI = m_AngleTheta = m_AngleZeta = NumericTraits<ScalarType>::Zero;
  
  m_Directions.Fill( 0.0 );
  m_Directions(0,1) = 1; m_Directions(1,0) = 1; m_Directions(2,2) = 1;

  m_spacing.Fill( 0.0 );
  m_spacing(0,0) = 1.0; m_spacing(1,1)= 1.0; m_spacing(2,2) = 1.0;
  
  m_RefSpacing.Fill( 0.0 );
  m_RefSpacing(0,0) = 1.0; m_RefSpacing(1,1)= 1.0; m_RefSpacing(2,2) = 1.0;
  
  m_SpacingDirection.Fill( 1.0 );
  m_SpacingDirection[1] = -1.0;
}


// Constructor with default arguments
template<class TScalarType>
GIBUB3DTransform<TScalarType>::
GIBUB3DTransform(unsigned int spaceDimension, 
                         unsigned int parametersDimension) :
  Superclass(spaceDimension,parametersDimension)
{
  m_isComputedZYX = false;
  m_AngleFI = m_AngleTheta = m_AngleZeta = NumericTraits<ScalarType>::Zero;
  
  m_Directions.Fill( 0.0 );
  m_Directions(0,1) = 1; m_Directions(1,0) = 1; m_Directions(2,2) = 1;

  m_spacing.Fill( 0.0 );
  m_spacing(0,0) = 1.0; m_spacing(1,1)= 1.0; m_spacing(2,2) = 1.0;
  
  m_RefSpacing.Fill( 0.0 );
  m_RefSpacing(0,0) = 1.0; m_RefSpacing(1,1)= 1.0; m_RefSpacing(2,2) = 1.0;
  
  m_SpacingDirection.Fill( 1.0 );
  m_SpacingDirection[1] = -1.0;
}

// Constructor with default arguments
template<class TScalarType>
GIBUB3DTransform<TScalarType>::
GIBUB3DTransform(const MatrixType & matrix,
                         const OutputPointType & offset) :
  Superclass(matrix, offset)
{
  m_isComputedZYX = true;
  
  m_Directions.Fill( 0.0 );
  m_Directions(0,1) = 1; m_Directions(1,0) = 1; m_Directions(2,2) = 1;
  

  m_spacing.Fill( 0.0 );
  m_spacing(0,0) = 1.0; m_spacing(1,1)= 1.0; m_spacing(2,2) = 1.0;
  
  m_RefSpacing.Fill( 0.0 );
  m_RefSpacing(0,0) = 1.0; m_RefSpacing(1,1)= 1.0; m_RefSpacing(2,2) = 1.0;
  
  m_SpacingDirection.Fill( 1.0 );
  m_SpacingDirection[1] = -1.0;
}


// Destructor
template<class TScalarType>
GIBUB3DTransform<TScalarType>::
~GIBUB3DTransform()
{
}



// Set Rotational Part
template <class TScalarType>
void
GIBUB3DTransform<TScalarType>
::SetGibUbRotation(ScalarType angleFI,ScalarType angleTheta,ScalarType angleZeta)
{
  m_AngleFI    = m_SpacingDirection[0] * angleZeta;
  m_AngleTheta = m_SpacingDirection[1] * angleTheta;
  m_AngleZeta  = m_SpacingDirection[2] * angleFI;
  
  this->ComputeGibUbMatrix();
  this->ComputeGibUbOffset();
}

// Set Angles
template <class TScalarType>
void
GIBUB3DTransform<TScalarType>
::SetVarGibUbRotation(ScalarType angleFI,ScalarType angleTheta,ScalarType angleZeta)
{
  m_AngleFI    = m_SpacingDirection[0] * angleZeta;
  m_AngleTheta = m_SpacingDirection[1] * angleTheta;
  m_AngleZeta  = m_SpacingDirection[2] * angleFI;
}

// Compute angles from the rotation matrix
template <class TScalarType>
void
GIBUB3DTransform<TScalarType>
::ComputeMatrixGibUbParameters(void)
{
//   if(m_ComputeZYX)
//     {
//     m_AngleY = -vcl_asin(this->GetMatrix()[2][0]);
//     double C = vcl_cos(m_AngleY);
//     if(vcl_fabs(C)>0.00005)
//       {
//       double x = this->GetMatrix()[2][2] / C;
//       double y = this->GetMatrix()[2][1] / C;
//       m_AngleX = vcl_atan2(y,x);
//       x = this->GetMatrix()[0][0] / C;
//       y = this->GetMatrix()[1][0] / C;
//       m_AngleZ = vcl_atan2(y,x);
//       }
//     else
//       {
//       m_AngleX = NumericTraits< ScalarType >::Zero;
//       double x = this->GetMatrix()[1][1];
//       double y = -this->GetMatrix()[0][1];
//       m_AngleZ = vcl_atan2(y,x);
//       }
//     }
//   else
//     {
//     m_AngleX = vcl_asin(this->GetMatrix()[2][1]);
//     double A = vcl_cos(m_AngleX);
//     if(vcl_fabs(A)>0.00005)
//       {
//       double x = this->GetMatrix()[2][2] / A;
//       double y = -this->GetMatrix()[2][0] / A;
//       m_AngleY = vcl_atan2(y,x);
// 
//       x = this->GetMatrix()[1][1] / A;
//       y = -this->GetMatrix()[0][1] / A;
//       m_AngleZ = vcl_atan2(y,x);
//       }
//     else
//       {
//       m_AngleZ = NumericTraits< ScalarType >::Zero;
//       double x = this->GetMatrix()[0][0];
//       double y = this->GetMatrix()[1][0];
//       m_AngleY = vcl_atan2(y,x);
//       }
//     }
  this->ComputeGibUbMatrix();
}


// Compute the matrix
template <class TScalarType>
void
GIBUB3DTransform<TScalarType>
::ComputeGibUbMatrix( void )
{
   const ScalarType c1 = vcl_cos( m_AngleZeta );
   const ScalarType c2 = vcl_cos( m_AngleTheta );
   const ScalarType c3 = vcl_cos( m_AngleFI );
   
   const ScalarType s1 = vcl_sin( m_AngleZeta );
   const ScalarType s2 = vcl_sin( m_AngleTheta );
   const ScalarType s3 = vcl_sin( m_AngleFI );
//   const ScalarType one = NumericTraits< ScalarType >::One;
//   const ScalarType zero = NumericTraits< ScalarType >::Zero;

  // Row 0
  m_GibUbMatrix(0,0) = c1*c3 - c2*s1*s3;
  m_GibUbMatrix(0,1) = -c3*s1 - c1*c2*s3;
  m_GibUbMatrix(0,2) = s2 * s3;

  // Row 1
  m_GibUbMatrix(1,0) = c2*c3*s1 + c1*s3;
  m_GibUbMatrix(1,1) = c1*c2*c3 - s1*s3;
  m_GibUbMatrix(1,2) = -c3 * s2;

  // Row 2
  m_GibUbMatrix(2,0) = s1 * s2;
  m_GibUbMatrix(2,1) = c1 * s2;
  m_GibUbMatrix(2,2) = c2;
  
  MatrixType invSpacing;
  invSpacing.Fill( 0.0 );
  
  invSpacing(0,0) = 1.0 / m_spacing(0,0);
  invSpacing(1,1) = 1.0 / m_spacing(1,1);
  invSpacing(2,2) = 1.0 / m_spacing(2,2);
  
  MatrixType newMatrix = m_spacing * ((m_Directions * m_GibUbMatrix * invSpacing) * m_Directions );
  
  this->SetVarMatrix( newMatrix );
}

// Compute the offset
template <class TScalarType>
void
GIBUB3DTransform<TScalarType>
::ComputeGibUbOffset( void )
{
  const MatrixType & matrix = this->GetGibUbMatrix();
  
  OffsetType offset;
  for(unsigned int i=0; i<itkGetStaticConstMacro(OutputSpaceDimension); i++)
    {
    offset[i] = m_GibUbTranslation[i] + m_GibUbCenter[i];
    for(unsigned int j=0; j<itkGetStaticConstMacro(InputSpaceDimension); j++)
      {
      offset[i] -= matrix[i][j] * m_GibUbCenter[j];
      }
    }

  m_GibUbOffset = offset;
  
  OutputVectorType directionOffset = m_RefSpacing * (m_Directions * m_GibUbOffset);
  OutputVectorType spacingOffset;
  spacingOffset[0] =  m_SpacingDirection[0] * directionOffset[0];
  spacingOffset[1] =  m_SpacingDirection[1] * directionOffset[1];
  spacingOffset[2] =  m_SpacingDirection[2] * directionOffset[2];
  
  this->SetOffset( spacingOffset );
}

template <class TScalarType>
void
GIBUB3DTransform<TScalarType>
::SetGibUbTranslation(const OutputVectorType & translation)
{
  m_GibUbTranslation = translation;
//  this->ComputeGibUbOffset();
  OutputVectorType directionTraslation = m_RefSpacing * (m_Directions * m_GibUbTranslation);
  OutputVectorType spacingTraslation;
  spacingTraslation[0] = m_SpacingDirection[0] * directionTraslation[0];
  spacingTraslation[1] = m_SpacingDirection[1] * directionTraslation[1];
  spacingTraslation[2] = m_SpacingDirection[2] * directionTraslation[2];
  this->SetTranslation( spacingTraslation );
}

template <class TScalarType>
void
GIBUB3DTransform<TScalarType>
::SetGibUbCenter(const InputPointType & center)
{
  m_GibUbCenter = center;
 // this->ComputeGibUbOffset();
  
  InputPointType directionCenter = m_spacing * (m_Directions * m_GibUbCenter);
  InputPointType spacingCenter;
  spacingCenter[0] = m_SpacingDirection[0] * directionCenter[0];
  spacingCenter[1] = m_SpacingDirection[1] * directionCenter[1];
  spacingCenter[2] = m_SpacingDirection[2] * directionCenter[2];
  
  this->SetCenter(spacingCenter);
}


template <class TScalarType>
void
GIBUB3DTransform<TScalarType>
::SetParameters( const ParametersType & parameters )
{
  Superclass::SetParameters( parameters );
}

// Set Parameters
template <class TScalarType>
void
GIBUB3DTransform<TScalarType>
::SetGibUbParameters( const GibUbParametersType & parameters )
{
  itkDebugMacro( << "Setting GIBUB parameters " << parameters );

  // Set angles with parameters
  m_AngleFI    = m_SpacingDirection[0] * parameters[2];
  m_AngleTheta = m_SpacingDirection[1] * parameters[1];
  m_AngleZeta  = m_SpacingDirection[2] * parameters[0];
    
  // Transfer the translation part
  OutputVectorType newTranslation;
  newTranslation[0] = parameters[3];
  newTranslation[1] = parameters[4];
  newTranslation[2] = parameters[5];

  m_spacing(0,0) = parameters[6];
  m_spacing(1,1) = parameters[7];
  m_spacing(2,2) = parameters[8];  
  
  this->ComputeGibUbMatrix();
  
  this->SetVarGibUbTranslation(newTranslation);
  
  this->ComputeGibUbOffset();
  

  // Modified is always called since we just have a pointer to the
  // parameters and cannot know if the parameters have changed.
  this->Modified();

  itkDebugMacro(<<"After setting parameters ");
}



// Get Parameters
template <class TScalarType>
const typename GIBUB3DTransform<TScalarType>::GibUbParametersType &
GIBUB3DTransform<TScalarType>
::GetGibUbParameters( void ) const
{
  this->m_GibUbParameters[0] = m_SpacingDirection[2] * m_AngleZeta;
  this->m_GibUbParameters[1] = m_SpacingDirection[1] * m_AngleTheta;
  this->m_GibUbParameters[2] = m_SpacingDirection[0] * m_AngleFI;
  this->m_GibUbParameters[3] = m_GibUbTranslation[0];
  this->m_GibUbParameters[4] = m_GibUbTranslation[1];
  this->m_GibUbParameters[5] = m_GibUbTranslation[2];
  this->m_GibUbParameters[6] = m_spacing(0,0);
  this->m_GibUbParameters[7] = m_spacing(1,1);
  this->m_GibUbParameters[8] = m_spacing(2,2);

  return this->m_GibUbParameters;
}

template <class TScalarType>
const typename GIBUB3DTransform<TScalarType>::ParametersType &
GIBUB3DTransform<TScalarType>
::GetParameters( void ) const
{
  return this->Superclass::GetParameters();
}


// Get jacobian
template<class TScalarType>
const typename GIBUB3DTransform<TScalarType>::JacobianType &
GIBUB3DTransform<TScalarType>::
GetJacobian( const InputPointType & p ) const
{
  return this->Superclass::GetJacobian( p );
}

// Get an inverse of this transform
template<class TScalarType>
bool
GIBUB3DTransform<TScalarType>
::GetInverse(Self* inverse) const
{
  return this->Superclass::GetInverse(inverse);
}
   
// Return an inverse of this transform
template<class TScalarType>
typename GIBUB3DTransform<TScalarType>::InverseTransformBasePointer
GIBUB3DTransform<TScalarType>
::GetInverseTransform() const
{
  Pointer inv = New();
  return this->GetInverse(inv) ? inv.GetPointer() : NULL;
}


// Print self
template<class TScalarType>
void
GIBUB3DTransform<TScalarType>::
PrintSelf(std::ostream &os, Indent indent) const
{

  Superclass::PrintSelf(os,indent);
  
  os << indent << "Euler's GibUb angles: AngleFi=" << m_AngleFI
     << " AngleTheta=" << m_AngleTheta
     << " AngleZeta=" << m_AngleZeta
     << std::endl;
  os << indent << "m_isComputedZYX = " << m_isComputedZYX << std::endl;
}


} // namespace

#endif
