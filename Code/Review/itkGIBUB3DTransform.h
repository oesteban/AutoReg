/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGIBUB3DTransform.h,v $
  Language:  C++
  Date:      $Date: 2010-03-16 08:49:50 $
  Version:   $Revision: 1.15 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGIBUB3DTransform_h
#define __itkGIBUB3DTransform_h

#include <iostream>
#include "itkEuler3DTransform.h"
#include "itkExceptionObject.h"
#include "itkMatrix.h"
#include "itkVersor.h"
#include "itkVector.h"


namespace itk
{

/** \class GIBUB3DTransform
 * \brief GIBUB3DTransform of a vector space (e.g. space coordinates)
 *
 * This transform applies a rotation about a specific coordinate or
 * centre of rotation followed by a translation.
 *
 * \ingroup Transforms
 */
template < class TScalarType=double >    // Data type for scalars 
class ITK_EXPORT GIBUB3DTransform : 
        public Euler3DTransform< TScalarType >
{
public:
  /** Standard class typedefs. */
  typedef GIBUB3DTransform        Self;
  typedef Euler3DTransform< TScalarType > Superclass;
  typedef SmartPointer<Self>              Pointer;
  typedef SmartPointer<const Self>        ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro( GIBUB3DTransform, Euler3DTransform );

  /** New macro for creation of through a Smart Pointer */
  itkNewMacro( Self );

  /** Dimension of the space. */
  itkStaticConstMacro(SpaceDimension, unsigned int, 3);
  itkStaticConstMacro(InputSpaceDimension, unsigned int, 3);
  itkStaticConstMacro(OutputSpaceDimension, unsigned int, 3);
  itkStaticConstMacro(ParametersDimension, unsigned int, 6);
  itkStaticConstMacro(GIBUBParametersDimension, unsigned int, 9 );

  typedef typename Superclass::ParametersType             ParametersType;
  typedef typename Superclass::ParametersType             GibUbParametersType;
  typedef typename Superclass::ParametersValueType        ParametersValueType;
  typedef typename Superclass::ParametersValueType        GibUbParametersValueType;
  typedef typename Superclass::JacobianType               JacobianType;
  typedef typename Superclass::ScalarType                 ScalarType;
  typedef typename Superclass::InputVectorType            InputVectorType;
  
  typedef typename Superclass::OutputVectorType           OutputVectorType;
  typedef typename Superclass::InputCovariantVectorType  
                                                     InputCovariantVectorType;
  typedef typename Superclass::OutputCovariantVectorType  
                                                     OutputCovariantVectorType;

  typedef typename Superclass::InputVnlVectorType         InputVnlVectorType;
  typedef typename Superclass::OutputVnlVectorType        OutputVnlVectorType;
  typedef typename Superclass::InputPointType             InputPointType;
  typedef typename Superclass::OutputPointType            OutputPointType;
  typedef typename Superclass::MatrixType                 MatrixType;
  typedef typename Superclass::InverseMatrixType          InverseMatrixType;
  typedef typename Superclass::CenterType                 CenterType;
  typedef typename Superclass::TranslationType            TranslationType;
  typedef typename Superclass::TranslationValueType       TranslationValueType;
  typedef typename Superclass::OffsetType                 OffsetType;

  
  typedef itk::Vector< TScalarType, itkGetStaticConstMacro(InputSpaceDimension) > SpacingType;
  /** Base inverse transform type. This type should not be changed to the
   * concrete inverse transform type or inheritance would be lost. */
  typedef typename Superclass::InverseTransformBaseType InverseTransformBaseType;
  typedef typename InverseTransformBaseType::Pointer    InverseTransformBasePointer;

  /** Set the transformation from a container of parameters
   * This is typically used by optimizers.  There are nine parameters. The first
   * three represent the angles of rotation (in radians) around each one of the
   * axes (X,Y,Z), the next three parameters represent the coordinates of the
   * center of rotation and the last three parameters represent the
   * translation. */
  void SetParameters( const ParametersType & parameters );
  
  void SetGibUbParameters( const GibUbParametersType & parameters );
  
  /** Get the parameters that uniquely define the transform
   * This is typically used by optimizers. There are nine parameters. The first
   * three represent the angles of rotation (in radians) around each one of the
   * axes (X,Y,Z), the next three parameters represent the coordinates of the
   * center of rotation and the last three parameters represent the
   * translation. */
  const ParametersType & GetParameters( void ) const;
  
  const GibUbParametersType & GetGibUbParameters( void ) const;

  /** This method computes the Jacobian matrix of the transformation.
   * given point or vector, returning the transformed point or
   * vector. The rank of the Jacobian will also indicate if the 
   * transform is invertible at this point. */
  const JacobianType & GetJacobian(const InputPointType  &point ) const;

  /** Get an inverse of this transform. */
  bool GetInverse(Self* inverse) const;

  /** Return an inverse of this transform. */
  virtual InverseTransformBasePointer GetInverseTransform() const;
  
  
  void ComputeGibUbMatrix( void );
  void ComputeMatrixGibUbParameters( void );
  void ComputeGibUbTranslation(void);
  
  void SetVarGibUbRotation(ScalarType angleFI,ScalarType angleTheta,ScalarType angleZeta);
  void SetGibUbRotation(ScalarType angleFI,ScalarType angleTheta,ScalarType angleZeta);
  
  void SetVarGibUbTranslation(const OutputVectorType & translation)
    { m_GibUbTranslation = translation; }
    
  void SetGibUbTranslation(const OutputVectorType & translation);
  
  void ComputeGibUbOffset(void);
  void SetVarGibUbOffset(const OutputVectorType & offset)
    { m_GibUbOffset = offset; }
    
  const MatrixType & GetGibUbMatrix() const
    { return m_GibUbMatrix; }
  
  void SetVarGibUbCenter(const InputPointType & center)
    {
    m_GibUbCenter = center;
    }  
  void SetGibUbCenter(const InputPointType & center);

  void SetImageSpacing( const SpacingType & spacing )
   { m_spacing(0,0) = spacing[0]; m_spacing(1,1) = spacing[1];m_spacing(2,2) = spacing[2]; }
   
  void SetReferenceSpacing( const SpacingType & spacing )
   { m_RefSpacing(0,0) = spacing[0]; m_RefSpacing(1,1) = spacing[1];m_RefSpacing(2,2) = spacing[2]; }
  
protected:
  GIBUB3DTransform();
  GIBUB3DTransform(unsigned int SpaceDimension,
                           unsigned int ParametersDimension);
  GIBUB3DTransform(const MatrixType & matrix,
                           const OutputPointType & offset);
  ~GIBUB3DTransform();

  /**
   * Print contents of an GIBUB3DTransform
   */
  void PrintSelf(std::ostream &os, Indent indent) const;

private:
  GIBUB3DTransform(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  ScalarType                     m_AngleFI; 
  ScalarType                     m_AngleTheta; 
  ScalarType                     m_AngleZeta;

  InputPointType                 m_GibUbCenter;
  OutputVectorType               m_GibUbTranslation;
  OutputVectorType               m_GibUbOffset;

  MatrixType                     m_spacing;
  MatrixType                     m_RefSpacing;
  
  MatrixType                     m_Directions;
  MatrixType                     m_GibUbMatrix;
  
  OutputVectorType               m_SpacingDirection;
  
  mutable  GibUbParametersType   m_GibUbParameters;
  
  bool        m_isComputedZYX;
  
}; //class GIBUB3DTransform


}  // namespace itk

/** Define instantiation macro for this template. */
#define ITK_TEMPLATE_GIBUB3DTransform(_, EXPORT, x, y) namespace itk { \
  _(1(class EXPORT GIBUB3DTransform< ITK_TEMPLATE_1 x >)) \
  namespace Templates { typedef GIBUB3DTransform< ITK_TEMPLATE_1 x > \
                                            GIBUB3DTransform##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkGIBUB3DTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkGIBUB3DTransform.txx"
#endif

#endif /* __itkGIBUB3DTransform_h */
