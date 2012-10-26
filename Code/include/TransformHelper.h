/*
    Copyright (c) 2010, Oscar Esteban - oesteban@die.upm.es
                        with Biomedical Image Technology, UPM (BIT-UPM)
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
        * Neither the name of the BIT-UPM, nor the names of its contributors
        may be used to endorse or promote products derived from this software
        without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY Oscar Esteban ''AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL OSCAR ESTEBAN BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TRANSFORMHELPER_H
#define TRANSFORMHELPER_H

#ifndef DEG
#define DEG M_PI/180.0
#endif

#include <boost/any.hpp>
#include <boost/type_traits.hpp>
#include <vector>
#include <stdio.h>
#include <iostream>

#include <itkObject.h>
#include <itkSmartPointer.h>
#include <itkArray.h>
#include <itkTransform.h>
#include <itkMatrixOffsetTransformBase.h>
#include <itkCenteredEuler3DTransform.h>
#include <itkEuler3DTransform.h>
#include <itkVersorRigid3DTransform.h>
#include <itkVersorTransform.h>
#include <itkScaleSkewVersor3DTransform.h>
#include <itkScaleVersor3DTransform.h>
#include <itkSimilarity3DTransform.h>
#include <itkAffineTransform.h>
#include <itkCenteredAffineTransform.h>
#include <itkBSplineDeformableTransform.h>

#include <itkCenteredTransformInitializer.h>
#include <itkBSplineDeformableTransformInitializer.h>

template <class TFixedImage, class TMovingImage, class TScalarType>
class TransformHelper : public itk::Object
{
  public:
    typedef TransformHelper                             Self;
    typedef itk::Object                                 Superclass;
    typedef itk::SmartPointer< Self >                   Pointer;
    typedef itk::SmartPointer< const Self >             ConstPointer;

    itkTypeMacro( TransformHelper, Object );

    itkNewMacro( Self );
    
    /** Some convenient typedefs. */
    typedef TFixedImage                            FixedImageType;
    typedef typename FixedImageType::Pointer       FixedImagePointer;
    typedef typename FixedImageType::ConstPointer  FixedImageConstPointer;
    typedef TMovingImage                           MovingImageType;
    typedef typename MovingImageType::Pointer      MovingImagePointer;
    typedef typename MovingImageType::ConstPointer MovingImageConstPointer;

    itkStaticConstMacro(FixedImageDimension, unsigned int, TFixedImage::ImageDimension);
    itkStaticConstMacro(MovingImageDimension, unsigned int, TMovingImage::ImageDimension);
    
    /** Get the size of the input space */
    unsigned int GetFixedImageDimension(void) const {return FixedImageDimension;}

    /** Get the size of the output space */
    unsigned int GetMovingImageDimension(void) const {return MovingImageDimension;}

    /** Type of the scalar representing coordinate and vector elements. */
    typedef  TScalarType     ScalarType;
    
    
    enum TRANSFORM_TYPE {
      EULER_3D,
      EULER_3D_CENTERED,
      VERSOR_RIGID_3D,
      VERSOR_NON_RIGID_3D,
      AFFINE_7P,
      AFFINE_9P,
      AFFINE_15P,
      BSPLINE
    };

    enum INITIALIZATION_TYPE {
      INIT_NONE,
      INIT_GEOMETRY,
      INIT_MOMENTS,
      INIT_AFFINE,
      INIT_PRINCIPAL_AXIS,
      INIT_TRANSFORM
    };
    
    typedef itk::Transform< ScalarType, FixedImageDimension, FixedImageDimension>            GenericType;
    typedef itk::MatrixOffsetTransformBase
                          < ScalarType, FixedImageDimension, FixedImageDimension>            RigidType;
    typedef itk::CenteredEuler3DTransform
                          < ScalarType >            CenteredEulerType;
    typedef itk::Euler3DTransform
                          < ScalarType >            EulerType;
    typedef itk::VersorRigid3DTransform
                          < ScalarType >            VersorRigidType;
    typedef itk::VersorTransform
                          < ScalarType >            VersorType;
    typedef itk::ScaleSkewVersor3DTransform
                          < ScalarType >            Affine15PType;
    typedef itk::ScaleVersor3DTransform
                          < ScalarType >            Affine9PType;
    typedef itk::Similarity3DTransform
                          < ScalarType >            Affine7PType;
    typedef itk::AffineTransform
                          < ScalarType, FixedImageDimension >                                AffineType;
    typedef itk::CenteredAffineTransform
                              < ScalarType, FixedImageDimension >                            CenteredAffineType;
    //typedef itk::BSplineDeformableTransform

    typedef std::vector< boost::any >                                                        ParametersType;
    typedef itk::Array< double >                                                             ScalesType;
    
    template <class TRegistration>
    void Connect( TRegistration* r);
    
    void SetTransformType( TRANSFORM_TYPE type );
    
    void SetInitializationType( INITIALIZATION_TYPE type )
    { m_InitType = type; }
    
    void SetInitializationImages( FixedImageConstPointer fixed, MovingImageConstPointer moving, bool swapLevel = false )
    { m_FixedImage = fixed; m_MovingImage = moving; m_SwapLevel=swapLevel; }

    void SetInitialization( INITIALIZATION_TYPE type, FixedImageConstPointer fixed, MovingImageConstPointer moving, bool swapLevel = false )
    { m_InitType = type; m_FixedImage = fixed; m_MovingImage = moving; m_SwapLevel=swapLevel; }
       
    void SetInitialization( FixedImageConstPointer fixed, MovingImageConstPointer moving, typename GenericType::Pointer transform )
    { m_FixedImage = fixed; m_MovingImage = moving; m_InitTransform = transform; }
    
    void SetInitialTransform( GenericType* transform )
    { m_InitTransform = transform; m_InitType = INIT_TRANSFORM; }

    
    typename GenericType::Pointer GetTransform()
    { return m_Transform; }
    
    typename GenericType::Pointer GetInverseTransform();
    
  
    template <class TInputImage, class TReferenceImage>
    TReferenceImage* ResampleImage( typename TInputImage::ConstPointer input, typename TReferenceImage::ConstPointer ref, bool useInverse = false, std::string fileName = "" );
    
    void SaveTransformToFile( std::string fileName );

    itkGetConstMacro(NumberOfParameters, unsigned int);
    itkGetConstMacro(Type, TRANSFORM_TYPE );

  protected:
    TransformHelper();
    ~TransformHelper() {};

    virtual void PrintSelf(std::ostream& os, itk::Indent indent) const;
  private:
    TRANSFORM_TYPE                       m_Type;
    INITIALIZATION_TYPE                  m_InitType;
    ParametersType                       m_Parameters;
    bool                                 m_SwapLevel;
    unsigned int                         m_NumberOfParameters;
    
    FixedImageConstPointer               m_FixedImage;
    MovingImageConstPointer              m_MovingImage;
    typename GenericType::Pointer        m_InitTransform;
    
    typename GenericType::Pointer        m_Transform;
    typename RigidType::Pointer          m_RigidTransform;
    typename CenteredEulerType::Pointer  m_CenteredEulerTransform;
    typename EulerType::Pointer          m_EulerTransform;
    typename VersorRigidType::Pointer    m_VersorRigidTransform;
    typename VersorType::Pointer         m_VersorTransform;
    typename Affine15PType::Pointer      m_Affine15PTransform;
    typename Affine9PType::Pointer       m_Affine9PTransform;
    typename Affine7PType::Pointer       m_Affine7PTransform;
    typename AffineType::Pointer         m_AffineTransform;
    typename CenteredAffineType::Pointer m_CenteredAffineTransform;

    
    template <class TTransform>
    typename TTransform::Pointer InitializeTransform();

    template <class TTransform>
    typename TTransform::Pointer Instantiate();
};

#include "TransformHelper.txx"

#endif // TRANSFORMHELPER_H
