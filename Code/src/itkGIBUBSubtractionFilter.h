/* --------------------------------------------------------------------------------------
 * File:    itkGIBUBSubtractionFilter.h
 * Date:    22/08/2011
 * Author:  Berta Martí bmarti@ub.es
 * Version: 0.1
 * License: BSD
 * --------------------------------------------------------------------------------------

 Copyright (c) 2011, Berta Martí - bmarti@ub.es
 with Grupo de Imágenes Biomédicas, Universidad de Barcelona (GIB-UB)
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

 THIS SOFTWARE IS PROVIDED BY Berta Martí ''AS IS'' AND ANY
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


#ifndef ITKGIBUBSUBTRACTIONFILTER_H_
#define ITKGIBUBSUBTRACTIONFILTER_H_

#include <itkBinaryFunctorImageFilter.h>

namespace itk
{

/** \class SubtractImageFilter
 * \brief Implements an operator for pixel-wise subtraction of two images.
 *
 * Output = Input1 - Input2.
 *
 * This class is parametrized over the types of the two
 * input images and the type of the output image.
 * Numeric conversions (castings) are done by the C++ defaults.
 *
 * \ingroup IntensityImageFilters Multithreaded
 */
namespace Function {

template< class TInput1, class TInput2=TInput1, class TOutput=TInput1>
class Sub2
{
public:
  Sub2() {}
  ~Sub2() {}
  bool operator!=( const Sub2 & ) const
    {
    return false;
    }
  bool operator==( const Sub2 & other ) const
    {
    return !(*this != other);
    }
  inline TOutput operator()( const TInput1 & A, const TInput2 & B) const
    {
		if( A != (TInput1)0 && B != (TInput2)0 )
		{

			if ((A - B)*100/A <=100 && (A - B)*100/A >= -100)
			{
				return (TOutput)((A - B)*100/A);
			}
			else {
				if ((A - B)*100/A >100)
				{
					return (TOutput)(100);
				}
				else{
					return (TOutput)(-100);
				}
			}
		}
		else
		{
			return (TOutput)(0);
		}
	}
};
}

template <class TInputImage1, class TInputImage2=TInputImage1, class TOutputImage=TInputImage1>
class ITK_EXPORT GIBUBSubtractionFilter :
    public
BinaryFunctorImageFilter<TInputImage1,TInputImage2,TOutputImage,
                         Function::Sub2<
  typename TInputImage1::PixelType,
  typename TInputImage2::PixelType,
  typename TOutputImage::PixelType>   >
{
public:
  /** Standard class typedefs. */
  typedef GIBUBSubtractionFilter                         Self;
  typedef BinaryFunctorImageFilter<
    TInputImage1,TInputImage2,TOutputImage,
    Function::Sub2< typename TInputImage1::PixelType,
                    typename TInputImage2::PixelType,
                    typename TOutputImage::PixelType> >  Superclass;
  typedef SmartPointer<Self>                             Pointer;
  typedef SmartPointer<const Self>                       ConstPointer;


  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(GIBUBSubtractionFilter,
               BinaryFunctorImageFilter);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(Input1Input2OutputAdditiveOperatorsCheck,
    (Concept::AdditiveOperators<typename TInputImage1::PixelType,
                                typename TInputImage2::PixelType,
                                typename TOutputImage::PixelType>));
  /** End concept checking */
#endif

protected:
  GIBUBSubtractionFilter() {}
  virtual ~GIBUBSubtractionFilter() {}

private:
  GIBUBSubtractionFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented


};

} // end namespace itk


#endif /* ITKGIBUBSUBTRACTIONFILTER_H_ */
