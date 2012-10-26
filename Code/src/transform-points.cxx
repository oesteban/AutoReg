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

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkTransformFileReader.h>
#include <itkVersorRigid3DTransform.h>


int main( int argc, char *argv[] )
{
    if ( argc < 2 )
    {
        std::cerr << "Missing Parameters " << std::endl;
        std::cerr << "Usage: " << argv[0];
        std::cerr << " itkTransform";
        std::cerr <<  std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::VersorRigid3DTransform< double > Transform;
    // Read FocusDet transform from file
    itk::TransformFileReader::Pointer readtfm = itk::TransformFileReader::New();
    std::string tfFileName (argv[1]);
    readtfm->SetFileName( tfFileName );
    readtfm->Update();


    Transform::Pointer tf_test;

    typedef itk::TransformFileReader::TransformListType * TFList;
    TFList transforms = readtfm->GetTransformList();
    itk::TransformFileReader::TransformListType::const_iterator it;
    it = transforms->begin();

    if ( !strcmp( (*it)->GetNameOfClass(), "VersorRigid3DTransform"))
    {
        tf_test = static_cast< Transform * > ( (*it).GetPointer() ) ;
    }
    
    tf_test->GetInverse( tf_test );

    typedef itk::Image< float, 3u >          Image;
    typedef itk::ImageFileReader<Image>      Reader;
    Reader::Pointer r = Reader::New();
    r->SetFileName(argv[2]);
    r->Update();
    Image::Pointer im = r->GetOutput();
    Image::SizeType size = im->GetLargestPossibleRegion().GetSize();
    
    std::vector< Image::PointType > innerPoints;
    std::vector< Image::PointType > outerPoints;
//     innerPoints.resize(8);
//     outerPoints.resize(8);
    
    for (unsigned int i = 0; i<2; i++)
      for (unsigned int j = 0; j<2; j++)
	for (unsigned int k = 0; k<2; k++)
    {
      Image::PointType p;
      Image::IndexType idx;
      idx[0] = 0 + (size[0]-1)*k;
      idx[1] = 0 + (size[1]-1)*j;
      idx[2] = 0 + (size[2]-1)*i;
      im->TransformIndexToPhysicalPoint(idx, p);
      innerPoints.push_back( p );
    }


    for ( unsigned int i = 0; i < innerPoints.size(); i++)
    {
        outerPoints.push_back(tf_test->TransformPoint( innerPoints[i] ));
        std::cout << "   " << (i+1) << "   ";
	std::cout<< std::showpoint <<  std::setw(15) << std::setprecision(7) << (double) innerPoints[i][0];
        std::cout<< std::showpoint <<  std::setw(15) << std::setprecision(7) <<  (double) innerPoints[i][1];
        std::cout<< std::showpoint <<  std::setw(15) << std::setprecision(7) <<  (double) innerPoints[i][2];
        std::cout<< std::showpoint <<  std::setw(15) << std::setprecision(7) <<  (double) outerPoints[i][0];
        std::cout<< std::showpoint <<  std::setw(15) << std::setprecision(7) <<  (double) outerPoints[i][1];
        std::cout<< std::showpoint <<  std::setw(15) << std::setprecision(7) <<  (double) outerPoints[i][2] << std::endl;

    }
}