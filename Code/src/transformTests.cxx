/* ==============================================================================================

* ================================================================================================
*/
//

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>
#include <itkTransformFileReader.h>
#include <itkVersorRigid3DTransform.h>
#include <itkEuler3DTransform.h>
#include <itkGIBUB3DTransform.h>
#include <itkOrientImageFilter.h>
#include <itkTranslationTransform.h>

#include <itkContinuousIndex.h>
#include <itkImageRegionIterator.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#define DEG M_PI/180.0


int main( int argc, char *argv[] )
{
    if ( argc < 2 )
    {
        std::cerr << "Missing Parameters " << std::endl;
        std::cerr << "Usage: " << argv[0];
        std::cerr << "case_id goldTransformsList testTransformFile fixedImage";
        std::cerr <<  std::endl;
        return EXIT_FAILURE;
    }

    //typedef  unsigned char  PixelRange;
    //const unsigned int Dimension = 2;
    typedef  float  PixelRange;
    const int Dimension = 3;

    typedef itk::Image< PixelRange, Dimension >          Image;
    typedef itk::ImageRegionIterator< Image >            Iterator;
    typedef itk::VersorRigid3DTransform< double >        VersorTransform;
    typedef itk::Euler3DTransform< double >              EulerTransform;
    typedef itk::ImageFileReader< Image  >               Reader;
    typedef itk::ImageFileWriter< Image >                Writer;
    typedef itk::ResampleImageFilter< Image, Image>      Resampler;
    typedef itk::TranslationTransform<double,Dimension>  Translation;

    
    bool lesionMask = true;
//     bool outputImages = (argc>=4 && atoi(argv[3])==1 );
//     bool landmarks = (argc>=3 && atoi(argv[2])== 1);
    unsigned int case_id = atoi(argv[1]);
    std::string dir( argv[3] );
    std::string tfFileName = dir + "/tforms.txt";

    
    // ----------------------------------------------------------------
    // LECTURA DE LAS IMÁGENES
    // ----------------------------------------------------------------
    Reader::Pointer r_fixed  = Reader::New();
    Reader::Pointer r_moving = Reader::New();
    Reader::Pointer r_mask   = Reader::New();
    Writer::Pointer w  = Writer::New();
    
    r_fixed->SetFileName(  argv[4] );
    r_moving->SetFileName(  argv[5] );
    r_mask->SetFileName( argv[6] );
    try {
        r_fixed->Update();
        r_moving->Update();
    }
    catch ( itk::ExceptionObject & err )
    {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        return EXIT_FAILURE;
    }

    try {
        r_mask->Update();
    } catch (...)
    {
        lesionMask = false;
    }

    Image::Pointer i_bas = r_moving->GetOutput();
    i_bas->Update();
    Image::Pointer i_ref = r_fixed->GetOutput();
    i_ref->Update();
    
    Image::SizeType      s_ref  = i_ref->GetLargestPossibleRegion().GetSize();
    Image::PointType center_ref;
    double center_index_ref[3];
    center_index_ref[0] = (s_ref[0]-1) * 0.5;
    center_index_ref[1] = (s_ref[1]-1) * 0.5;
    center_index_ref[2] = (s_ref[2]-1) * 0.5;
//    center_index_ref[0]=127; center_index_ref[1]=127; center_index_ref[2]=57;
    i_ref->TransformContinuousIndexToPhysicalPoint<double>( center_index_ref , center_ref );

    Image::SizeType      s_bas  = i_bas->GetLargestPossibleRegion().GetSize();
    Image::PointType center_bas;
    double center_index_bas[3];
    center_index_bas[0] = (s_bas[0]-1) * 0.5;
    center_index_bas[1] = (s_bas[1]-1) * 0.5;
    center_index_bas[2] = (s_bas[2]-1) * 0.5;
//    center_index_bas[0]=center_index_bas[1]=64; center_index_bas[2]=27;
    i_bas->TransformContinuousIndexToPhysicalPoint<double>( center_index_bas , center_bas );

    Image::PointType corner_bas;
    double corner_index_bas[3];
    corner_index_bas[0] = s_bas[0]-1;
    corner_index_bas[1] = s_bas[1]-1;
    corner_index_bas[2] = s_bas[2]-1;
    i_bas->TransformContinuousIndexToPhysicalPoint<double>( corner_index_bas , corner_bas );
    Image::IndexType origin_index_bas; origin_index_bas.Fill( 0 );
    Image::PointType origin_bas;
    i_bas->TransformIndexToPhysicalPoint( origin_index_bas , origin_bas );
   
    itk::Vector<double, 3u> centers_diff = center_bas - center_ref;
    itk::Vector<double, 3u> diff_ref = center_ref - i_ref->GetOrigin();
    itk::Vector<double, 3u> diff_bas = center_bas - i_bas->GetOrigin();
//    i_ref->SetOrigin( i_ref->GetOrigin() + diff );
//    i_bas->SetOrigin( i_bas->GetOrigin() - diff );

    Image::PointType new_center_ref = center_ref;
//    i_ref->TransformContinuousIndexToPhysicalPoint<double>( center_index_ref , new_center_ref );
    Image::PointType new_center_bas;
//    i_bas->TransformContinuousIndexToPhysicalPoint<double>( center_index_bas , new_center_bas );
    new_center_bas = center_bas - diff_ref; 
    
    // TRANSFORMS CONFIGURATION ---------------------------------------------------------
    typedef itk::TransformFileReader::TransformListType * TFList;
    
    // Create and generate a Gold Registration transform.
    EulerTransform::Pointer tf_gold;
    itk::TransformFileReader::Pointer tf_gold_reader = itk::TransformFileReader::New();
    tf_gold_reader->SetFileName( argv[2] );
    tf_gold_reader->Update();
    TFList tfs_gold_list = tf_gold_reader->GetTransformList();
    itk::TransformFileReader::TransformListType::const_iterator it_gold = tfs_gold_list->begin();

    std::advance( it_gold, case_id );
    
    if ( !strcmp( (*it_gold)->GetNameOfClass(), "Euler3DTransform"))
    {
        tf_gold = static_cast< EulerTransform * > ( (*it_gold).GetPointer() ) ;
    }
    else
    {
      std::cout << "No transform found" << std::endl;
    }
    EulerTransform::Pointer tf_gold_t = EulerTransform::New();
    tf_gold_t->SetIdentity();
    tf_gold_t->SetOffset( diff_ref );

    EulerTransform::Pointer tf_gold_i = EulerTransform::New();
    tf_gold_i->SetCenter( tf_gold->GetCenter() );
    tf_gold->GetInverse( tf_gold_i );

    std::cout << "Oref =" << center_ref << std::endl;
/*
    Image::PointType new_center = tf_gold_i->TransformPoint( center_ref );

    std::cout << "Tgold(Oref)=" << new_center << std::endl;
   
    std::cout << "Tgold-1(Oref)=" << tf_gold->TransformPoint( new_center ) << "?=" << center_ref << std::endl;
*/
    EulerTransform::Pointer tf_gold_t_i = EulerTransform::New();
    tf_gold_t_i->SetCenter( tf_gold_t->GetCenter() );
    tf_gold_t->GetInverse( tf_gold_t_i );
/*
    Image::PointType o2 = tf_gold_t_i->TransformPoint( new_center );
    std::cout << "Tgoldt(Tgold(Oref))=" << o2 << std::endl;

    Image::PointType o1 = tf_gold_t->TransformPoint( o2 );
    std::cout << "Tgoldt-1(Tgoldt(Tgold(Oref)))=" << o1 << "?=" << new_center << std::endl;

    Image::PointType oref = tf_gold->TransformPoint( o1 );
    std::cout << "Tgold-1( Tgoldt-1 ( Tgoldt ( Tgold (Oref) ) ) )=" << oref << "?=" << center_ref << std::endl;
*/
    tf_gold->Compose( tf_gold_t, true );

/*
    std::cout << "Ref origin2?=" << tf_gold->TransformPoint( o2 ) << std::endl;


    Resampler::Pointer res_gold_moving = Resampler::New();
    res_gold_moving->SetInput( i_ref );
    res_gold_moving->SetTransform( tf_gold );
    res_gold_moving->SetUseReferenceImage( true );
    res_gold_moving->SetReferenceImage( i_bas );
    w->SetInput( res_gold_moving->GetOutput() );
    w->SetFileName( dir +  "/result_gold_moving.hdr" );
    w->Update();
*/

    EulerTransform::Pointer tf_gold_inverse = EulerTransform::New();
    tf_gold_inverse->SetCenter( tf_gold->GetCenter() );
    tf_gold->GetInverse(tf_gold_inverse);

    EulerTransform::Pointer tf_gold_t_inv = EulerTransform::New();
    tf_gold_t_inv->SetIdentity();
    tf_gold_t_inv->SetOffset( diff_bas );

  /*  
    std::cout << "Gold_tf_center=" << tf_gold->GetCenter() << std::endl;
    std::cout << "origin_bas=" << origin_bas << std::endl;
    std::cout << "corner_bas=" << corner_bas << std::endl;
    std::cout << "center_ref=" << center_ref << std::endl;
    std::cout << "center_bas=" << center_bas << std::endl;
    std::cout << "diff_bas = " << diff_bas << std::endl;
    std::cout << "diff_ref = " << diff_ref << std::endl;
    std::cout << "centers_dif = " << centers_diff << std::endl;
*/

    tf_gold_inverse->Compose( tf_gold_t_inv );
  

    // Read FocusDet transform from file
    VersorTransform::Pointer tf_test_inverse;
    itk::TransformFileReader::Pointer tf_test_reader = itk::TransformFileReader::New();
    tf_test_reader->SetFileName( tfFileName );
    tf_test_reader->Update();
    TFList tfs_test_list = tf_test_reader->GetTransformList();
    itk::TransformFileReader::TransformListType::const_iterator it_test;
    it_test = tfs_test_list->begin();

    if ( !strcmp( (*it_test)->GetNameOfClass(), "VersorRigid3DTransform"))
    {
        tf_test_inverse = static_cast< VersorTransform * > ( (*it_test).GetPointer() ) ;
    }
   
    VersorTransform::Pointer tf_test = VersorTransform::New();
    tf_test->SetCenter( tf_test_inverse->GetCenter() );
    tf_test_inverse->GetInverse( tf_test );
    
    // END TRANSFORMS CONFIGURATION -----------------------------------------------------


    Resampler::Pointer res_gold = Resampler::New();
    res_gold->SetInput( i_bas );
    res_gold->SetTransform( tf_gold_inverse );
    res_gold->SetUseReferenceImage( true );
    res_gold->SetReferenceImage( i_ref );
    w->SetInput( res_gold->GetOutput() );
    w->SetFileName( dir +  "/result_gold_ref.hdr" );
    w->Update();

    Resampler::Pointer res_test = Resampler::New();
    res_test->SetInput( i_bas );
    res_test->SetTransform( tf_test_inverse );
    res_test->SetUseReferenceImage( true );
    res_test->SetReferenceImage( i_ref );
    w->SetInput( res_test->GetOutput() );
    w->SetFileName( dir +  "/result_test_ref.hdr" );
    w->Update();



    double totalMSEError = 0.0;
    int totalMSEVoxels = 0;

    for ( int i = -1; i<2; i++ )
    {
        for (int j= -1; j<2; j++)
        {
            for (int k = -1; k<2; k++)
            {
	      if ( i== 0 || j==0 || k==0) continue;
	      
	      Image::PointType p_ref = new_center_ref;
              p_ref[0] += 60 * i;
              p_ref[1] += 60 * j;
              p_ref[2] += 60 * k;
//               std::cout << "Pref[" << totalMSEVoxels<< "]="<< p_ref << std::endl;
	      
	      //Image::PointType p;
              Image::PointType p_moving = tf_gold_inverse->TransformPoint( p_ref );
// 	      p_moving += diff;
// 	      std::cout << "Pmoving[" << totalMSEVoxels<< "]="<< p_moving << std::endl;
	      
	      Image::PointType p_test = tf_test->TransformPoint( p_moving );
// 	      std::cout << "Ptest[" << totalMSEVoxels<< "]="<< p_test << std::endl;

	      double error = p_ref.EuclideanDistanceTo<Image::PointType::CoordRepType>( p_test );
// 	      std::cout << "Error[" << totalMSEVoxels<< "]="<< error << std::endl;
	      totalMSEError+=error;
	      totalMSEVoxels++;
            }
        }
    }
    std::cout << "Results for RM" << std::setw(2) << std::setfill('0') << case_id << ":" << std::endl;
    std::cout << "TRE= " << ( totalMSEError / totalMSEVoxels ) << " mm." << std::endl;


    std::vector < Image::PointType > lesion_points;
    if ( lesionMask )
    {
        Image::Pointer lesion_mask = r_mask->GetOutput();
        Iterator it (lesion_mask, lesion_mask->GetLargestPossibleRegion() );

        for ( it = it.Begin(); !it.IsAtEnd(); ++it )
        {
	  if( it.Get()>0.0)
	  {
	    Image::PointType point;
	    lesion_mask->TransformIndexToPhysicalPoint( it.GetIndex(), point );
	    lesion_points.push_back( point );
	  }
        }
    }
    
    if ( !lesion_points.empty() )
    {
        double totalError = 0.0;
        unsigned int totalVoxels = 0;
        std::vector< Image::PointType >::const_iterator it = lesion_points.begin();

	while ( it != lesion_points.end() )
        {
	  Image::PointType p_ref = *it;
              //Image::PointType p;
              Image::PointType p_moving = tf_gold_inverse->TransformPoint( p_ref );
	      //p_moving += diff;
	      
	    Image::PointType p_test = tf_test->TransformPoint( p_moving );

            double error = p_ref.EuclideanDistanceTo<Image::PointType::CoordRepType>( p_test );

            totalError+=error;
            totalVoxels++;
            ++it;
        }
        std::cout << "WI= " << ( totalError / totalVoxels ) << " mm." << std::endl;
    }
// 
//     if (outputImages)
//     {
//       // Gold standard registered generation
//       // Base image grid
//       Resampler::Pointer res_gold = Resampler::New();
//       res_gold->SetInput( i_bas );
//       res_gold->SetTransform( tf_gold );
//       res_gold->SetUseReferenceImage( true );
//       res_gold->SetReferenceImage( i_bas );
//       w->SetInput( res_gold->GetOutput() );
//       w->SetFileName( dir.str() +  "result_" + case_id_str + "_gold_bas.hdr" );
//       w->Update();
// 
//       res_gold->SetUseReferenceImage(false);
//       res_gold->SetOutputOrigin( new_ref_origin );
//       res_gold->SetOutputDirection( i_ref->GetDirection() );
//       res_gold->SetOutputSpacing( i_ref->GetSpacing() );
//       res_gold->SetSize( i_ref->GetLargestPossibleRegion().GetSize() );
//       res_gold->Update();
//       w->SetFileName( dir.str() +  "result_" + case_id_str +"_gold_ref.hdr" );
//       w->Update();
// 
//       // Test registered generation
//       // Base image grid
//       Resampler::Pointer resample_test_ref = Resampler::New();
//       resample_test_ref->SetInput( i_bas );
//       resample_test_ref->SetTransform( tf_test );
//       resample_test_ref->SetUseReferenceImage(false);
//       resample_test_ref->SetOutputOrigin( new_bas_origin );
//       resample_test_ref->SetOutputDirection( i_bas->GetDirection() );
//       resample_test_ref->SetOutputSpacing( i_bas->GetSpacing() );
//       resample_test_ref->SetSize( i_bas->GetLargestPossibleRegion().GetSize() );
//       resample_test_ref->Update();
//       Image::Pointer test_bas_im = resample_test_ref->GetOutput();
// 
//       Writer::Pointer w_test = Writer::New();
//       w_test->SetInput( resample_test_ref->GetOutput() );
//       w_test->SetFileName( dir.str() +  "result_" + case_id_str + "_test_bas.hdr" );
//       w_test->Update();
// 
//       resample_test_ref->SetUseReferenceImage(true);
//       resample_test_ref->SetReferenceImage( i_ref );
//       resample_test_ref->Update();
//       w_test->SetFileName( dir.str() +  "result_" + case_id_str + "_test_ref.hdr" );
//       w_test->Update();
//     }
    
    return EXIT_SUCCESS;
}

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

// -------------------------------------------------------------------------------
// transformITKImage.cxx
// v1.0 - 2011-01-27
//
// This code performs a geometric transform on the input image to an output image
// using the transform number case_id on the file transformsItkFile.
// This file must have the ITKTransform list format.
// -------------------------------------------------------------------------------
/*
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>
#include <itkTransformFileReader.h>
#include <itkVersorRigid3DTransform.h>
#include <itkGIBUB3DTransform.h>
#include <itkOrientedImage.h>
#include <itkOrientImageFilter.h>

#include <itkContinuousIndex.h>
#include <itkImageRegionIterator.h>
#include <itkTransformFileWriter.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#define DEG M_PI/180.0


namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;


int main( int argc, char *argv[] )
{
    typedef short  PixelRange;
    const int Dimension = 3;

    typedef itk::Image< PixelRange, Dimension >          Image;
    typedef itk::ImageRegionIterator< Image >            Iterator;
    typedef itk::VersorRigid3DTransform< double >        Transform;
    typedef itk::ImageFileReader< Image  >               Reader;
    typedef itk::ImageFileWriter< Image >                Writer;
    typedef itk::ResampleImageFilter< Image, Image>      Resampler;
    typedef itk::OrientImageFilter< Image, Image >       Orienter;


    std::string interictalFileName;
    std::string ictalFileName;
    std::string referenceFileName;
    std::string transformsFileName;
    std::string outputPrefix;
    bool output_images = false;
    bool output_ref = true;
    bool output_std = false;
    bool output_bas = false;
    bool output_act = false;

	bpo::options_description desc("Usage");
	desc.add_options()
			("help", "show help message")
			("interictal-image", bpo::value< std::string >(&interictalFileName)->required(), "interictal image file")
			("ictal-image", bpo::value< std::string >(&ictalFileName)->required(), "ictal image file")
			("reference-image,R", bpo::value< std::string >(&referenceFileName)->required(), "reference image file")
			("transforms-file,F"
					"", bpo::value< std::string >(&transformsFileName)->required(), "itk::VersorRigid3DTranform transforms file in ITK format. First transform=coregistration. Second transform=registration.")
			("output-prefix,o", bpo::value< std::string >(&outputPrefix), "prefix for output files")
			("output-all,a", bpo::bool_switch(&output_images)->default_value(false), "output all images")
			("output-ref,r", bpo::bool_switch(&output_ref)->default_value(true), "output images in reference space")
			("output-bas,b", bpo::bool_switch(&output_bas)->default_value(false), "output images in basal space")
			("output-act,i", bpo::bool_switch(&output_act)->default_value(false), "output images in activated space")
			("output-std,s", bpo::bool_switch(&output_std)->default_value(false), "output images in standard space [128 128 54] [3.32 3.32 3.32]");

	bpo::positional_options_description pdesc;
	pdesc.add("interictal-image", 1);
	pdesc.add("ictal-image", 2);

	bpo::variables_map vmap;
	bpo::store(bpo::command_line_parser(argc, argv).options(desc).positional(pdesc).run(), vmap);

	if (vmap.count("help")) {
	    std::cout << desc << std::endl;
	    return 1;
	}

	try {
		bpo::notify(vmap);
	} catch ( boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::program_options::required_option> > &err){
		std::cout << "Error: " << err.what() << std::endl;
		std::cout << desc << std::endl;
		return EXIT_FAILURE;
	}


    if (!bfs::exists( transformsFileName ) )
    {
      std::cout << "Error: Transforms file not found (supplied path=" << transformsFileName << ")." << std::endl;
      return EXIT_FAILURE;
    }


    // ----------------------------------------------------------------
    // LECTURA DE LAS IMÁGENES
    // ----------------------------------------------------------------
    Reader::Pointer r  = Reader::New();
    Reader::Pointer r2  = Reader::New();
    Reader::Pointer r3  = Reader::New();
    Writer::Pointer w  = Writer::New();


    r->SetFileName(  interictalFileName );
    r2->SetFileName( ictalFileName );
    r3->SetFileName( referenceFileName );
    try {
        r->Update();
        r2->Update();
        r3->Update();
    }
    catch ( itk::ExceptionObject & err )
    {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        return EXIT_FAILURE;
    }

    Image::Pointer i_bas = r->GetOutput();
    Image::Pointer i_act = r2->GetOutput();
    Image::Pointer i_ref = r3->GetOutput();

    Image::PointType o_bas = i_bas->GetOrigin();
    Image::PointType o_act = i_act->GetOrigin();
    Image::PointType o_ref = i_ref->GetOrigin();

    Image::SizeType      s_bas  = i_bas->GetLargestPossibleRegion().GetSize();
    Image::PointType center_bas;
    double center_index_bas[3];
    center_index_bas[0] = (s_bas[0]-1) * 0.5;
    center_index_bas[1] = (s_bas[1]-1) * 0.5;
    center_index_bas[2] = (s_bas[2]-1) * 0.5;
    i_bas->TransformContinuousIndexToPhysicalPoint<double>( center_index_bas, center_bas );

    Image::SizeType      s_ref  = i_ref->GetLargestPossibleRegion().GetSize();
    Image::PointType center_ref;
    double center_index_ref[3];
    center_index_ref[0] = (s_ref[0]-1) * 0.5;
    center_index_ref[1] = (s_ref[1]-1) * 0.5;
    center_index_ref[2] = (s_ref[2]-1) * 0.5;
    i_ref->TransformContinuousIndexToPhysicalPoint<double>( center_index_ref, center_ref );

    Image::SizeType      s_act  = i_act->GetLargestPossibleRegion().GetSize();
    Image::PointType center_act;
    double center_index_act[3];
    center_index_act[0] = (s_act[0]-1) * 0.5;
    center_index_act[1] = (s_act[1]-1) * 0.5;
    center_index_act[2] = (s_act[2]-1) * 0.5;
    i_act->TransformContinuousIndexToPhysicalPoint<double>( center_index_act, center_act );

    Image::PointType o_zero; o_zero.Fill(0.0);

    Transform::InputVectorType ov_ii2i   = o_bas - o_act;
    Transform::InputVectorType ov_rm2i   = o_ref - o_act;
    Transform::InputVectorType ov_rm2ii  = o_ref - o_bas;

    Transform::InputVectorType v_rm2ii   = center_ref - center_bas;
    Transform::InputVectorType v_rm2i    = center_ref - center_act;
    Transform::InputVectorType v_ii2i    = center_bas - center_act;

    // TRANSFORMS CONFIGURATION ---------------------------------------------------------
    Transform::Pointer tf_coreg;
    Transform::Pointer tf_reg;

    itk::TransformFileReader::Pointer tf_gold_reader = itk::TransformFileReader::New();
    tf_gold_reader->SetFileName( transformsFileName );
    tf_gold_reader->Update();

    typedef itk::TransformFileReader::TransformListType * TFList;
    TFList tfs_gold_list = tf_gold_reader->GetTransformList();
    itk::TransformFileReader::TransformListType::const_iterator it_gold = tfs_gold_list->begin();

    Transform::Pointer tf_reg_inv = Transform::New();
    tf_reg_inv->SetIdentity();

    try {
		if ( !strcmp( (*it_gold)->GetNameOfClass(), "VersorRigid3DTransform"))
		{
			tf_coreg = static_cast< Transform * > ( (*it_gold).GetPointer() ) ;
		}
		else
		{
		  std::cout << "No co-registration transform found" << std::endl;
		}

		it_gold++;


		if ( !strcmp( (*it_gold)->GetNameOfClass(), "VersorRigid3DTransform"))
		{
			tf_reg_inv = static_cast< Transform * > ( (*it_gold).GetPointer() ) ;
		}
		else
		{
		  std::cout << "No registration transform found" << std::endl;
		}
    } catch( itk::ExceptionObject & err ) {

    	std::cout << "Error reading transforms" << std::endl << err << std::endl;
    }

    Transform::MatrixType M_reg_inv = tf_reg_inv->GetInverseMatrix();
    Transform::OffsetType O_reg = tf_reg_inv->GetOffset();

    Transform::MatrixType M_coreg = tf_coreg->GetMatrix();
    Transform::OffsetType O_coreg = tf_coreg->GetOffset();

    Transform::MatrixType M_final = M_reg_inv * M_coreg;
    Transform::OffsetType O_final = M_reg_inv * ( O_coreg - O_reg - (o_bas-o_zero) ) + (o_act-o_zero);

    Transform::Pointer tf_coreg_i = Transform::New();
    tf_coreg_i->SetMatrix( M_final );
    tf_coreg_i->SetOffset( O_final );


    if (output_images || output_ref ) {
		Resampler::Pointer res_bas = Resampler::New();
		res_bas->SetUseReferenceImage( true );
		res_bas->SetReferenceImage( i_ref );
		res_bas->SetInput( i_bas );
		res_bas->SetTransform( tf_coreg );
		res_bas->Update();

		w->SetInput( res_bas->GetOutput() );
		w->SetFileName( outputPrefix + "interictal.nii.gz" );
		w->Update();

		Resampler::Pointer res_act = Resampler::New();
		res_act->SetUseReferenceImage( true );
		res_act->SetReferenceImage( i_ref );
		res_act->SetInput( i_act );
		res_act->SetTransform( tf_coreg_i );
		res_act->Update();

		w->SetInput( res_act->GetOutput() );
		w->SetFileName( outputPrefix + "ictal.nii.gz" );
		w->Update();
    }

    if (output_images || output_bas ) {
		Resampler::Pointer res_bas = Resampler::New();
		res_bas->SetUseReferenceImage( false );
		res_bas->SetOutputParametersFromImage( i_bas );
		res_bas->SetSize(s_bas);
		res_bas->SetOutputOrigin( o_bas + v_rm2ii );
		res_bas->SetInput( i_bas );
		res_bas->SetTransform( tf_coreg );
		res_bas->Update();

		w->SetInput( res_bas->GetOutput() );
		w->SetFileName( outputPrefix + "interictal_bas.nii.gz" );
		w->Update();

	    Resampler::Pointer res_act = Resampler::New();
		res_act->SetUseReferenceImage( false );
		res_act->SetOutputParametersFromImage( i_bas );
		res_act->SetSize( s_bas );
		res_act->SetOutputOrigin( o_bas + v_rm2ii );
		res_act->SetInput( i_act );
		res_act->SetTransform( tf_coreg_i );
		res_act->Update();

		w->SetInput( res_act->GetOutput() );
		w->SetFileName( outputPrefix + "ictal_bas.nii.gz" );
		w->Update();
    }

    if (output_images || output_act ) {
		Resampler::Pointer res_bas = Resampler::New();
		res_bas->SetUseReferenceImage( false );
		res_bas->SetOutputParametersFromImage( i_act );
		res_bas->SetSize(s_act);
		res_bas->SetOutputOrigin( o_bas + v_rm2i );
		res_bas->SetInput( i_bas );
		res_bas->SetTransform( tf_coreg );
		res_bas->Update();

		w->SetInput( res_bas->GetOutput() );
		w->SetFileName( outputPrefix + "interictal_act.nii.gz" );
		w->Update();

		Resampler::Pointer res_act = Resampler::New();
		res_act->SetUseReferenceImage( false );
		res_act->SetOutputParametersFromImage( i_act );
		res_act->SetSize( s_act );
		res_act->SetOutputOrigin( o_act + v_rm2i );
		res_act->SetInput( i_act );
		res_act->SetTransform( tf_coreg_i );
		res_act->Update();

		w->SetInput( res_act->GetOutput() );
		w->SetFileName( outputPrefix + "ictal_act.nii.gz" );
		w->Update();
    }

    if (output_images || output_std ) {
    	Image::SizeType s_std; s_std.Fill(128); s_std[2]=54; // Size = [128 128 54]
    	Image::SpacingType sp_std; sp_std.Fill(3.32);        // Spacing = [ 3.32 3.32 3.32
    	Image::PointType o_std; o_std.Fill(0.0);
    	Image::DirectionType d_std = i_ref->GetDirection();

    	Image::PointType center_std;
    	center_std[0] = (s_std[0]-1)*0.5*3.32;
    	center_std[1] = (s_std[1]-1)*0.5*3.32;
    	center_std[2] = (s_std[2]-1)*0.5*3.32;

    	center_std =  d_std * center_std;

    	Transform::InputVectorType v_rm2std   = center_ref - center_std;

		Resampler::Pointer res_bas = Resampler::New();
		res_bas->SetUseReferenceImage( false );
		res_bas->SetSize(s_std);
		res_bas->SetOutputOrigin( o_std + v_rm2std );
		res_bas->SetOutputDirection( d_std );
		res_bas->SetOutputSpacing( sp_std );
		res_bas->SetInput( i_bas );
		res_bas->SetTransform( tf_coreg );
		res_bas->Update();

		w->SetInput( res_bas->GetOutput() );
		w->SetFileName( outputPrefix + "interictal_std.nii.gz" );
		w->Update();

	    Resampler::Pointer res_act = Resampler::New();
	    res_act->SetUseReferenceImage( false );
	    res_act->SetSize(s_std);
	    res_act->SetOutputOrigin( o_std + v_rm2std );
	    res_act->SetOutputDirection( d_std );
	    res_act->SetOutputSpacing( sp_std );
	    res_act->SetInput( i_act );
	    res_act->SetTransform( tf_coreg_i );
	    res_act->Update();

		w->SetInput( res_act->GetOutput() );
		w->SetFileName( outputPrefix + "ictal_std.nii.gz" );
		w->Update();
    }

    return EXIT_SUCCESS;
}
*/
