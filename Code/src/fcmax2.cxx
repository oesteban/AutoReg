#include "fcmax2.h"


double ComputeNormalizationFactor(InternalImageType::Pointer im_i, InternalImageType::Pointer im_ii, InternalImageType::Pointer mask) {
	double factor = 0.0;

	//Inital variables
	int num;
	double nint = 100.0;

	//Dimensions of target images
	int MMM = im_i->GetLargestPossibleRegion().GetNumberOfPixels();
	
	//Histogram length 
	int NINT=100;
	int P=10*NINT-1;					

	//Initialize the vectors of histogram (l and h)
	VectorDouble l;
	l.resize(P);
	VectorDouble h;
	h.resize(P);

	//Initialize the vector that will contains the ratio between Target1/Target2
	VectorDouble x;
	x.resize(MMM);

/*
	//Profile ratios to get the inital factor (I/II)
	typedef itk::ImageRegionIterator<InputImageType> Iterator;
	Iterator it(mask, mask->GetLargestPossibleRegion());

	for (it = it.Begin(); !it.IsAtEnd(); ++it) {
		InputImageType::IndexType idx = it.GetIndex();
		unsigned int i = mask->ComputeOffset(idx);
		if ( (it.Get() * im_ii->GetPixel(idx)) > 0){
			x[i] = ( im_i->GetPixel(idx) / im_ii->GetPixel(idx) );

			if (x[i]>0.) {
				num = (int) floor( x[i]*nint+.5);
				if (num<P) h[num]+=1;
			}
		}
	}
*/

	for (int i=0;i<MMM;i++) {		
		if (mask->GetBufferPointer()[i]*im_ii->GetBufferPointer()[i]>0) {
			x[i]=(im_i->GetBufferPointer()[i]/im_ii->GetBufferPointer()[i]);

			if (x[i]>0.) {
				num =(int)floor(x[i]*nint+.5);
				if (num<P)	h[num]=h[num]+1.;
			}
		} 
	} 

	//Smooth histogram
	for (int i=2;i<P-2;i++) {
		l[i]=(h[i-2]+h[i-1]+h[i]+h[i+1]+h[i+2])/5.;
	}
   
	// Get the maximum of l vector and its position
	double maxValue = *std::max_element(l.begin(), l.end());
	VectorDouble::iterator iterator = max_element(l.begin(), l.end());
	int positionMax = std::distance(l.begin(), iterator);

    //Average
	double N=0.0;
	double sx=0.0;
	for (int i=positionMax-NINT/20; i<=positionMax+NINT/20; i++) {		
		N += h[i];
		sx += i*h[i]/nint;
	}
   
	double average = sx/N;

#ifndef NDEBUG
	std::cout << std::endl << "Maximum value: " << maxValue << " Average: " << average << " Total number of values: " << N << std::endl;
#endif

	//Fiting a parabola  y=ax2+bx+c 
	//Init variables
	double sx1=0.0;
	double sx2=0.0;
	double sx3=0.0;
	double sx4=0.0;
	double sy=0.0;
	double sxy=0.0;
	double sx2y=0.0;	
	
	double ix;
	for(int i=positionMax-NINT/10;i<positionMax+NINT/10;i++){
		ix=(float)(i);
		sx1+=ix;
		sx2+=ix*ix;
		sx3+=ix*ix*ix;
		sx4+=ix*ix*ix*ix;
		sy+=h[i];
		sxy+=ix*h[i];
	    sx2y+=ix*ix*h[i];	
	}

	double n_p=nint/5.0;
	
	double a, b, c, denominador;
	determinant(sx4,sx3,sx2,sx3,sx2,sx1,sx2,sx1,n_p,&denominador);
	determinant(sx2y,sx3,sx2,sxy,sx2,sx1,sy,sx1,n_p,&a);
	determinant(sx4,sx2y,sx2,sx3,sxy,sx1,sx2,sy,n_p,&b);
	determinant(sx4,sx3,sx2y,sx3,sx2,sxy,sx2,sx1,sy,&c);
	
	a/=denominador;
	b/=denominador;
	c/=denominador;

	double new_f=-(b)/(2.*a)/nint;
	double new_h=-b*b/(4.*a)+c;

	//New factor that comes from parabola fitted
	factor=new_f;

#ifndef NDEBUG
	std::cout << "New Factor: " << new_f << std::endl; 
#endif	

 	return factor;
}
   
//**************************************calcular determinant*****************************************/
void determinant(double a1,double a2,double a3,double a4,double a5,double a6,double a7,double a8,double a9,double *det)
{
	*det=a1*a5*a9+a4*a8*a3+a2*a6*a7-a3*a5*a7-a1*a6*a8-a2*a4*a9;
}
