#include "gaborCLP.h"
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkConvolutionImageFilter.h>
#include <itkGaborImageSource.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>
#include <itkEuler2DTransform.h>
#include <itkEuler3DTransform.h>
#include <fstream>

#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

//To check the image voxel type
void GetImageType( std::string fileName ,
                   itk::ImageIOBase::IOPixelType &pixelType ,
                   itk::ImageIOBase::IOComponentType &componentType,
                   unsigned int &dimension
                   )
{
  typedef itk::Image< unsigned char , 3 > ImageType ;
  itk::ImageFileReader< ImageType >::Pointer imageReader ;
  imageReader = itk::ImageFileReader< ImageType >::New() ;
  imageReader->SetFileName( fileName.c_str() ) ;
  imageReader->UpdateOutputInformation() ;
  pixelType = imageReader->GetImageIO()->GetPixelType() ;
  componentType = imageReader->GetImageIO()->GetComponentType() ;
  dimension = imageReader->GetImageIO()->GetNumberOfDimensions() ;
}


template<class PixelType, unsigned int dimension>
void SaveImage(typename itk::Image< PixelType , dimension >::Pointer image,
               std::string outputVolume,
               std::string suffix,
               std::ofstream &listFilesStream
              )
{
  typedef itk::Image< PixelType , dimension > ImageType ;
  typedef itk::ImageFileWriter< ImageType > WriterType ;
  typename WriterType::Pointer writer = WriterType::New() ;
  std::string outFileName = outputVolume + suffix + ".nrrd";
  std::cout << "Writing output volume: " << outFileName << std::endl;
  writer->SetFileName( outFileName.c_str() ) ;
  writer->SetInput( image ) ;
  writer->SetUseCompression( 1 ) ;
  writer->Update() ;
  if(listFilesStream.is_open())
  {
    listFilesStream << outFileName << std::endl;
  }
}

template<class PixelType, class RealType, unsigned int dimension>
void Process(typename itk::ResampleImageFilter<itk::Image< PixelType , dimension >,
                                               itk::Image< PixelType , dimension >, RealType>::Pointer resampler,
             typename itk::GaborImageSource<itk::Image< PixelType , dimension > >::OutputImageType::Pointer gaborOutput,
             typename itk::Image< PixelType , dimension >::Pointer image,
             std::string outputVolume,
             std::string suffix,
             std::ofstream &listFilesStream,
             std::string outputFilter
             )
{
  typedef itk::Image< PixelType , dimension > ImageType;
  typedef itk::BSplineInterpolateImageFunction< ImageType, RealType, double > InterpolatorType ;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New() ;
  interpolator->SetSplineOrder(3);
  resampler->SetInterpolator( interpolator );
  resampler->SetInput( gaborOutput );
  resampler->SetOutputSpacing( gaborOutput->GetSpacing() );
  resampler->SetOutputOrigin( gaborOutput->GetOrigin() );
  resampler->SetSize( gaborOutput->GetLargestPossibleRegion().GetSize() );
  resampler->Update();
  //Convolve gabor filter with input image
  std::cout<<"Convolution between gabor filter ("<<  suffix <<") and input image"<<std::endl;
  typedef itk::ConvolutionImageFilter<ImageType> ConvolutionFilterType;
  typename ConvolutionFilterType::Pointer convoluter
      = ConvolutionFilterType::New();
  convoluter->SetInput( image );
  convoluter->SetKernelImage( resampler->GetOutput() );
  convoluter->NormalizeOn();
  SaveImage<PixelType,dimension>(convoluter->GetOutput(),outputVolume,suffix,listFilesStream);
  std::ofstream _nothing;
  if(!outputFilter.empty())
  {
    SaveImage<PixelType,dimension>(resampler->GetOutput(),outputFilter,suffix, _nothing);
  }
}

template<class PixelType, class RealType, unsigned int Dimension>
struct Rotation
{
};

template<class PixelType, class RealType>
struct Rotation<PixelType, RealType, 3>
{
  typedef typename itk::Image< PixelType , 3 > ImageType ;
  typedef typename itk::GaborImageSource<ImageType>::OutputImageType GaborOutputType;
  typedef typename itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
  static void SetRotation(typename GaborOutputType::Pointer gaborOutput,
                          typename ImageType::Pointer image,
                          unsigned int rotations,
                          itk::Point<RealType,3> center,
                          std::string outputVolume,
                          std::string suffix,
                          std::ofstream &listFilesStream,
                          std::string outputFilter
                          )
  {
    typedef itk::ResampleImageFilter<ImageType, ImageType, RealType> ResamplerType;
    typename ResamplerType::Pointer resampler = ResamplerType::New();
    typedef itk::Euler3DTransform<RealType> Euler3DType;
    typename Euler3DType::Pointer euler3D = Euler3DType::New();
    for( unsigned int i = 0 ; i < rotations; i++)
    {
      float theta1 = itk::Math::pi/float(rotations)*float(i) ;
      float theta1_deg = 180.0/float(rotations)*float(i) ;
      for( unsigned int j = 0 ; j < rotations; j++)
      {
        float theta2 = itk::Math::pi/float(rotations)*float(j) ;
        float theta2_deg = 180.0/float(rotations)*float(j) ;
        for( unsigned int k = 0 ; k < rotations; k++)
        {
          float theta3 = itk::Math::pi/float(rotations)*float(k) ;
          float theta3_deg = 180.0/float(rotations)*float(k) ;
          euler3D->SetRotation( theta1, theta2, theta3 );
          std::ostringstream rotation_suffix ;
          rotation_suffix << suffix << "_" << theta1_deg << "_" << theta2_deg << "_" << theta3_deg;
          euler3D->SetCenter( center );
          resampler->SetTransform(euler3D);
          Process<PixelType, RealType, 3>(resampler, gaborOutput, image, outputVolume,rotation_suffix.str(),listFilesStream, outputFilter);
        }
      }
    }
  }
};

template<class PixelType, class RealType>
struct Rotation<PixelType, RealType, 2>
{
  typedef typename itk::Image< PixelType , 2 > ImageType ;
  typedef typename itk::GaborImageSource<ImageType>::OutputImageType GaborOutputType;
  typedef typename itk::ResampleImageFilter<ImageType, ImageType, RealType> ResamplerType;
  static void SetRotation(typename GaborOutputType::Pointer gaborOutput,
                          typename ImageType::Pointer image,
                          unsigned int rotations,
                          itk::Point<RealType,2> center,
                          std::string outputVolume,
                          std::string suffix,
                          std::ofstream &listFilesStream,
                          std::string outputFilter
                          )
  {
    typename ResamplerType::Pointer resampler = ResamplerType::New();
    typedef itk::Euler2DTransform<RealType> Euler2DType;
    typename Euler2DType::Pointer euler2D = Euler2DType::New();
    for( unsigned int i = 0 ; i < rotations; i++)
    {
      float theta = itk::Math::pi/float(rotations)*float(i) ;
      float theta_deg = 180.0/float(rotations)*float(i) ;
      euler2D->SetRotation( theta );
      euler2D->SetCenter( center );
      resampler->SetTransform( euler2D );
      std::ostringstream rotation_suffix;
      rotation_suffix << suffix << "_" << theta_deg;
      Process<PixelType, RealType, 2>(resampler, gaborOutput, image, outputVolume, rotation_suffix.str(),listFilesStream, outputFilter);
    }
  }
};

template< class PixelType, unsigned int dimension >
int GaborImageFilter( int argc , char * argv[] )
{
  PARSE_ARGS ;
  if(filterSize.size() != dimension )
  {
    std::cerr << "Size dimension must match image dimension" << std::endl;
    return EXIT_FAILURE ;
  }
  typedef itk::Image< PixelType , dimension > ImageType ;
  typedef double RealType ;
  typedef itk::ImageFileReader< ImageType > ReaderType ;
  typename ReaderType::Pointer reader = ReaderType::New() ;
  reader->SetFileName( inputVolume ) ;
  reader->Update() ;

  // Create gabor bank
  typedef itk::GaborImageSource<ImageType> GaborSourceType;
  typename GaborSourceType::ArrayType mean;
  typename ImageType::SizeType size ;
  for( size_t i =0 ; i < dimension ; i++)
  {
    size[i] = filterSize[i];
  }

  typename ImageType::PointType origin ;
  origin.Fill(0.0);
  typename ImageType::SpacingType spacing ;
  spacing.Fill(1);
  typename ImageType::DirectionType direction;
  direction.SetIdentity();
  typename itk::ContinuousIndex<RealType,dimension> cIndex;
  typename ImageType::PointType center;
  typename ImageType::Pointer refImage = ImageType::New();
  refImage->SetRegions(size);
  refImage->SetOrigin(origin);
  refImage->SetSpacing(spacing);
  refImage->SetDirection(direction);
  for( unsigned int i = 0; i < dimension; i++ )
  {
    cIndex[i] = size[i]/2.0;
  }
  refImage->TransformContinuousIndexToPhysicalPoint(cIndex,center);
  for( unsigned int i = 0 ; i < dimension ; i++)
  {
    mean[i]=center[i];
  }
  std::cout<<"Mean:"<<mean<<std::endl;
  std::cout<<"origin:"<<origin<<std::endl;
  std::ofstream listFilesStream;
  if( !listFiles.empty())
  {
    listFilesStream.open(listFiles.c_str());
  }
  for( size_t f = 0 ; f < frequencies.size() ; f++)
  {
    for( size_t s=0; s < sigma.size(); s++)
    {
      typename GaborSourceType::ArrayType sigma_c;
      for( unsigned int i = 0; i < dimension; i++ )
      {
        sigma_c[i] = sigma[s];
      }
      std::cout<<"sigma:"<<sigma_c<<std::endl;
      typename GaborSourceType::Pointer gabor = GaborSourceType::New();
      std::cout<< "Creating Filter (sigma="<<sigma[s]<<",freq="<<frequencies[f]<<")"<<std::endl;
      gabor->SetDirection(direction);
      gabor->SetFrequency(frequencies[f]);
      gabor->SetMean(mean);
      gabor->SetOrigin(origin);
      gabor->SetSigma(sigma_c);
      gabor->SetSize(size);
      gabor->SetSpacing(spacing);
      gabor->SetCalculateImaginaryPart(imaginary);
      gabor->Update();
      std::ostringstream suffix ;
      suffix << "_f_" << frequencies[f] << "_s_" << sigma[s] ;
      Rotation<PixelType, RealType, dimension>::SetRotation(gabor->GetOutput(), reader->GetOutput(), rotations, center, outputVolume, suffix.str(), listFilesStream, outputFilter );
    }
  }
  if( !listFiles.empty())
  {
    listFilesStream.close();
  }
  return EXIT_SUCCESS ;
}

template<unsigned int dimension>
int doIt(itk::ImageIOBase::IOComponentType componentType, int argc , char * argv[]
       )
{
  //templated over the input image voxel type
  switch( componentType )
  {
     case itk::ImageIOBase::UCHAR:
        //return GaborImageFilter< unsigned char, dimension>( argc, argv) ;
        //break ;
     case itk::ImageIOBase::CHAR:
//        return GaborImageFilter< char, dimension>( argc, argv) ;
//        break ;
     case itk::ImageIOBase::USHORT:
//        return GaborImageFilter< unsigned short, dimension>( argc, argv) ;
//        break ;
     case itk::ImageIOBase::SHORT:
//        return GaborImageFilter< short, dimension>( argc, argv) ;
//        break ;
     case itk::ImageIOBase::UINT:
//        return GaborImageFilter< unsigned int, dimension>( argc, argv) ;
//        break ;
     case itk::ImageIOBase::INT:
//        return GaborImageFilter< int, dimension>( argc, argv) ;
//        break ;
     case itk::ImageIOBase::ULONG:
//        return GaborImageFilter< unsigned long, dimension>( argc, argv) ;
//        break ;
     case itk::ImageIOBase::LONG:
//        return GaborImageFilter< long, dimension>( argc, argv) ;
//        break ;
     case itk::ImageIOBase::FLOAT:
//        return GaborImageFilter< float, dimension>( argc, argv) ;
//        break ;
     case itk::ImageIOBase::DOUBLE:
        return GaborImageFilter< double, dimension>( argc, argv) ;
        break ;
     case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
     default:
        std::cerr << "unknown component type" << std::endl ;
        break ;
  }
}

int main( int argc , char * argv[] )
{
   PARSE_ARGS;
   itk::ImageIOBase::IOPixelType pixelType ;
   itk::ImageIOBase::IOComponentType componentType ;
   unsigned int dimension=0;
   GetImageType( inputVolume , pixelType , componentType, dimension ) ;
   if( dimension == 2)
   {
     doIt<2>(componentType, argc , argv );
   }
   else if( dimension == 3)
   {
     doIt<3>(componentType, argc , argv );
   }
   else
   {
     std::cerr << "Input dimensions not supported: " << dimension << std::endl;
     return EXIT_FAILURE;
   }

   return EXIT_FAILURE ;
}



