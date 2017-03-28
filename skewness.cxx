#include "skewnessCLP.h"
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkSkewnessImageFilter.h"
#include <itkImageRegionIteratorWithIndex.h>

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

template< class PixelType, unsigned int dimension >
int SkewnessImageFilter( std::string inputVolume , std::string outputVolume, int radius, bool bias )
{

  typedef itk::Image< PixelType , dimension > ImageType ;
  typedef itk::ImageFileReader< ImageType > ReaderType ;
  typename ReaderType::Pointer reader = ReaderType::New() ;
  reader->SetFileName( inputVolume ) ;
  reader->Update() ;
  typedef itk::SkewnessImageFilter<ImageType> SkewnessFilterType ;
  typename SkewnessFilterType::Pointer skewness= SkewnessFilterType::New();
  skewness->SetInput(reader->GetOutput());
  skewness->SetBias(bias);
  skewness->SetRadius(radius);
  skewness->Update();
  typedef itk::ImageFileWriter< typename SkewnessFilterType::OutputImageType > WriterType ;
  typename WriterType::Pointer writer = WriterType::New() ;
  writer->SetFileName( outputVolume ) ;
  writer->SetInput( skewness->GetOutput() ) ;
  writer->SetUseCompression( 1 ) ;
  writer->Update() ;
  return EXIT_SUCCESS ;
}

template<unsigned int dimension>
int doIt(itk::ImageIOBase::IOComponentType componentType,
       std::string inputVolume ,
       std::string outputVolume,
       int radius,
       bool bias
       )
{
  //templated over the input image voxel type
  switch( componentType )
  {
     case itk::ImageIOBase::UCHAR:
        return SkewnessImageFilter< unsigned char, dimension>( inputVolume , outputVolume, radius, bias ) ;
        break ;
     case itk::ImageIOBase::CHAR:
        return SkewnessImageFilter< char, dimension>( inputVolume , outputVolume, radius, bias ) ;
        break ;
     case itk::ImageIOBase::USHORT:
        return SkewnessImageFilter< unsigned short, dimension>( inputVolume , outputVolume, radius, bias ) ;
        break ;
     case itk::ImageIOBase::SHORT:
        return SkewnessImageFilter< short, dimension>( inputVolume , outputVolume, radius, bias ) ;
        break ;
     case itk::ImageIOBase::UINT:
        return SkewnessImageFilter< unsigned int, dimension>( inputVolume , outputVolume, radius, bias ) ;
        break ;
     case itk::ImageIOBase::INT:
        return SkewnessImageFilter< int, dimension>( inputVolume , outputVolume, radius, bias ) ;
        break ;
     case itk::ImageIOBase::ULONG:
        return SkewnessImageFilter< unsigned long, dimension>( inputVolume , outputVolume, radius, bias ) ;
        break ;
     case itk::ImageIOBase::LONG:
        return SkewnessImageFilter< long, dimension>( inputVolume , outputVolume, radius, bias ) ;
        break ;
     case itk::ImageIOBase::FLOAT:
        return SkewnessImageFilter< float, dimension>( inputVolume , outputVolume, radius, bias ) ;
        break ;
     case itk::ImageIOBase::DOUBLE:
        return SkewnessImageFilter< double, dimension>( inputVolume , outputVolume, radius, bias ) ;
        break ;
     case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
     default:
        std::cerr << "unknown component type" << std::endl ;
        break ;
  }
}

int main( int argc , char * argv[] )
{
   PARSE_ARGS ;
   itk::ImageIOBase::IOPixelType pixelType ;
   itk::ImageIOBase::IOComponentType componentType ;
   unsigned int dimension=0;
   GetImageType( inputVolume , pixelType , componentType, dimension ) ;
   if( dimension == 2)
   {
     doIt<2>(componentType, inputVolume , outputVolume, radius, bias );
   }
   else if( dimension == 3)
   {
     doIt<3>(componentType, inputVolume , outputVolume, radius, bias );
   }
   else
   {
     std::cerr << "Input dimensions: " << dimension << std::endl;
     return EXIT_FAILURE;
   }

   return EXIT_FAILURE ;
}



