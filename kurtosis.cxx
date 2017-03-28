#include "kurtosisCLP.h"
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkKurtosisImageFilter.h"
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
int KurtosisImageFilter( std::string inputVolume , std::string outputVolume, int radius, bool bias, bool excess )
{

  typedef itk::Image< PixelType , dimension > ImageType ;
  typedef itk::ImageFileReader< ImageType > ReaderType ;
  typename ReaderType::Pointer reader = ReaderType::New() ;
  reader->SetFileName( inputVolume ) ;
  reader->Update() ;
  typedef itk::KurtosisImageFilter<ImageType> KurtosisFilterType ;
  typename KurtosisFilterType::Pointer kurt= KurtosisFilterType::New();
  kurt->SetInput(reader->GetOutput());
  kurt->SetBias(bias);
  kurt->SetRadius(radius);
  kurt->SetExcess(excess);
  typedef itk::ImageFileWriter< typename KurtosisFilterType::OutputImageType > WriterType ;
  typename WriterType::Pointer writer = WriterType::New() ;
  writer->SetFileName( outputVolume ) ;
  writer->SetInput( kurt->GetOutput() ) ;
  writer->SetUseCompression( 1 ) ;
  writer->Update() ;
  return EXIT_SUCCESS ;
}

template<unsigned int dimension>
int doIt(itk::ImageIOBase::IOComponentType componentType,
       std::string inputVolume ,
       std::string outputVolume,
       int radius,
       bool bias,
       bool excess
       )
{
  //templated over the input image voxel type
  switch( componentType )
  {
     case itk::ImageIOBase::UCHAR:
        return KurtosisImageFilter< unsigned char, dimension>( inputVolume , outputVolume, radius, bias, excess) ;
        break ;
     case itk::ImageIOBase::CHAR:
        return KurtosisImageFilter< char, dimension>( inputVolume , outputVolume, radius, bias, excess) ;
        break ;
     case itk::ImageIOBase::USHORT:
        return KurtosisImageFilter< unsigned short, dimension>( inputVolume , outputVolume, radius, bias, excess) ;
        break ;
     case itk::ImageIOBase::SHORT:
        return KurtosisImageFilter< short, dimension>( inputVolume , outputVolume, radius, bias, excess) ;
        break ;
     case itk::ImageIOBase::UINT:
        return KurtosisImageFilter< unsigned int, dimension>( inputVolume , outputVolume, radius, bias, excess) ;
        break ;
     case itk::ImageIOBase::INT:
        return KurtosisImageFilter< int, dimension>( inputVolume , outputVolume, radius, bias, excess) ;
        break ;
     case itk::ImageIOBase::ULONG:
        return KurtosisImageFilter< unsigned long, dimension>( inputVolume , outputVolume, radius, bias, excess) ;
        break ;
     case itk::ImageIOBase::LONG:
        return KurtosisImageFilter< long, dimension>( inputVolume , outputVolume, radius, bias, excess) ;
        break ;
     case itk::ImageIOBase::FLOAT:
        return KurtosisImageFilter< float, dimension>( inputVolume , outputVolume, radius, bias, excess) ;
        break ;
     case itk::ImageIOBase::DOUBLE:
        return KurtosisImageFilter< double, dimension>( inputVolume , outputVolume, radius, bias, excess) ;
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
     doIt<2>(componentType, inputVolume , outputVolume, radius, bias, excess );
   }
   else if( dimension == 3)
   {
     doIt<3>(componentType, inputVolume , outputVolume, radius, bias, excess );
   }
   else
   {
     std::cerr << "Input dimensions: " << dimension << std::endl;
     return EXIT_FAILURE;
   }

   return EXIT_FAILURE ;
}



