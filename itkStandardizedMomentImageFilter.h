/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkStandardizedMomentImageFilter_h
#define itkStandardizedMomentImageFilter_h

#include "itkBoxImageFilter.h"
#include "itkImage.h"
#include "itkNumericTraits.h"
#include "itkConstNeighborhoodIterator.h"

namespace itk
{
/** \class StandardizedMomentImageFilter
 * \brief Calculate the local noise in an image.
 *
 * Computes an image where a given pixel is the standard deviation of
 * the pixels in a neighborhood about the corresponding input pixel.
 * This serves as an estimate of the local noise (or texture) in an
 * image. Currently, this noise estimate assume a piecewise constant
 * image.  This filter should be extended to fitting a (hyper) plane
 * to the neighborhood and calculating the standard deviation of the
 * residuals to this (hyper) plane.
 *
 * \sa Image
 * \sa Neighborhood
 * \sa NeighborhoodOperator
 * \sa NeighborhoodIterator
 *
 * \ingroup IntensityImageFilters
 * \ingroup ITKImageFilterBase
 *
 * \wiki
 * \wikiexample{Statistics/StandardizedMomentImageFilter,Compute the local noise in an image}
 * \endwiki
 */
template< typename TInputImage >
class StandardizedMomentImageFilter:
  public BoxImageFilter< TInputImage,
        itk::Image<typename NumericTraits< TInputImage >::RealType,TInputImage::ImageDimension > >
{
public:
  /** Extract dimension from input and output image. */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
//NumericTraits< typename TInputImage::PixelType>::RealType
  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage  InputImageType;

  /** Standard class typedefs. */
  typedef StandardizedMomentImageFilter                                  Self;
  typedef SmartPointer< Self >                              Pointer;
  typedef SmartPointer< const Self >                        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

//  /** Run-time type information (and related methods). */
//  itkTypeMacro(StandardizedMomentImageFilter, BoxImageFilter);

//  /** Image typedef support. */
  typedef typename InputImageType::PixelType                 InputPixelType;
  typedef typename NumericTraits< InputPixelType >::RealType InputRealType;
  typedef typename NumericTraits< InputPixelType >::RealType OutputPixelType;
  typedef typename itk::Image<OutputPixelType,TInputImage::ImageDimension >    OutputImageType;
  typedef BoxImageFilter< InputImageType, OutputImageType > Superclass;

  typedef typename InputImageType::RegionType  InputImageRegionType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef typename InputImageType::SizeType InputSizeType;

#ifdef ITK_USE_CONCEPT_CHECKING
  // Begin concept checking
  itkConceptMacro( InputHasNumericTraitsCheck,
                   ( Concept::HasNumericTraits< InputPixelType > ) );
   //End concept checking
#endif
  itkSetMacro(Bias,bool);
  itkGetMacro(Bias,bool);
protected:
  StandardizedMomentImageFilter();
  virtual ~StandardizedMomentImageFilter() {}

  /** StandardizedMomentImageFilter can be implemented as a multithreaded filter.
   * Therefore, this implementation provides a ThreadedGenerateData()
   * routine which is called for each processing thread. The output
   * image data is allocated automatically by the superclass prior to
   * calling ThreadedGenerateData().  ThreadedGenerateData can only
   * write to the portion of the output image specified by the
   * parameter "outputRegionForThread"
   *
   * \sa BoxImageFilter::ThreadedGenerateData(),
   *     BoxImageFilter::GenerateData() */
  InputRealType GetStandardizedMoment(ConstNeighborhoodIterator< InputImageType > bit, unsigned int neighborhoodSize, unsigned int moment );
private:
  StandardizedMomentImageFilter(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;
  bool m_Bias;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkStandardizedMomentImageFilter.hxx"
#endif

#endif
