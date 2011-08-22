/** This is a test for the parametric image source fitting classes. It
    generates a 2D Gaussian image with a known sigma, scale, and mean,
    and then uses the ITK registration framework with the adapter class
    itkImageToParametricImageSource to fit a 2D Gaussian to the generated
    image. */
    

#define EXIT_SUCCESS 0
#define EXIT_ERROR 1000


#include <itkAmoebaOptimizer.h>
#include <itkConjugateGradientOptimizer.h>
#include <itkGaussianImageSource.h>
#include <itkImageToParametricImageSourceMetric.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkMeanSquaresImageToImageMetric.h>
#include <itkParametricGaussianImageSource.h>


typedef double
  PixelType;
typedef itk::Image<PixelType, 2>
  GaussianImageType;
typedef itk::GaussianImageSource<GaussianImageType>
  GaussianImageSourceType;
typedef itk::ParametricGaussianImageSource<GaussianImageType>
  ParametricGaussianImageSourceType;
typedef ParametricGaussianImageSourceType::ParametersType
  ParametersType;
typedef itk::LinearInterpolateImageFunction<GaussianImageType, double>
  InterpolatorType;
typedef itk::ImageToParametricImageSourceMetric<GaussianImageType, ParametricGaussianImageSourceType>
  MetricType;
typedef itk::MeanSquaresImageToImageMetric<GaussianImageType, GaussianImageType>
  DelegateMetricType;
typedef itk::AmoebaOptimizer
  AmoebaOptimizerType;
typedef itk::ConjugateGradientOptimizer
  ConjugateGradientOptimizerType;


int main(int argc, char* argv[]) {

  // Parameters common to original and fit Gaussian image
  GaussianImageType::SizeValueType size[2] = {64, 64};
  double spacing[2] = {1.0, 1.0};
  double origin[2]  = {0.0, 0.0};

  // Parameters specific to original Gaussian image
  GaussianImageSourceType::ArrayType trueSigma, trueMean;
  trueSigma[0] = 20.0;
  trueSigma[1] = 30.0;
  trueMean[0]  = 37.0;
  trueMean[1]  = 27.0;
  double trueScale = 4.5;

  ParametersType originalParameters(5);
  originalParameters[0] = trueSigma[0];
  originalParameters[1] = trueSigma[1];
  originalParameters[2] = trueMean[0];
  originalParameters[3] = trueMean[1];
  originalParameters[4] = trueScale;

  // The original image
  GaussianImageSourceType::Pointer originalImage = 
    GaussianImageSourceType::New();
  originalImage->SetSize(size);
  originalImage->SetSpacing(spacing);
  originalImage->SetOrigin(origin);  
  originalImage->SetSigma(trueSigma);
  originalImage->SetMean(trueMean);
  originalImage->SetScale(trueScale);
  originalImage->GetOutput()->Update();

  // The Gaussian image source for fitting
  ParametricGaussianImageSourceType::Pointer fittingImageSource =
    ParametricGaussianImageSourceType::New();
  fittingImageSource->SetSize(size);
  fittingImageSource->SetSpacing(spacing);
  fittingImageSource->SetOrigin(origin);

  // Set up the interpolator
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  // Set up the delegate metric
  DelegateMetricType::Pointer delegate = DelegateMetricType::New();

  // Set up the metric
  MetricType::Pointer metric = MetricType::New();
  metric->SetInterpolator(interpolator);
  metric->SetDelegateMetric(delegate);
  metric->SetFixedImage(originalImage->GetOutput());
  metric->SetFixedImageRegion(originalImage->GetOutput()->GetLargestPossibleRegion());
  metric->SetMovingImageSource(fittingImageSource);
  
  // Set up the optimizer and starting parameters
  ParametersType startingParameters(fittingImageSource->GetNumberOfParameters());
  startingParameters[0] = trueSigma[0] - 5.0;
  startingParameters[1] = trueSigma[1] + 10.0;
  startingParameters[2] = trueMean[0] - 10.0;
  startingParameters[3] = trueMean[1] + 7.0;
  startingParameters[4] = trueScale / 2.0;

  for (unsigned int i = 0; i < fittingImageSource->GetNumberOfParameters(); i++) {
    metric->EnableParameter(i);
  }

  AmoebaOptimizerType::Pointer amoebaOptimizer =
    AmoebaOptimizerType::New();
  amoebaOptimizer->SetCostFunction(metric);
  amoebaOptimizer->SetInitialPosition(startingParameters);
  amoebaOptimizer->SetMaximumNumberOfIterations(1000);
  amoebaOptimizer->StartOptimization();

  // Check the results
  int status = EXIT_SUCCESS;
  ParametersType solution = amoebaOptimizer->GetCurrentPosition();
  for (unsigned int i = 0; i < fittingImageSource->GetNumberOfParameters(); i++) {
    if (fabs(solution[i] - originalParameters[i]) > 1e-2) {
      std::cout << "[FAIL] - AmoebaOptimizer" << std::endl;
      status = EXIT_ERROR;
      break;
    }
  }
  if (status == EXIT_SUCCESS) {
    std::cout << "[PASS] - AmoebaOptimizer" << std::endl;
  }
  std::cout << "Expected parameters: " << originalParameters 
            << ", " << "Recovered parameters: " << solution << std::endl;

  ConjugateGradientOptimizerType::Pointer gradientOptimizer = 
    ConjugateGradientOptimizerType::New();
  gradientOptimizer->SetCostFunction(metric);
  gradientOptimizer->SetInitialPosition(startingParameters);

  // Start the optimization
  gradientOptimizer->StartOptimization();

  // Check the results
  solution = gradientOptimizer->GetCurrentPosition();
  for (unsigned int i = 0; i < fittingImageSource->GetNumberOfParameters(); i++) {
    if (fabs(solution[i] - originalParameters[i]) > 1e-2) {
      std::cout << "[FAIL] - ConjugateGradientOptimizer" << std::endl;
      status = EXIT_ERROR;
      break;
    }
  }
  if (status == EXIT_SUCCESS) {
    std::cout << "[PASS] - ConjugateGradientOptimizer" << std::endl;
  }
  std::cout << "Expected parameters: " << originalParameters 
            << ", " << "Recovered parameters: " << solution << std::endl;

  return status;
}
