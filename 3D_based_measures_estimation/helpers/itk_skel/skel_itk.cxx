#include "itkImage.h"
#include "itkMetaImageIO.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkBinaryThinningImageFilter3D.h"
#include "itkMedialThicknessImageFilter3D.h"
#include <iostream>
#include <string>

int main(int argc, char * argv[])
{
	
	// set filenames
	std::string inDir = argv[1];
	std::string coralName = argv[2];
	std::string filenameIn = inDir + "/" + coralName + ".mhd";
	std::string filenameOut = inDir + "/" + coralName + "_skel.mhd";
	std::cout << filenameOut << std::endl;
	// Read a .mhd/.zraw file into an image
	typedef unsigned char PixelType;
	typedef itk::Image< PixelType, 3 > ImageType;
	typedef float ThickType;
	typedef itk::Image<ThickType, 3 > ThickImageType;
	typedef itk::ImageFileReader < ImageType > VolumeReaderType;
	VolumeReaderType::Pointer reader = VolumeReaderType::New();
	
	reader->SetFileName(filenameIn);
	reader->Update();
	ImageType::Pointer image;
	image = reader->GetOutput();

	ImageType::RegionType region = image->GetLargestPossibleRegion();
	// use binary thinning to create voxel skeleton
	using FilterType = itk::BinaryThinningImageFilter3D<ImageType, ImageType>;
	FilterType::Pointer thinFilter = FilterType::New();
	std::cout << "Begin skel" << std::endl;
	thinFilter->SetInput(image);
	thinFilter->Update();
	ImageType::Pointer voxelImage = thinFilter->GetOutput();

	// get medial thickness
	using FilterTypeThick = itk::MedialThicknessImageFilter3D<ImageType, ThickImageType>;
	FilterTypeThick::Pointer thickFilter = FilterTypeThick::New();
	std::cout << "begin thickness" << std::endl;
	thickFilter->SetInput(image);
	thickFilter->Update();
	ThickImageType::Pointer thickImage = thickFilter->GetOutput();
	std::cout << thickImage->GetLargestPossibleRegion().GetNumberOfPixels() << std::endl;
	
	using IteratorThickType = itk::ImageRegionIterator<ThickImageType>;
	IteratorThickType thickIter(thickImage, thickImage->GetLargestPossibleRegion());
	
	using IteratorVoxType = itk::ImageRegionIterator<ImageType>;
	IteratorVoxType skelIter(voxelImage, voxelImage->GetLargestPossibleRegion());
	std::cout << "begin checking stuff" << std::endl;
	thickIter.GoToBegin();
	skelIter.GoToBegin();

	while(!thickIter.IsAtEnd()){
		
		if (!skelIter.Get()) {
			thickIter.Set(-1);
		}
		++skelIter;
		++thickIter;
	}

	
	using ThickWriterType = itk::ImageFileWriter<ThickImageType>;
	ThickWriterType::Pointer writerThick = ThickWriterType::New();

	writerThick->SetFileName(filenameOut);
	writerThick->SetInput(thickImage);
	writerThick->Update();
	writerThick->Write();

}
