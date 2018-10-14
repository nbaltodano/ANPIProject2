#ifndef ANPI_RESISTORGRID_HPP
#define ANSPI_RESISTORGRID_HPP

#include <iostream>
#include <cstdlib>
#include <string>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <AnpiConfig.hpp>
#include <Matrix.hpp>
#include <Exception.hpp>

namespace anpi
{
	//Pack a pair of indices od the nodes of a resistor
	struct indexPair
	{
		//Row of the first node
		std::size_t row1;
		//Column of the first node
		std::size_t col1;
		//Row of the second node
		std::size_t row2;
		//Column of the second node
		std::size_t col2;
	};

	struct Nodo
    {
        std::size_t  _row;
        std::size_t  _col;

        Nodo (size_t row, size_t col)
        {
            this ->_row = row;
            this ->_col = col;
        }
    };

	class ResistorGrid
	{
		private:
		//Matrix of the current equation system
		Matrix<float> A;
		//Vector of the current equation system
		std::vector<float> b;
		//Raw map data
		Matrix<float> rawMap;
		const size_t lowResist = 1;
		const size_t highResist = 1000000;
		cv::Mat<float> rawMap_CV;

		public:
			Matrix<float> _X, _Y;
		

        ResistorGrid()
        {
            const size_t  _totalSize = 2 * rawMap.rows() * rawMap.cols() - (rawMap.rows() + rawMap.cols()) +1;
            A.allocate(_totalSize, _totalSize);
            A.fill(0);
            b.resize(_totalSize);
        }
		/∗∗
		∗ Construct the grid from the given file
		∗ @return true if successful or false otherwise
		∗/
		bool build(const std::string filename);
		/∗∗
		∗ Compute the internal data to navigate between the given nodes
		∗/
		bool navigate(const index Pair& nodes);

		std::size_t minor(const std::size_t num1, const std::size_t num2)
		{
			if ((std::abs(int(num1))) >= (std::abs(int(num2))))
				return num2;
			else
				return num1;
		}

		std::size_t ResistorGrid::nodesToIndex (cost std::size_t row1,
												cons std::size_t col1,
												cons std::size_t row2,
												cons std::size_t col2)
        {
        	
        	if ((((std::abs(int(row2 - row1)) == 1) && (std:abs(int(col2 - col1))) ==0)) || ((std::abs(int(row2 - row1)) ==0) && (std::abs(int(col2 - col1)) == 1)))
        	{
        		size_t resistRow, resistCol, finalIndex;
        		if  ((std::abs(int(row2 - row1)) ==0) && (std::abs(int(col2 - col1)) == 1))
        		{
        			resistRow = 2*row1;
        			resistCol = minor(col1, col2);
        		}
        		else
        		{
        			resistRow = row1 + row2;
        			resistCol = col1;
        		}
        		size_t finalIndex, n;
        		n = rawMap.cols();
        		if (col1 == col2)  //Vertical
        			finalIndex = resistRow * n + resistRow + 1 - (resistRow - resistCol) + (resistRow - 3)/2;
        		else 
        			finalIndex = resistRow/2 * 2*n + resistCol + resistRow/2;
        	}
        	else
        		throw anpi::Exception ("Nodes are not next to each other");
		}

		indexPair ResistorGrid::indexToNodes (const std::size_t _index)
		{
			if(rawMap.cols() * (rawMap.rows() +1) + (rawMap.cols() + 1) *rawMap.rows() - 1)
			{
				int(size_t) nC = rawMap.cols();
				int(size_t) nR = rawMap.rows();
				size_t row_node, col_node
				indexPair result;
				if(_index == 0)
				{
					result.row1 = 0;
					result.col1 = 0;
					result.row2 = 0;
					result.col2 = 0;
				}
				else
				{
					for(size_t i = 0, i<= rawMap.cols()+1, i++)
						for (size_t j = 0, rawMap.rows()+1, j++)
					{
							row_node = nodesToIndex(i, j, i, j+1);
							col_node = nodesToIndex(i, j, i+1, j);

						if (row_node >= _index)
						{
							result.row1 = i;
							result.col1 = j;
							result.row2 = i;
							result.col2 = j+1;
							break;
						}
						else if (col_node >= _index)
						{
							result.row1 = i;
							result.col1 = j;
							result.row2 = i+1;
							result.col2 = j;
							break;
						}
					}
					return result;		
				}
			}

			else
			{
				throw anpi::Exception("Index out of range");
			}

		}

		size_t assign_rest(const std::size_t _idex)
		{
			const indexPair nodeToFind = indexToNodes(_idex);
			size_t resistVal;
			if(rawMap(nodeToFind.row1, nodeToFind.col1)==0 || rawMap(nodeToFind.row2, nodeToFind.col2) == 0)
				resistVal = highResist;
			else
				resistVal = lowResist;
			return resistVal;
		}

		bool build(const std::string &filename)
		{
			// Build the name of the image in the data path
      		std::string mapPath = std::string( ANPI_DATA_PATH ) + "/" + filename;

      		// Read the image using the OpenCV
      		cv::imread(mapPath.c_str(),
                 		CV_LOAD_IMAGE_GRAYSCALE).convertTo(rawMap_CV,CV_32FC1);
      		rawMap_CV /= 255.0f; // normalize image range to 0 .. 255

      		if(rawMap_CV.cols == 0 || rawMap_CV.rows == 0 || rawMap_CV.data == NULL) 
      		{
        		throw anpi::Exception("Problem creating the map");
     		}

      		// Convert the OpenCV matrix into an anpi matrix
      		// We have to use the std::allocator to avoid an exact stride
      		anpi::Matrix<float,std::allocator<float> > amapTmp(static_cast<const size_t>(rawMap_CV.rows),
            		                                             static_cast<const size_t>(rawMap_CV.cols),
               			                                          rawMap_CV.ptr<float>());
      		// And transform it to a SIMD-enabled matrix
      		anpi::Matrix<float> amap(amapTmp);
      		rawMap = amap;

      		// Initialize the other data
      		const size_t totalVariables = 2*rawMap.rows()*rawMap.cols()-(rawMap.rows()+rawMap.cols());
      		A_.allocate(totalVariables,totalVariables);
      		A_.fill(static_cast<float>(0));
      		_X.allocate(rawMap.rows(),rawMap.cols());
      		_X.fill(static_cast<float >(0));
      		_Y.allocate(rawMap.rows(),rawMap.cols());
      		_Y.fill(static_cast<float >(0));
      		return true;
		}



	};
}
