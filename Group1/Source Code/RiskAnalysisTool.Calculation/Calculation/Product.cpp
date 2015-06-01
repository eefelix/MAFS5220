#include <pch.h>
#include "Product.hpp"

std::map< std::string, boost::shared_ptr<QuantLib::StochasticProcess1D> >
	Product::symToProcess = std::map< std::string, boost::shared_ptr<QuantLib::StochasticProcess1D> >();