#pragma once

#include <memory>
#include <boost/shared_ptr.hpp>

namespace RiskAnalysisTool {
namespace Calculation {
    template<typename T>
    inline boost::shared_ptr<T> make_boost_ptr(std::shared_ptr<T>& ptr)
    {
        return boost::shared_ptr<T>(ptr.get(), [ptr](T*) mutable
        {
            ptr.reset();
        });
    }

    template<typename T>
    inline std::shared_ptr<T> make_std_ptr(boost::shared_ptr<T>& ptr)
    {
        return std::shared_ptr<T>(ptr.get(), [ptr](T*) mutable
        {
            ptr.reset();
        });
    }
}
}