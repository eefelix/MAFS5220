#include <pch.h>
#include <boost/thread/thread.hpp>
#include <boost/thread.hpp>
#include <ql/types.hpp>


namespace QuantLib {

    Integer sessionId() {
#if defined _WINDLL || defined _WINDOWS

#include <Windows.h>
		
		return (Integer)::GetCurrentThreadId();
#else
        reutrn 0;
#endif
    }
}