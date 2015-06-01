#pragma once

#include <ql/methods/montecarlo/sample.hpp>
#include <vector>
#include <iostream>
namespace QuantLib
{
	template <typename __Sample_Type>
	class FakeUniformRng
	{
	private:
		std::vector<__Sample_Type> sequence_;
	public:
		typedef QuantLib::Sample<__Sample_Type> sample_type;

		explicit FakeUniformRng(const std::vector<__Sample_Type> &sequence)
		{
			sequence_ = sequence;
			cur_ = 0;
		}

		sample_type next() const
		{
			cur_ %= sequence_.size();
			return sample_type(sequence_[cur_], 1.0);
		}

		QuantLib::Real nextReal() const;
		unsigned long nextInt32() const;
	private:
		mutable size_t cur_;
	};

	template <>
	QuantLib::Real FakeUniformRng<QuantLib::Real>::nextReal() const
	{
		return next().value;
	}

	template <>
	unsigned long FakeUniformRng<unsigned long>::nextInt32() const
	{
		return next().value;
	}
}