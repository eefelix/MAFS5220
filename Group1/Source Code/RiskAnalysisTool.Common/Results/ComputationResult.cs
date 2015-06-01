using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;

namespace RiskAnalysisTool.Results
{
    [DataContract]
    [KnownType(typeof(CreditVaRResult))]
    [KnownType(typeof(UCvaResult))]
    [KnownType(typeof(BvaResult))]
    [KnownType(typeof(TvaResult))]
    public abstract class ComputationResult
    {
    }
}
