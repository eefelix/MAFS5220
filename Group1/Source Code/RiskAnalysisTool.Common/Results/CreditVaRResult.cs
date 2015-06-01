using System;
using System.Linq;
using System.Runtime.Serialization;

namespace RiskAnalysisTool.Results
{
    [DataContract]
    public class CreditVaRResult : ComputationResult
    {
        [DataMember]
        public double CreditES { get; set; }

        [DataMember]
        public double CreditVaR { get; set; }

        public override string ToString()
        {
            return string.Format("Credit VaR: ", CreditVaR);
        }

    }
}
