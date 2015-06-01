using System;
using System.Linq;
using System.Runtime.Serialization;

namespace RiskAnalysisTool.Results
{
    [DataContract]
    public class UCvaResult : ComputationResult
    {
        [DataMember]
        public double UCva { get; set; }

        public override string ToString()
        {
            return string.Format("UCVA: ", UCva);
        }
    }
}
