using System;
using System.Linq;
using System.Runtime.Serialization;

namespace RiskAnalysisTool.Results
{
    [DataContract]
    public class TvaResult : ComputationResult
    {
        [DataMember]
        public double Cva { get; set; }

        [DataMember]
        public double Dva { get; set; }

        [DataMember]
        public double Fva { get; set; }

        public override string ToString()
        {
            return string.Format("FVA: ", Fva);
        }
    }
}
