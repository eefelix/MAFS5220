using System;
using System.Linq;
using System.Runtime.Serialization;

namespace RiskAnalysisTool.Results
{
    [DataContract]
    public class BvaResult : ComputationResult
    {
        public double Bva
        {
            get { return Dva - Cva; }
        }

        [DataMember]
        public double Cva { get; set; }

        [DataMember]
        public double Dva { get; set; }

        public override string ToString()
        {
            return string.Format("CVA: {0}, DVA: {1}, BVA = {2}", Cva, Dva, Bva);
        }
    }
}
