using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;

namespace RiskAnalysisTool.Instruments
{
    [DataContract]
    public class InterestRateSwap : Swap
    {
        public InterestRateSwap()
        {
            this.FixedRate = 0.02;
            this.FloatingSpread = 0;
            this.Nominal = 1000000;
        }
        [DataMember]
        public double FixedRate { get; set; }

        [DataMember]
        public double FloatingSpread { get; set; }

        [DataMember]
        public double Nominal { get; set; }

        [DataMember]
        public double BlackVolatility { get; set; }

        [DataMember]
        public bool IsPayFloating { get; set; }

    }
}
