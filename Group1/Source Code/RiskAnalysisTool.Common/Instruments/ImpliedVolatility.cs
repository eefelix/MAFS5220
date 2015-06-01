using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;

namespace RiskAnalysisTool.Instruments
{
    [DataContract]
    public class ImpliedVolatility : Instrument
    {
        [DataMember]
        public double Value { get; set; }
    }
}
