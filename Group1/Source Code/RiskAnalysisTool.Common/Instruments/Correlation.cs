using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;

namespace RiskAnalysisTool.Instruments
{
    [DataContract]
    public class Correlation : Instrument
    {
        public Correlation()
        {
            this.Symbol1 = " ";
            this.Symbol2 = " ";
            this.Value = 0;
        }
        [DataMember]
        public string Symbol1 { get; set; }

        [DataMember]
        public string Symbol2 { get; set; }

        [DataMember]
        public double Value { get; set; }
    }
}
