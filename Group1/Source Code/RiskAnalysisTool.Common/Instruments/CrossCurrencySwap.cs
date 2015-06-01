using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;

namespace RiskAnalysisTool.Instruments
{
    [DataContract]
    public class CrossCurrencySwap : Swap
    {
        public CrossCurrencySwap()
        {
            this.PayNominal = 1000000;
            this.PayRate = 0.01;
            this.ReceiveNominal = 1000000;
            this.ReceiveRate = 0.01;
            this.ExchangeVolatility = 0.2;
        }

        [DataMember]
        public double PayNominal { get; set; }

        [DataMember]
        public double PayRate { get; set; }

        [DataMember]
        public double ReceiveNominal { get; set; }

        [DataMember]
        public double ReceiveRate { get; set; }

        [DataMember]
        public double ExchangeVolatility { get; set; }



    }
}
