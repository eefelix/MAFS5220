using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using RiskAnalysisTool.Time;

namespace RiskAnalysisTool.Instruments
{
    [DataContract]
    public class EquitySwap : Swap
    {
        public EquitySwap()
        {
            this.Amount = 10000;
            this.StartPrice = 50;
            this.SpotPrice = 50;
            this.FixedRate = 0.02;
            this.DividendYield = 0.01;
            this.Volatility = 0.2;
        }
        public int Amount { get; set; }

        public double DividendYield { get; set; }

        public double FixedRate { get; set; }

        public double SpotPrice { get; set; }

        public double StartPrice { get; set; }

        public double Volatility { get; set; }
    }
}
