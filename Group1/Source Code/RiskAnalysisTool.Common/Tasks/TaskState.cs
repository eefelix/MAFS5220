using System;
using System.Collections.Generic;
using System.Linq;

namespace RiskAnalysisTool.Tasks
{
    public enum TaskState
    {
        Pending = 0,
        Running,
        Complete,
        Failed
    }
}
