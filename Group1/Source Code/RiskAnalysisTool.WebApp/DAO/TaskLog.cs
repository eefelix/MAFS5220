using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.ComponentModel.DataAnnotations.Schema;
using System.Linq;
using System.Text;

namespace RiskAnalysisTool.Worker.DAO
{
    public class TaskLog
    {
        [Key]
        public virtual int ID { get; set; }

        [Index]
        public virtual string State { get; set; }

        [Index]
        public virtual DateTime CreationTime { get; set; }

        public virtual string TaskContent { get; set; }
    }
}
