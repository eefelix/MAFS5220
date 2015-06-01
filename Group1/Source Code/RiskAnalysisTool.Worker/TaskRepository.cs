using System;
using System.Collections.Generic;
using System.Data.Entity;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using RiskAnalysisTool.Tasks;
using RiskAnalysisTool.Worker.DAO;

namespace RiskAnalysisTool.Worker
{
    class TaskRepository : DbContext
    {
        public TaskRepository()
            : base()
        {
        }

        protected override void OnModelCreating(DbModelBuilder modelBuilder)
        {
            base.OnModelCreating(modelBuilder);
            modelBuilder.Entity<TaskLog>().ToTable("Tasks");
        }

        public DbSet<TaskLog> Task
        {
            get { return this.Set<TaskLog>(); }
        }
    }
}
