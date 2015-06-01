//using System;
//using System.Collections.Generic;
//using System.Linq;
//using System.Net;
//using System.Net.Http;
//using System.Web.Http;
//using Newtonsoft.Json;
//using RiskAnalysisTool.Tasks;
//using RiskAnalysisTool.Worker;
//using RiskAnalysisTool.Worker.DAO;

//namespace RiskAnalysisTool.WebApp.Controllers
//{
//    public class TaskController : ApiController
//    {
//        // GET: api/Task
//        public IEnumerable<ComputationTask> Get()
//        {
//            using (TaskRepository repository = new TaskRepository())
//            {
//                var q = from t in repository.Task
//                        orderby t.CreationTime
//                        select t;

//                q = q.Select(e =>
//                {
//                    ComputationTask task = JsonConvert.DeserializeObject<ComputationTask>(e.TaskContent);
//                    task.ID = e.ID;
//                });

//                return q.ToArray();
//            }
//        }

//        // GET: api/Task/5
//        public ComputationTask Get(int id)
//        {
//            using (TaskRepository repository = new TaskRepository())
//            {
//                var q = from t in repository.Task
//                        where t.ID == id
//                        orderby t.CreationTime
//                        select t;
               
//                q = q.Select(e =>
//                {
//                    ComputationTask task = JsonConvert.DeserializeObject<ComputationTask>(e.TaskContent);
//                    task.ID = e.ID;
//                });

//                return q.FirstOrDefault();
//            }
//        }

//        // POST: api/Task
//        public void Post([FromBody]ComputationTask value)
//        {
//            using (TaskRepository repository = new TaskRepository())
//            {
//                TaskLog taskLog = repository.Task.Create();
//                taskLog.CreationTime = DateTime.UtcNow;
//                taskLog.State = TaskState.Pending.ToString();
//                taskLog.TaskContent = JsonConvert.SerializeObject(value);
//                repository.SaveChanges();
//            }
//        }

//        // PUT: api/Task/5
//        public void Put(int id, [FromBody]ComputationTask value)
//        {
//            using (TaskRepository repository = new TaskRepository())
//            {
//                var q = from t in repository.Task
//                        where t.ID == TaskState.Running.ToString()
//                            && t.State != task
//                        orderby t.CreationTime
//                        select new TaskLog { ID = t.ID, State = t.State, CreationTime = t.CreationTime };

//                TaskLog tasklog = q.FirstOrDefault();
//                if (tasklog == null)
//                {
//                    this.NotFound();
//                }

//                tasklog.State = TaskState.Pending.ToString();
//                repository.SaveChanges();
//            }
//        }

//        // DELETE: api/Task/5
//        public void Delete(int id)
//        {
//            using (TaskRepository repository = new TaskRepository())
//            {
//                var q = from t in repository.Task
//                        where t.ID == TaskState.Running.ToString()
//                            && t.State != task
//                        orderby t.CreationTime
//                        select new TaskLog { ID = t.ID, State = t.State, CreationTime = t.CreationTime };

//                TaskLog tasklog = q.FirstOrDefault();
//                if (tasklog == null)
//                {
//                    return;
//                }

//                repository.Task.Remove(tasklog);
//                repository.SaveChanges();
//            }
//        }
//    }
//}
