using System;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Threading.Tasks;
using RiskAnalysisTool.Instruments;
using RiskAnalysisTool.MobileApp.Infrastructure;
using RiskAnalysisTool.Tasks;
using XLabs;
using XLabs.Forms.Mvvm;

namespace RiskAnalysisTool.MobileApp.ViewModels
{
    public class TaskListViewModel : ViewModel
    {
        //private readonly AsyncCommand _addTaskCommand;

        //public AsyncCommand AddTaskCommand
        //{
        //    get { return _addTaskCommand; }
        //}

        //private readonly AsyncCommand<ComputationTask> _editTaskCommand;

        //public AsyncCommand<ComputationTask> EditTaskCommand
        //{
        //    get { return _editTaskCommand; }
        //}

        //private readonly AsyncCommand<ComputationTask> _deleteTaskCommand;

        //public AsyncCommand<ComputationTask> DeleteTaskCommand
        //{
        //    get { return _deleteTaskCommand; }
        //}

        //private ObservableCollection<ComputationTask> _tasks;

        //public TaskListViewModel()
        //{
        //    _tasks = new ObservableCollection<ComputationTask>();
        //    _addTaskCommand = new AsyncCommand(AddTask);
        //    _editTaskCommand = new AsyncCommand<ComputationTask>(EditTask);
        //    _deleteTaskCommand = new AsyncCommand<ComputationTask>(DeleteTask);
        //}

        //public ObservableCollection<ComputationTask> Tasks
        //{
        //    get { return _tasks; }
        //    set { SetProperty(ref _tasks, value); }
        //}

        //public async Task AddTask()
        //{
        //    ComputationTask task = await this.InvokeAsync<TaskDetailViewModel, ComputationTask>();
        //    ServiceHelper.CallService<ComputationTask, object>("api/Task", "POST", task);
        //}
        //public async Task DeleteTask(ComputationTask task)
        //{
        //    await ServiceHelper.CallService<ComputationTask, object>("api/Task/" + task.ID, "DELETE", task);
        //}

        //public async Task EditTask(ComputationTask task)
        //{
        //    ComputationTask newTask = await this.InvokeAsync<TaskDetailViewModel, ComputationTask>(
        //        (vm, v) =>
        //        {
        //            vm.LoadModel(task);
        //        });

        //    await ServiceHelper.CallService<ComputationTask, object>("api/Task", "POST", newTask);
        //}

    }
}