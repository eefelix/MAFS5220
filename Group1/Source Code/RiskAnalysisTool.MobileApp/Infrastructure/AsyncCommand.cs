using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Threading.Tasks;
using System.Windows.Input;
using XLabs;

namespace RiskAnalysisTool.MobileApp.Infrastructure
{
    public class AsyncCommand<T> : RelayCommand<T>, ICommand
    {
        public AsyncCommand(Func<T, Task> task)
            : base((arg) => RunTask(task, (T)arg))
        {
        }

        public AsyncCommand(Func<T, Task> task, Predicate<T> canExecute)
            : base((arg) => RunTask(task, (T)arg), canExecute)
        {
        }

        private static async Task RunTask(Func<T, Task> task, T argument)
        {
            try
            {
                await task(argument);
            }
            catch (Exception ex)
            {
                Debug.WriteLine(ex);
            }
            finally
            {
            }
        }
    }

    public class AsyncCommand : RelayCommand, ICommand
    {
        public AsyncCommand(Func<Task> task)
            : base(() => RunTask(task))
        {
        }

        public AsyncCommand(Func<Task> task, Func<bool> canExecute)
            : base(() => RunTask(task), canExecute)
        {
        }

        private static async Task RunTask(Func<Task> task)
        {
            try
            {
                await task();
            }
            catch (Exception ex)
            {
                Debug.WriteLine(ex);
            }
            finally
            {
            }
        }
    }
}
