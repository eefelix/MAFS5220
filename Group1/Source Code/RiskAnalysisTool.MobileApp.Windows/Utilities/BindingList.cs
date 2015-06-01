using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace System.ComponentModel
{
    public class BindingList<T> : ObservableCollection<T>
    {
        public BindingList()
            : base()
        {
        }

        public BindingList(IList<T> list)
            : base(list)
        {
        }

        public void ResetItem(int index)
        {
            this.SetItem(index, this[index]);
        }
    }
}
